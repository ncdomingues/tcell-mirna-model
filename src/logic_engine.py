# -*- coding: utf-8 -*-
"""
A small multi-valued logical (Boolean-network-style) simulation engine,
built to read the exact rule notation used in Table 5 of the thesis:
    &  AND       |  OR       !  NOT       NAME:k  "NAME has reached value >= k"

This intentionally avoids external logical-modelling packages (GINsim, R's
BoolNet/MaBoSS) -- it is a from-scratch, inspectable reimplementation of the
same threshold multi-valued semantics, in plain Python/numpy.
"""
import re
import random
from collections import defaultdict

TOKEN_RE = re.compile(r"""
    (?P<WS>\s+)
  | (?P<LPAREN>\()
  | (?P<RPAREN>\))
  | (?P<AND>&)
  | (?P<OR>\|)
  | (?P<NOT>!)
  | (?P<IDENT>[A-Za-z][A-Za-z0-9\-]*(?::\d+)?)
""", re.VERBOSE)


def tokenize(formula):
    pos = 0
    tokens = []
    while pos < len(formula):
        m = TOKEN_RE.match(formula, pos)
        if not m:
            raise ValueError(f"Cannot tokenize formula at position {pos}: {formula!r}")
        pos = m.end()
        kind = m.lastgroup
        if kind == "WS":
            continue
        tokens.append((kind, m.group()))
    return tokens


class Parser:
    """Recursive-descent parser. Precedence (low to high): OR, AND, NOT, atom.
    This matches Python's own `not` > `and` > `or` binding, i.e. the natural
    reading of expressions like 'A & B | C & !D' as '(A&B) | (C&!D)'.
    """

    def __init__(self, tokens):
        self.tokens = tokens
        self.i = 0

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else (None, None)

    def advance(self):
        tok = self.tokens[self.i]
        self.i += 1
        return tok

    def parse(self):
        node = self.parse_or()
        if self.i != len(self.tokens):
            raise ValueError(f"Unexpected trailing tokens: {self.tokens[self.i:]}")
        return node

    def parse_or(self):
        node = self.parse_and()
        while self.peek()[0] == "OR":
            self.advance()
            rhs = self.parse_and()
            node = ("OR", node, rhs)
        return node

    def parse_and(self):
        node = self.parse_not()
        while self.peek()[0] == "AND":
            self.advance()
            rhs = self.parse_not()
            node = ("AND", node, rhs)
        return node

    def parse_not(self):
        if self.peek()[0] == "NOT":
            self.advance()
            return ("NOT", self.parse_not())
        return self.parse_atom()

    def parse_atom(self):
        kind, text = self.peek()
        if kind == "LPAREN":
            self.advance()
            node = self.parse_or()
            kind2, _ = self.advance()
            if kind2 != "RPAREN":
                raise ValueError("Expected closing parenthesis")
            return node
        if kind == "IDENT":
            self.advance()
            if ":" in text:
                name, thr = text.split(":")
                return ("REF", name, int(thr))
            return ("REF", text, None)
        raise ValueError(f"Unexpected token: {kind} {text!r}")


def parse_formula(formula):
    return Parser(tokenize(formula)).parse()


def collect_refs(ast, out):
    tag = ast[0]
    if tag == "REF":
        out.add(ast[1])
    elif tag == "NOT":
        collect_refs(ast[1], out)
    else:  # AND / OR
        collect_refs(ast[1], out)
        collect_refs(ast[2], out)


def eval_ast(ast, state, aliases):
    tag = ast[0]
    if tag == "REF":
        name = aliases.get(ast[1], ast[1])
        threshold = ast[2] if ast[2] is not None else 1
        return state.get(name, 0) >= threshold
    if tag == "NOT":
        return not eval_ast(ast[1], state, aliases)
    if tag == "AND":
        return eval_ast(ast[1], state, aliases) and eval_ast(ast[2], state, aliases)
    if tag == "OR":
        return eval_ast(ast[1], state, aliases) or eval_ast(ast[2], state, aliases)
    raise ValueError(f"Unknown AST node: {ast}")


class LogicalModel:
    """A multi-valued logical network built from (node, value, formula) rules."""

    def __init__(self, rules, aliases=None):
        self.aliases = aliases or {}
        self.compiled = defaultdict(list)  # node -> [(value, ast), ...] desc by value
        all_targets = set()
        all_refs = set()
        for node, value, formula, _desc in rules:
            ast = parse_formula(formula)
            self.compiled[node].append((value, ast))
            all_targets.add(node)
            collect_refs(ast, all_refs)
        for node in self.compiled:
            self.compiled[node].sort(key=lambda pair: -pair[0])

        resolved_refs = {self.aliases.get(r, r) for r in all_refs}
        self.rule_nodes = sorted(all_targets)
        self.input_nodes = sorted(resolved_refs - all_targets)
        self.all_nodes = sorted(all_targets | resolved_refs)

    def evaluate_node(self, node, state):
        for value, ast in self.compiled.get(node, []):
            if eval_ast(ast, state, self.aliases):
                return value
        return 0

    def async_run(self, initial_state, n_steps, rng):
        """Run one asynchronous trajectory. Returns (final_state, history),
        history is a list of state-dict snapshots (one per sweep)."""
        state = dict(initial_state)
        for n in self.all_nodes:
            state.setdefault(n, 0)
        history = [dict(state)]
        order_pool = list(self.rule_nodes)
        steps_per_sweep = len(order_pool)
        for step in range(n_steps):
            node = order_pool[step % steps_per_sweep]
            if step % steps_per_sweep == 0:
                rng.shuffle(order_pool)
            state[node] = self.evaluate_node(node, state)
            if step % steps_per_sweep == steps_per_sweep - 1:
                history.append(dict(state))
        return state, history

    def ensemble_run(self, initial_state, n_sweeps=8, n_runs=200, seed=0):
        """Run many independent asynchronous trajectories from the same
        initial condition and return:
          - final_states: list of dicts (one per run)
          - traj: dict node -> list of average-active-fraction per sweep
        """
        rng = random.Random(seed)
        n_steps = n_sweeps * len(self.rule_nodes)
        finals = []
        sweep_sums = None
        n_sweep_points = n_sweeps + 1
        for run in range(n_runs):
            final_state, history = self.async_run(initial_state, n_steps, rng)
            finals.append(final_state)
            if sweep_sums is None:
                sweep_sums = {n: [0.0] * n_sweep_points for n in self.all_nodes}
            for t, snap in enumerate(history):
                for n in self.all_nodes:
                    if snap.get(n, 0) >= 1:
                        sweep_sums[n][t] += 1.0
        traj = {n: [v / n_runs for v in vals] for n, vals in sweep_sums.items()}
        return finals, traj
