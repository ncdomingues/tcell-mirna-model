# -*- coding: utf-8 -*-
"""
Unit tests for the logic engine (tokenizer, parser, evaluator, LogicalModel)
-- independent of whether the digitized biology is "correct", these confirm
the code implements the notation it claims to implement.

Run with:  python -m unittest tests.test_logic_engine -v
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from logic_engine import tokenize, parse_formula, eval_ast, collect_refs, LogicalModel
from rules_data import RULES, ALIASES


class TestTokenizer(unittest.TestCase):
    def test_simple_tokens(self):
        toks = tokenize("A & !B | (C:2)")
        kinds = [k for k, _ in toks]
        self.assertEqual(kinds, ["IDENT", "AND", "NOT", "IDENT", "OR", "LPAREN", "IDENT", "RPAREN"])

    def test_hyphenated_identifiers(self):
        toks = tokenize("AKT1-PKB & miR-34c-5p & !NFKBIE-IkB")
        idents = [t for k, t in toks if k == "IDENT"]
        self.assertEqual(idents, ["AKT1-PKB", "miR-34c-5p", "NFKBIE-IkB"])

    def test_threshold_suffix_stays_attached(self):
        toks = tokenize("IL2R:2 & STAT5:1")
        idents = [t for k, t in toks if k == "IDENT"]
        self.assertEqual(idents, ["IL2R:2", "STAT5:1"])

    def test_whitespace_ignored(self):
        a = tokenize("A&B|C")
        b = tokenize("A  &  B  |  C")
        self.assertEqual([t for _, t in a], [t for _, t in b])


class TestParserPrecedence(unittest.TestCase):
    """NOT > AND > OR, matching Python's own not/and/or binding."""

    def test_and_binds_tighter_than_or(self):
        # A & B | C  ==  (A&B) | C, not A & (B|C)
        ast = parse_formula("A & B | C")
        self.assertEqual(ast[0], "OR")
        self.assertEqual(ast[1], ("AND", ("REF", "A", None), ("REF", "B", None)))
        self.assertEqual(ast[2], ("REF", "C", None))

    def test_not_binds_tighter_than_and(self):
        ast = parse_formula("!A & B")
        self.assertEqual(ast[0], "AND")
        self.assertEqual(ast[1], ("NOT", ("REF", "A", None)))

    def test_parentheses_override_precedence(self):
        ast = parse_formula("A & (B | C)")
        self.assertEqual(ast[0], "AND")
        self.assertEqual(ast[2][0], "OR")

    def test_real_formula_matches_documented_reading(self):
        # IL2 value-2 rule; MODEL_NOTES.md documents this exact grouping
        formula = "(NFAT & NFKB) | (JUN & FOS) & STAT5 & !(STAT5 & STAT6) & !(NFKB & TBX21)"
        ast = parse_formula(formula)
        self.assertEqual(ast[0], "OR")
        self.assertEqual(ast[1], ("AND", ("REF", "NFAT", None), ("REF", "NFKB", None)))

    def test_documented_ambiguous_sirt1_reading(self):
        # MODEL_NOTES.md: SIRT1 parses as FOXO3 | (rest), not (FOXO3 & rest) --
        # this test locks in that documented (if surprising) reading so a
        # future refactor can't silently change it.
        ast = parse_formula("FOXO3 | (FOXO3 & CREB) & !TP53 & !miR-155-5p")
        self.assertEqual(ast[0], "OR")
        self.assertEqual(ast[1], ("REF", "FOXO3", None))


class TestEvaluator(unittest.TestCase):
    def test_and_or_not(self):
        ast = parse_formula("A & !B | C")
        self.assertTrue(eval_ast(ast, {"A": 1, "B": 0, "C": 0}, {}))
        self.assertFalse(eval_ast(ast, {"A": 1, "B": 1, "C": 0}, {}))
        self.assertTrue(eval_ast(ast, {"A": 0, "B": 1, "C": 1}, {}))

    def test_bare_reference_is_truthy_at_1(self):
        ast = parse_formula("A")
        self.assertFalse(eval_ast(ast, {"A": 0}, {}))
        self.assertTrue(eval_ast(ast, {"A": 1}, {}))
        self.assertTrue(eval_ast(ast, {"A": 2}, {}))

    def test_threshold_reference_requires_level(self):
        ast = parse_formula("A:2")
        self.assertFalse(eval_ast(ast, {"A": 1}, {}))
        self.assertTrue(eval_ast(ast, {"A": 2}, {}))

    def test_missing_node_defaults_to_zero(self):
        # B is absent from the state entirely -> treated as 0 -> !B is True
        ast = parse_formula("A & !B")
        self.assertTrue(eval_ast(ast, {"A": 1}, {}))
        ast2 = parse_formula("B")
        self.assertFalse(eval_ast(ast2, {"A": 1}, {}))

    def test_alias_resolution(self):
        ast = parse_formula("AKT")
        self.assertTrue(eval_ast(ast, {"AKT1-PKB": 1}, {"AKT": "AKT1-PKB"}))
        self.assertFalse(eval_ast(ast, {"AKT": 1}, {"AKT": "AKT1-PKB"}))  # must resolve, not match literal


class TestLogicalModelThresholdSemantics(unittest.TestCase):
    def setUp(self):
        rules = [
            ("X", 2, "A & B", "high"),
            ("X", 1, "A", "low"),
        ]
        self.model = LogicalModel(rules, {})

    def test_highest_satisfied_value_wins(self):
        self.assertEqual(self.model.evaluate_node("X", {"A": 1, "B": 1}), 2)
        self.assertEqual(self.model.evaluate_node("X", {"A": 1, "B": 0}), 1)
        self.assertEqual(self.model.evaluate_node("X", {"A": 0, "B": 0}), 0)

    def test_input_node_detection(self):
        self.assertEqual(self.model.rule_nodes, ["X"])
        self.assertEqual(self.model.input_nodes, ["A", "B"])


class TestFullRuleTableParses(unittest.TestCase):
    """Every formula actually transcribed from the thesis must parse without error."""

    def test_all_rules_parse(self):
        for node, value, formula, _desc in RULES:
            try:
                parse_formula(formula)
            except Exception as e:  # pragma: no cover
                self.fail(f"Formula for {node}={value} failed to parse: {formula!r} ({e})")

    def test_model_builds_without_error(self):
        model = LogicalModel(RULES, ALIASES)
        self.assertGreater(len(model.rule_nodes), 0)


class TestStructuralIntegrity(unittest.TestCase):
    """Sanity-checks the digitization itself against the thesis's Table 5,
    independent of whether the biology it encodes is correct."""

    def setUp(self):
        self.model = LogicalModel(RULES, ALIASES)

    def test_rule_row_count_matches_thesis_table(self):
        # Table 5 has 89 (node, value, formula) rows (legend row excluded)
        self.assertEqual(len(RULES), 89)

    def test_rule_governed_node_count(self):
        self.assertEqual(len(self.model.rule_nodes), 76)

    def test_expected_input_nodes(self):
        # Every node referenced in a formula but never assigned a rule --
        # these are the only nodes with no upstream regulation in Table 5.
        expected = {"APC-Antigen", "CD3", "CD45RA", "CGC", "CREBBP"}
        self.assertEqual(set(self.model.input_nodes), expected)

    def test_no_undefined_alias_targets(self):
        # every ALIASES value must itself resolve to a real rule-governed node
        for raw, resolved in ALIASES.items():
            self.assertIn(resolved, self.model.rule_nodes,
                           f"alias {raw!r} -> {resolved!r} does not point at a rule-governed node")

    def test_both_mirnas_are_rule_governed_not_inputs(self):
        self.assertIn("miR-34c-5p", self.model.rule_nodes)
        self.assertIn("miR-155-5p", self.model.rule_nodes)


if __name__ == "__main__":
    unittest.main(verbosity=2)
