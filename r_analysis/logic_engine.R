# -*- coding: utf-8 -*-
# Multi-valued logical (Boolean-network-style) simulation engine in R,
# independently reimplementing the exact semantics of src/logic_engine.py
# (tokenize -> recursive-descent parse -> AST -> evaluate), for the same
# rule notation: & = AND, | = OR, ! = NOT, NODE:k = "NODE has reached
# value >= k".
#
# This deliberately does NOT use the original 2017 approach found in
# Modelling/Base files/GinSimR_ND.R and its siblings, which built each
# rule's truth value via gsub() substring substitution directly on the
# formula TEXT (replacing e.g. "TP53" with "1" or "0" and then
# eval(parse(text=...))). That approach has a real, demonstrated failure
# mode: gsub() has no word-boundary awareness, so adjacent node names can
# concatenate into unintended numeric literals -- e.g. a rule containing
# "TGFBRIL6R" (a missing "|" between TGFBR and IL6R) silently evaluates as
# a two-digit number the interpreter coerces to TRUE, instead of raising an
# error. A proper tokenizer treats each identifier as an atomic unit and
# cannot make that mistake.

# ---------------------------------------------------------------- tokenizer
tokenize <- function(formula) {
  pattern <- "\\s+|\\(|\\)|&|\\||!|[A-Za-z][A-Za-z0-9\\-]*(:\\d+)?"
  matches <- regmatches(formula, gregexpr(pattern, formula, perl = TRUE))[[1]]
  tokens <- list()
  for (m in matches) {
    if (grepl("^\\s+$", m)) next
    kind <- if (m == "(") "LPAREN"
      else if (m == ")") "RPAREN"
      else if (m == "&") "AND"
      else if (m == "|") "OR"
      else if (m == "!") "NOT"
      else "IDENT"
    tokens[[length(tokens) + 1]] <- list(kind = kind, text = m)
  }
  tokens
}

# ---------------------------------------------------------------- parser
# Precedence (low to high): OR, AND, NOT, atom -- matches Python's own
# `not` > `and` > `or` binding, i.e. 'A & B | C & !D' reads as
# '(A&B) | (C&!D)'.
new_parser <- function(tokens) {
  list(tokens = tokens, i = 1)
}

p_peek <- function(p) if (p$i <= length(p$tokens)) p$tokens[[p$i]] else list(kind = NA, text = NA)

parse_formula <- function(formula) {
  tokens <- tokenize(formula)
  state <- new.env()
  state$tokens <- tokens
  state$i <- 1

  peek <- function() if (state$i <= length(state$tokens)) state$tokens[[state$i]] else list(kind = NA, text = NA)
  advance <- function() { tok <- state$tokens[[state$i]]; state$i <- state$i + 1; tok }

  parse_or <- function() {
    node <- parse_and()
    while (identical(peek()$kind, "OR")) {
      advance()
      rhs <- parse_and()
      node <- list(tag = "OR", left = node, right = rhs)
    }
    node
  }
  parse_and <- function() {
    node <- parse_not()
    while (identical(peek()$kind, "AND")) {
      advance()
      rhs <- parse_not()
      node <- list(tag = "AND", left = node, right = rhs)
    }
    node
  }
  parse_not <- function() {
    if (identical(peek()$kind, "NOT")) {
      advance()
      return(list(tag = "NOT", inner = parse_not()))
    }
    parse_atom()
  }
  parse_atom <- function() {
    tok <- peek()
    if (identical(tok$kind, "LPAREN")) {
      advance()
      node <- parse_or()
      close <- advance()
      if (!identical(close$kind, "RPAREN")) stop("Expected closing parenthesis in: ", formula)
      return(node)
    }
    if (identical(tok$kind, "IDENT")) {
      advance()
      if (grepl(":", tok$text, fixed = TRUE)) {
        parts <- strsplit(tok$text, ":", fixed = TRUE)[[1]]
        return(list(tag = "REF", name = parts[1], thr = as.integer(parts[2])))
      }
      return(list(tag = "REF", name = tok$text, thr = NA_integer_))
    }
    stop("Unexpected token parsing formula: ", formula)
  }

  ast <- parse_or()
  if (state$i <= length(state$tokens)) stop("Unexpected trailing tokens in: ", formula)
  ast
}

collect_refs <- function(ast, acc) {
  if (ast$tag == "REF") {
    acc[[ast$name]] <- TRUE
  } else if (ast$tag == "NOT") {
    acc <- collect_refs(ast$inner, acc)
  } else {
    acc <- collect_refs(ast$left, acc)
    acc <- collect_refs(ast$right, acc)
  }
  acc
}

eval_ast <- function(ast, state, aliases) {
  tag <- ast$tag
  if (tag == "REF") {
    name <- aliases[[ast$name]]
    if (is.null(name)) name <- ast$name
    thr <- if (is.na(ast$thr)) 1L else ast$thr
    val <- state[[name]]
    if (is.null(val)) val <- 0L
    return(val >= thr)
  }
  if (tag == "NOT") return(!eval_ast(ast$inner, state, aliases))
  if (tag == "AND") return(eval_ast(ast$left, state, aliases) && eval_ast(ast$right, state, aliases))
  if (tag == "OR") return(eval_ast(ast$left, state, aliases) || eval_ast(ast$right, state, aliases))
  stop("Unknown AST tag: ", tag)
}

# ---------------------------------------------------------------- model
# rules_df: data.frame(node, value, formula, description)
# aliases: named list, alias -> canonical
build_model <- function(rules_df, aliases = list()) {
  compiled <- list()
  all_targets <- character(0)
  all_refs <- character(0)

  for (i in seq_len(nrow(rules_df))) {
    node <- rules_df$node[i]
    value <- as.integer(rules_df$value[i])
    formula <- rules_df$formula[i]
    ast <- parse_formula(formula)
    if (is.null(compiled[[node]])) compiled[[node]] <- list()
    compiled[[node]][[length(compiled[[node]]) + 1]] <- list(value = value, ast = ast)
    all_targets <- union(all_targets, node)
    refs <- collect_refs(ast, list())
    all_refs <- union(all_refs, names(refs))
  }

  # sort each node's rule alternatives by decreasing value (multi-valued
  # threshold semantics: highest satisfied value wins)
  for (node in names(compiled)) {
    vals <- vapply(compiled[[node]], function(r) r$value, integer(1))
    compiled[[node]] <- compiled[[node]][order(-vals)]
  }

  resolved_refs <- unique(vapply(all_refs, function(r) {
    a <- aliases[[r]]
    if (is.null(a)) r else a
  }, character(1)))

  rule_nodes <- sort(all_targets)
  input_nodes <- sort(setdiff(resolved_refs, all_targets))
  all_nodes <- sort(union(all_targets, resolved_refs))

  list(compiled = compiled, aliases = aliases, rule_nodes = rule_nodes,
       input_nodes = input_nodes, all_nodes = all_nodes)
}

evaluate_node <- function(model, node, state) {
  rules <- model$compiled[[node]]
  if (is.null(rules)) return(0L)
  for (r in rules) {
    if (eval_ast(r$ast, state, model$aliases)) return(r$value)
  }
  0L
}

# One asynchronous trajectory. state0: named list/vector of node -> value.
# Returns list(final = named list, history = list of named-list snapshots,
# one per sweep including the initial state as sweep 0).
async_run <- function(model, state0, n_steps, rng_state = NULL) {
  state <- as.list(state0)
  for (n in model$all_nodes) if (is.null(state[[n]])) state[[n]] <- 0L

  history <- list(state)
  order_pool <- model$rule_nodes
  steps_per_sweep <- length(order_pool)

  for (step in seq_len(n_steps)) {
    idx0 <- (step - 1) %% steps_per_sweep
    if (idx0 == 0) order_pool <- sample(order_pool)  # reshuffle each sweep
    node <- order_pool[idx0 + 1]
    state[[node]] <- evaluate_node(model, node, state)
    if (idx0 == steps_per_sweep - 1) history[[length(history) + 1]] <- state
  }
  list(final = state, history = history)
}

# Monte Carlo ensemble. Returns list(finals = list of state lists,
# traj = data.frame(node, sweep, avg_active_fraction)).
ensemble_run <- function(model, state0, n_sweeps = 10, n_runs = 200, seed = 0) {
  set.seed(seed)
  n_steps <- n_sweeps * length(model$rule_nodes)
  n_sweep_points <- n_sweeps + 1
  finals <- vector("list", n_runs)
  sums <- matrix(0, nrow = length(model$all_nodes), ncol = n_sweep_points,
                  dimnames = list(model$all_nodes, NULL))

  for (run in seq_len(n_runs)) {
    res <- async_run(model, state0, n_steps)
    finals[[run]] <- res$final
    for (t in seq_along(res$history)) {
      snap <- res$history[[t]]
      for (nd in model$all_nodes) {
        v <- snap[[nd]]
        if (!is.null(v) && v >= 1) sums[nd, t] <- sums[nd, t] + 1
      }
    }
  }
  traj <- sums / n_runs
  list(finals = finals, traj = traj)
}

# Convenience: % of runs where `node` is active (>=1) at the final sweep.
activation_pct <- function(finals, node) {
  vals <- vapply(finals, function(f) { v <- f[[node]]; if (is.null(v)) 0 else v }, numeric(1))
  100 * mean(vals >= 1)
}

# Wilson 95% CI for a binomial proportion, k successes out of n, returns c(lo, hi) as percentages.
wilson_ci <- function(k, n, z = 1.96) {
  if (n == 0) return(c(0, 0))
  phat <- k / n
  denom <- 1 + z^2 / n
  centre <- phat + z^2 / (2 * n)
  adj <- z * sqrt((phat * (1 - phat) + z^2 / (4 * n)) / n)
  lo <- (centre - adj) / denom
  hi <- (centre + adj) / denom
  100 * c(max(0, lo), min(1, hi))
}
