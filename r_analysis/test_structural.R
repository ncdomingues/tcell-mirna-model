suppressMessages({
  library(readr)
})
source("logic_engine.R")
source("load_rules.R")

d <- load_rules_data("data")
stopifnot(nrow(d$rules_df) == 89)
model <- build_model(d$rules_df, d$aliases)

cat("n rules (rows):", nrow(d$rules_df), "\n")
cat("n rule_nodes:", length(model$rule_nodes), "\n")
cat("n input_nodes:", length(model$input_nodes), "\n")
cat("n all_nodes:", length(model$all_nodes), "\n")
cat("input_nodes:", paste(model$input_nodes, collapse = ", "), "\n")

# sanity: evaluate a known rule by hand -- TCR requires APC-Antigen & CD3 & CD28 & !IL2:2
state <- as.list(setNames(rep(0L, length(model$all_nodes)), model$all_nodes))
state[["APC-Antigen"]] <- 1L; state[["CD3"]] <- 1L; state[["CD28"]] <- 1L; state[["IL2"]] <- 0L
cat("TCR should activate (1):", evaluate_node(model, "TCR", state), "\n")
state[["IL2"]] <- 2L
cat("TCR should NOT activate with IL2>=2 (0):", evaluate_node(model, "TCR", state), "\n")
