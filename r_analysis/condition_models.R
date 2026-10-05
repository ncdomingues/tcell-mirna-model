# Shared model-building helpers, mirroring src/run_analysis.py
# (build_condition_models/base_initial_state) and src/validate_biology.py
# (build_polarisation_model/LIGANDS/BASE_INPUTS) exactly, so every R result
# is directly comparable to the already-validated Python output.

BASE_INPUTS <- c("APC-Antigen", "CD3", "CD45RA", "CGC", "CREBBP")
LIGANDS <- c("IL12", "IL4", "TGFB", "IL6")
MASTER_TFS <- c(TBX21 = "Th1", GATA3 = "Th2", RORC = "Th17", FOXP3 = "iTreg",
                 BCL6 = "Tfh", AHR = "Th22", `SPI1-PU1` = "Th9")

# Build one of the 4 standard miRNA conditions from the full rules_df,
# excluding any `extra_excluded` node rules (e.g. LIGANDS, to promote them
# to fixed exogenous inputs) and any miR-34c-5p rule swap requested.
build_condition <- function(rules_df, aliases, mirna_condition,
                             extra_excluded = character(0),
                             mir34c_override = NULL) {
  excluded <- switch(mirna_condition,
    "both" = character(0),
    "miR-34c-5p only" = "miR-155-5p",
    "miR-155-5p only" = "miR-34c-5p",
    "no miRNA" = c("miR-34c-5p", "miR-155-5p"),
    stop("unknown mirna_condition: ", mirna_condition)
  )
  excluded <- union(excluded, extra_excluded)
  df <- rules_df[!(rules_df$node %in% excluded), , drop = FALSE]

  if (!is.null(mir34c_override)) {
    # replace ALL existing miR-34c-5p rows with a single new Boolean rule
    df <- df[df$node != "miR-34c-5p", , drop = FALSE]
    df <- rbind(df, data.frame(node = "miR-34c-5p", value = 1L,
                                formula = mir34c_override,
                                description = "TF-combination hypothesis (2026 re-test)",
                                stringsAsFactors = FALSE))
  }
  build_model(df, aliases)
}

base_initial_state <- function(model, mirna_condition, forced_ligands = character(0)) {
  state <- as.list(setNames(rep(0L, length(model$all_nodes)), model$all_nodes))
  for (inp in BASE_INPUTS) if (!is.null(state[[inp]])) state[[inp]] <- 1L
  if (mirna_condition %in% c("both", "miR-34c-5p only") && !is.null(state[["miR-34c-5p"]])) {
    # miR-34c-5p/miR-155-5p are rule-governed unless swapped to override;
    # leave at 0 initial and let the async run compute it from rules.
  }
  for (lig in LIGANDS) {
    if (!is.null(state[[lig]])) state[[lig]] <- if (lig %in% forced_ligands) 1L else 0L
  }
  state
}

# Summarise an ensemble's finals into % active + Wilson CI for a set of nodes.
summarise_finals <- function(finals, nodes) {
  n <- length(finals)
  rows <- lapply(nodes, function(nd) {
    k <- sum(vapply(finals, function(f) { v <- f[[nd]]; if (is.null(v)) 0 else (v >= 1) }, logical(1)))
    ci <- wilson_ci(k, n)
    data.frame(node = nd, n = n, k = k, pct = round(100 * k / n, 1),
               ci_lo = round(ci[1], 1), ci_hi = round(ci[2], 1))
  })
  do.call(rbind, rows)
}
