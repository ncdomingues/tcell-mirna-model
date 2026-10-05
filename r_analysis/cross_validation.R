# Cross-implementation validation: run the independent R engine on the SAME
# 89-rule table and the SAME 4 miRNA conditions already reported in the
# thesis (Results 3.3.4 / Table 6), and check whether two independently
# written simulators (Python, this R port) agree on the resulting steady-
# state distributions. Agreement is the right bar here, not literal
# RNG-for-RNG reproduction -- each engine uses its own random number stream.
suppressMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr); library(jsonlite)
})
source("logic_engine.R")
source("load_rules.R")
source("condition_models.R")

set.seed(1)
N_SWEEPS <- 10
N_RUNS <- 300  # R is slower per-call than the Python engine; still a 1.7-4.9pp Wilson half-width at this n

d <- load_rules_data("data")
watchlist <- fromJSON("data/watchlist.json")
conditions <- c("no miRNA", "miR-155-5p only", "miR-34c-5p only", "both")

cat("Running R-engine ensembles (", N_RUNS, "runs x", N_SWEEPS, "sweeps x", length(conditions), "conditions)...\n")
r_results <- list()
for (cond in conditions) {
  model <- build_condition(d$rules_df, d$aliases, cond)
  state0 <- base_initial_state(model, cond)
  ens <- ensemble_run(model, state0, n_sweeps = N_SWEEPS, n_runs = N_RUNS, seed = 42)
  summ <- summarise_finals(ens$finals, watchlist)
  summ$condition <- cond
  r_results[[cond]] <- summ
  cat("  done:", cond, "\n")
}
r_df <- bind_rows(r_results)
write_tsv(r_df, "data/r_engine_results.tsv")

# ---- load Python's already-published results (output/full_results_with_ci.csv)
py_raw <- read_csv("../output/full_results_with_ci.csv", show_col_types = FALSE)
py_long <- py_raw %>%
  pivot_longer(-node, names_to = c("condition", ".value"), names_pattern = "(.*)_(pct|ci)") %>%
  mutate(
    condition = recode(condition,
      `no miRNA` = "no miRNA", `miR-155-5p only` = "miR-155-5p only",
      `miR-34c-5p only` = "miR-34c-5p only", `both` = "both"),
    py_pct = pct
  ) %>%
  separate(ci, into = c("py_ci_lo", "py_ci_hi"), sep = "-", convert = TRUE) %>%
  select(node, condition, py_pct, py_ci_lo, py_ci_hi) %>%
  filter(node %in% watchlist)

merged <- r_df %>%
  rename(r_pct = pct, r_ci_lo = ci_lo, r_ci_hi = ci_hi) %>%
  inner_join(py_long, by = c("node", "condition")) %>%
  mutate(
    ci_overlap = !(r_ci_hi < py_ci_lo | py_ci_hi < r_ci_lo),
    abs_diff = abs(r_pct - py_pct)
  )

write_tsv(merged, "data/cross_validation_merged.tsv")

n_total <- nrow(merged)
n_overlap <- sum(merged$ci_overlap)
cat(sprintf("\nCI overlap: %d / %d node-condition pairs (%.1f%%)\n",
            n_overlap, n_total, 100 * n_overlap / n_total))
cat(sprintf("Mean abs difference: %.2f pp, max: %.2f pp\n",
            mean(merged$abs_diff), max(merged$abs_diff)))
print(merged %>% filter(!ci_overlap) %>% arrange(desc(abs_diff)))

# ---- ggplot: R% vs Python% agreement scatter, faceted by condition
p <- ggplot(merged, aes(x = py_pct, y = r_pct, color = ci_overlap)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = r_ci_lo, ymax = r_ci_hi), width = 0, alpha = 0.4) +
  geom_errorbarh(aes(xmin = py_ci_lo, xmax = py_ci_hi), height = 0, alpha = 0.4) +
  geom_point(size = 1.8) +
  facet_wrap(~condition, ncol = 2) +
  scale_color_manual(values = c(`TRUE` = "#1F6F63", `FALSE` = "#B23B5E"),
                      labels = c(`TRUE` = "95% CIs overlap", `FALSE` = "CIs do not overlap"),
                      name = NULL) +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(
    title = "Independent R reimplementation vs. the Python ensemble (20-node watchlist x 4 conditions)",
    subtitle = sprintf("n=%d runs per condition (R) vs n=400 (Python); dashed line = perfect agreement", N_RUNS),
    x = "Python engine: % of runs active (95% Wilson CI)",
    y = "R engine: % of runs active (95% Wilson CI)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11))

ggsave("figures/cross_validation_scatter.png", p, width = 9, height = 8, dpi = 170, bg = "white")
cat("\nSaved figures/cross_validation_scatter.png\n")
