# Re-tests the miR-34c-5p transcriptional-regulator hypotheses actually
# explored by hand in 2016-2017 (Modelling/2017-mymodel/mir-34c tf
# variations/formulae_*.txt), against the FINAL published rule (Table 5,
# 2019 thesis) as the reference/baseline. The 2017 exploration varied
# several rules simultaneously between snapshots (its own miR-155-5p rule,
# an IL4 rule, a STAT3 rule and a GATA3 rule all drifted alongside the
# intended miR-34c-5p change -- confirmed by diffing the formulae_*.txt
# files directly), which confounds any single-factor conclusion. This
# re-test changes ONE thing at a time: only the miR-34c-5p rule is swapped
# in, against the otherwise-unmodified, already-validated 89-rule network,
# and run as a proper 300-run asynchronous ensemble with Wilson 95% CIs --
# neither of which the original ad hoc base-R plots(.) calls provided.
suppressMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr); library(jsonlite)
})
source("logic_engine.R")
source("load_rules.R")
source("condition_models.R")

set.seed(2)
N_SWEEPS <- 10
N_RUNS <- 300

d <- load_rules_data("data")
master_tfs <- fromJSON("data/master_tfs.json")

# --- the 7 historical hypotheses, transcribed verbatim from
# Modelling/2017-mymodel/mir-34c tf variations/formulae_<COMBO>.txt
hypotheses <- list(
  "Published (Table 5, 2019)"        = NULL,  # no override -- the real 2-value AND-rule
  "2017: GATA3 only"                 = "GATA3 & !STAT3",
  "2017: MYC only"                   = "MYC & !STAT3",
  "2017: TP53 or FOXO3"              = "(TP53 | FOXO3) & !STAT3",
  "2017: MYC and (TP53 or FOXO3)"    = "MYC & (TP53 | FOXO3) & !STAT3",
  "2017: GATA3 and (TP53 or FOXO3)"  = "GATA3 & (TP53 | FOXO3) & !STAT3",
  "2017: (GATA3 or MYC) and (TP53 or FOXO3)" = "(GATA3 | MYC) & (TP53 | FOXO3) & !STAT3",
  "2017: GATA3, MYC or FOS"          = "(GATA3 | MYC | FOS) & !STAT3"
)

watch_nodes <- unique(c(names(master_tfs), "miR-34c-5p", "STAT3", "GATA3", "MYC", "TP53", "FOXO3", "FOS"))

cat("Running", length(hypotheses), "miR-34c-5p regulatory hypotheses (both-miRNA context)...\n")
all_rows <- list()
for (hyp_name in names(hypotheses)) {
  override <- hypotheses[[hyp_name]]
  model <- build_condition(d$rules_df, d$aliases, "both", mir34c_override = override)
  state0 <- base_initial_state(model, "both")
  ens <- ensemble_run(model, state0, n_sweeps = N_SWEEPS, n_runs = N_RUNS, seed = 42)
  summ <- summarise_finals(ens$finals, watch_nodes)
  summ$hypothesis <- hyp_name
  all_rows[[hyp_name]] <- summ
  cat("  done:", hyp_name, "\n")
}
res <- bind_rows(all_rows)
write_tsv(res, "data/master_regulator_variations.tsv")

res$hypothesis <- factor(res$hypothesis, levels = rev(names(hypotheses)))

# ---- Figure A: master-TF activation heatmap across hypotheses
tf_res <- res %>% filter(node %in% names(master_tfs)) %>%
  mutate(th_label = paste0(node, " (", master_tfs[node], ")"))

p1 <- ggplot(tf_res, aes(x = th_label, y = hypothesis, fill = pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", pct)), size = 3, color = "black") +
  scale_fill_gradientn(colors = c("#F3F4EF", "#8AA6C7", "#1F6F63"), limits = c(0, 100), name = "% runs\nactive") +
  labs(
    title = "miR-34c-5p regulatory hypothesis vs. Th master-regulator outcome",
    subtitle = "300-run asynchronous ensembles; only the miR-34c-5p rule differs between rows\n(everything else fixed at the published, validated Table 5 network)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold", size = 11))

ggsave("figures/master_regulator_heatmap.png", p1, width = 10, height = 5.5, dpi = 170, bg = "white")

# ---- Figure B: does the hypothesis even make miR-34c-5p active? (necessary precondition)
mir_res <- res %>% filter(node == "miR-34c-5p")
p2 <- ggplot(mir_res, aes(x = hypothesis, y = pct)) +
  geom_col(fill = "#B23B5E") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, color = "grey30") +
  coord_flip(ylim = c(0, 100)) +
  labs(
    title = "How often is miR-34c-5p itself active under each hypothesis?",
    subtitle = "95% Wilson CI, 300 runs each",
    x = NULL, y = "% of runs with miR-34c-5p active"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 11))

ggsave("figures/master_regulator_mir34c_activity.png", p2, width = 8, height = 5, dpi = 170, bg = "white")

cat("\nSaved figures/master_regulator_heatmap.png and figures/master_regulator_mir34c_activity.png\n")
cat("\n=== miR-34c-5p activation by hypothesis ===\n")
print(mir_res %>% select(hypothesis, pct, ci_lo, ci_hi) %>% arrange(desc(pct)))
