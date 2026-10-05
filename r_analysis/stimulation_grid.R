# Different-input-settings analysis, mirroring the stimulation-condition
# folders actually explored in Modelling/Naldi 2010/{Th1,Th2,Th17,TReg}
# stimulation and Modelling/AJ 2015/{Th1,Th2,Th17,Treg} stimulation
# (filename pattern: which miRNA(s) x which cytokine cocktail on top of
# baseline IL2), but run here as a full 4x4=16-condition grid on the final
# 89-rule integrative model, with proper ensembles + Wilson CIs -- the
# original .zginml files in those folders are individual GINsim steady-
# state exports, one file per condition, with no aggregate comparison
# across the grid ever assembled.
#
# Ligand/master-TF convention matches src/validate_biology.py exactly
# (same LIGANDS list, same forced-exogenous-input mechanism), so this
# extends that already-published validation (Results Table 9, no-miRNA
# only) across all 4 miRNA conditions instead of just one.
suppressMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr); library(jsonlite)
})
source("logic_engine.R")
source("load_rules.R")
source("condition_models.R")

set.seed(3)
N_SWEEPS <- 10
N_RUNS <- 300

d <- load_rules_data("data")
master_tfs <- fromJSON("data/master_tfs.json")

stim_conditions <- list(
  "IL2 only"        = character(0),
  "Th1 (+IL12)"      = "IL12",
  "Th2 (+IL4)"       = "IL4",
  "Th17 (+IL6+TGFB)" = c("IL6", "TGFB"),
  "Treg (+TGFB)"     = "TGFB"
)
mirna_conditions <- c("no miRNA", "miR-155-5p only", "miR-34c-5p only", "both")

cat("Running", length(stim_conditions), "x", length(mirna_conditions),
    "=", length(stim_conditions) * length(mirna_conditions), "conditions...\n")

all_rows <- list()
for (stim_name in names(stim_conditions)) {
  forced <- stim_conditions[[stim_name]]
  for (mirna_cond in mirna_conditions) {
    model <- build_condition(d$rules_df, d$aliases, mirna_cond, extra_excluded = LIGANDS)
    state0 <- base_initial_state(model, mirna_cond, forced_ligands = forced)
    ens <- ensemble_run(model, state0, n_sweeps = N_SWEEPS, n_runs = N_RUNS, seed = 42)
    summ <- summarise_finals(ens$finals, names(master_tfs))
    summ$stimulation <- stim_name
    summ$mirna <- mirna_cond
    key <- paste(stim_name, mirna_cond, sep = " | ")
    all_rows[[key]] <- summ
    cat("  done:", key, "\n")
  }
}
res <- bind_rows(all_rows) %>%
  mutate(th_label = paste0(node, " (", master_tfs[node], ")"))
write_tsv(res, "data/stimulation_grid.tsv")

res$stimulation <- factor(res$stimulation, levels = names(stim_conditions))
res$mirna <- factor(res$mirna, levels = mirna_conditions)

# ---- Figure: faceted heatmap, one panel per stimulation, master-TF x miRNA condition
p <- ggplot(res, aes(x = mirna, y = th_label, fill = pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f", pct)), size = 2.9, color = "black") +
  facet_wrap(~stimulation, nrow = 1) +
  scale_fill_gradientn(colors = c("#F3F4EF", "#8AA6C7", "#1F6F63"), limits = c(0, 100), name = "% runs\nactive") +
  labs(
    title = "Th master-regulator activation across stimulation x miRNA conditions (5x4 grid, 300 runs each)",
    subtitle = "Classic polarising cytokines forced as fixed exogenous inputs (IL12/IL4/IL6+TGFB/TGFB), on top of TCR+IL2 baseline",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(face = "bold"))

ggsave("figures/stimulation_grid_heatmap.png", p, width = 15, height = 6, dpi = 170, bg = "white")
cat("\nSaved figures/stimulation_grid_heatmap.png\n")

# ---- validation check: does the expected master TF dominate under its own polarising condition (no miRNA)?
expect <- tibble::tribble(
  ~stimulation, ~expected_tf, ~expected_label,
  "Th1 (+IL12)", "TBX21", "Th1",
  "Th2 (+IL4)", "GATA3", "Th2",
  "Th17 (+IL6+TGFB)", "RORC", "Th17",
  "Treg (+TGFB)", "FOXP3", "iTreg"
)
check <- res %>% filter(mirna == "no miRNA") %>%
  inner_join(expect, by = c("stimulation" = "stimulation", "node" = "expected_tf")) %>%
  select(stimulation, expected_label, node, pct, ci_lo, ci_hi)
cat("\n=== Validation: does the expected master TF respond under its own polarising cytokine (no miRNA)? ===\n")
print(check)
write_tsv(check, "data/stimulation_grid_validation.tsv")
