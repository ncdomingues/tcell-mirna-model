# -*- coding: utf-8 -*-
"""
Reproduces the thesis's core comparison: naive CD4+ T-cell activation/
differentiation under TCR + IL2 stimulation, across four miRNA regulatory
conditions (neither / both / miR-34c-5p only / miR-155-5p only), using the
digitized rule table in rules_data.py and the engine in logic_engine.py.

Outputs (to ../output and ../figures):
  - steady_state_comparison.csv : per-node activation frequency by condition
  - phenotype_summary.csv       : master-TF activation frequency by condition
  - network_graph.png           : the regulatory network, coloured by category
  - trajectories.png            : average active-fraction over time, per TF
  - phenotype_bars.png          : master-TF activation bar chart by condition
"""
import os
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

from rules_data import RULES, ALIASES, SOURCE
from logic_engine import LogicalModel, parse_formula, collect_refs

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "output")
FIG_DIR = os.path.join(HERE, "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

N_SWEEPS = 10
N_RUNS = 400
SEED = 42

MASTER_TFS = {
    "TBX21": "Th1",
    "GATA3": "Th2",
    "RORC": "Th17",
    "FOXP3": "iTreg",
    "BCL6": "Tfh",
    "AHR": "Th22",
    "SPI1-PU1": "Th9",
}

WATCHLIST = [
    "NFAT", "NFKB", "STAT3", "STAT5", "TBX21", "GATA3", "RORC", "FOXP3",
    "BCL6", "AHR", "SPI1-PU1", "MYC", "TP53", "SP1", "MDM2", "IL17", "IL4",
    "IFNG", "IL2", "CDKN1A-p21",
]

# ---------------------------------------------------------------- build models
def build_condition_models():
    conditions = {}
    both = LogicalModel(RULES, ALIASES)
    conditions["both"] = (both, {})

    no_155 = [r for r in RULES if r[0] != "miR-155-5p"]
    conditions["miR-34c-5p only"] = (LogicalModel(no_155, ALIASES), {"miR-155-5p": 0})

    no_34c = [r for r in RULES if r[0] != "miR-34c-5p"]
    conditions["miR-155-5p only"] = (LogicalModel(no_34c, ALIASES), {"miR-34c-5p": 0})

    neither = [r for r in RULES if r[0] not in ("miR-34c-5p", "miR-155-5p")]
    conditions["no miRNA"] = (LogicalModel(neither, ALIASES), {"miR-34c-5p": 0, "miR-155-5p": 0})

    return conditions


def base_initial_state(model):
    state = {n: 0 for n in model.all_nodes}
    # Sustained antigen presentation + co-receptor + coactivator availability
    # in vitro (TCR + IL2 stimulation of naive CD4 T cells).
    for inp in ("APC-Antigen", "CD3", "CD45RA", "CGC", "CREBBP"):
        if inp in state:
            state[inp] = 1
    return state


def run_all():
    conditions = build_condition_models()
    results = {}
    for label, (model, overrides) in conditions.items():
        init = base_initial_state(model)
        init.update(overrides)
        finals, traj = model.ensemble_run(init, n_sweeps=N_SWEEPS, n_runs=N_RUNS, seed=SEED)
        results[label] = {"model": model, "finals": finals, "traj": traj}
    return conditions, results


# ---------------------------------------------------------------- summaries
def activation_frequency_table(results, nodes):
    order = ["no miRNA", "miR-155-5p only", "miR-34c-5p only", "both"]
    data = {}
    for label in order:
        finals = results[label]["finals"]
        n = len(finals)
        data[label] = [
            round(100.0 * sum(1 for f in finals if f.get(node, 0) >= 1) / n, 1)
            for node in nodes
        ]
    df = pd.DataFrame(data, index=nodes)
    return df[order]


def phenotype_table(results):
    return activation_frequency_table(results, list(MASTER_TFS.keys())).rename(
        index=MASTER_TFS
    )


# ---------------------------------------------------------------- figures
CATEGORY_COLORS = {
    "input": "#7C8579",
    "proximal_signalling": "#35577A",
    "transcription_factor": "#1F6F63",
    "cytokine": "#A3763B",
    "receptor": "#8AA6C7",
    "mirna": "#B23B5E",
    "other": "#9AA08F",
}

TF_NODES = {"NFAT", "NFKB", "TP53", "MDM2", "SP1", "MYC", "STAT1", "STAT3",
            "STAT4", "STAT5", "STAT6", "TBX21", "GATA3", "RORC", "FOXP3",
            "BCL6", "SPI1-PU1", "MAF", "AHR", "PRDM1-BLIMP1", "CREB",
            "FOXO1", "FOXO3", "SIRT1", "CTNNB1"}
CYTOKINE_NODES = {"IL2", "IL4", "IL6", "IL9", "IL10", "IL12", "IL17", "IL21",
                   "IL22", "IL23", "IFNG", "TGFB", "TNF"}
RECEPTOR_NODES = {"IL2R", "IL4R", "IL6R", "IL10R", "IL12R", "IL21R", "IL23R",
                   "IFNGR", "TGFBR", "TNFR", "IL1R"}
MIRNA_NODES = {"miR-34c-5p", "miR-155-5p"}
INPUT_NODES = {"APC-Antigen", "CD3", "CD45RA", "CGC", "CREBBP"}


def categorize(node):
    if node in INPUT_NODES:
        return "input"
    if node in MIRNA_NODES:
        return "mirna"
    if node in TF_NODES:
        return "transcription_factor"
    if node in CYTOKINE_NODES:
        return "cytokine"
    if node in RECEPTOR_NODES:
        return "receptor"
    return "proximal_signalling"


def plot_network(model, path):
    G = nx.DiGraph()
    for node in model.all_nodes:
        G.add_node(node, category=categorize(node))
    for node, value, formula, _desc in RULES:
        ast = parse_formula(formula)
        refs = set()
        collect_refs(ast, refs)
        for r in refs:
            r = ALIASES.get(r, r)
            G.add_edge(r, node)

    pos = nx.spring_layout(G, k=0.55, seed=7, iterations=200)
    fig, ax = plt.subplots(figsize=(13, 11))
    fig.patch.set_facecolor("#F3F4EF")
    ax.set_facecolor("#F3F4EF")

    colors = [CATEGORY_COLORS[categorize(n)] for n in G.nodes()]
    sizes = [900 if n in MIRNA_NODES else (420 if n in TF_NODES else 220) for n in G.nodes()]

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#CBD0C1", width=0.6,
                            arrows=True, arrowsize=6, alpha=0.55, connectionstyle="arc3,rad=0.03")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=colors, node_size=sizes,
                            linewidths=0.6, edgecolors="#1B2321", alpha=0.95)
    labels = {n: n for n in G.nodes() if n in TF_NODES | MIRNA_NODES | CYTOKINE_NODES}
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax, font_size=8,
                             font_family="DejaVu Sans", font_color="#1B2321")

    handles = [plt.Line2D([0], [0], marker='o', color='w', label=k.replace("_", " "),
                           markerfacecolor=v, markersize=10) for k, v in CATEGORY_COLORS.items()]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=9)
    ax.set_title("Digitized regulatory network: CD4+ T-cell activation & differentiation\n"
                  "extended with miR-34c-5p / miR-155-5p (rebuilt from PhD thesis Table 5)",
                  fontsize=12, color="#1B2321", loc="left")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_trajectories(results, path):
    tfs = ["TBX21", "GATA3", "RORC", "FOXP3", "STAT3", "NFKB"]
    conditions = ["no miRNA", "both"]
    colors = {"no miRNA": "#7C8579", "both": "#1F6F63"}

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    fig.patch.set_facecolor("#F3F4EF")
    for ax, tf in zip(axes.flat, tfs):
        ax.set_facecolor("#F3F4EF")
        for cond in conditions:
            traj = results[cond]["traj"][tf]
            ax.plot(range(len(traj)), traj, label=cond, color=colors[cond], linewidth=2)
        ax.set_title(tf, fontsize=11, color="#1B2321")
        ax.set_ylim(-0.03, 1.03)
        ax.grid(alpha=0.25, color="#CBD0C1")
        for spine in ax.spines.values():
            spine.set_color("#CBD0C1")
    axes[0, 0].set_ylabel("average active fraction")
    axes[1, 0].set_ylabel("average active fraction")
    for ax in axes[1]:
        ax.set_xlabel("asynchronous sweep")
    axes[0, -1].legend(loc="upper right", frameon=False, fontsize=9)
    fig.suptitle("Node activity over time, ensemble of 400 asynchronous runs\n"
                  "(TCR + IL2 stimulation, naive CD4+ T cell)", fontsize=12, color="#1B2321")
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_phenotype_bars(pheno_df, path):
    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor("#F3F4EF")
    ax.set_facecolor("#F3F4EF")
    x = np.arange(len(pheno_df.index))
    width = 0.19
    colors = ["#7C8579", "#8AA6C7", "#B23B5E", "#1F6F63"]
    for i, cond in enumerate(pheno_df.columns):
        ax.bar(x + (i - 1.5) * width, pheno_df[cond], width, label=cond, color=colors[i])
    ax.set_xticks(x)
    ax.set_xticklabels(pheno_df.index, fontsize=10)
    ax.set_ylabel("% of runs with master TF active")
    ax.set_title("Th-subtype master-regulator activation by miRNA condition", fontsize=12, color="#1B2321", loc="left")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=0.25, color="#CBD0C1")
    for spine in ax.spines.values():
        spine.set_color("#CBD0C1")
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    conditions, results = run_all()

    ss_table = activation_frequency_table(results, WATCHLIST)
    ss_table.to_csv(os.path.join(OUT_DIR, "steady_state_comparison.csv"))
    print("\n=== Node activation frequency (%), by condition ===")
    print(ss_table)

    pheno_df = phenotype_table(results)
    pheno_df.to_csv(os.path.join(OUT_DIR, "phenotype_summary.csv"))
    print("\n=== Master-TF activation frequency (%), by condition ===")
    print(pheno_df)

    both_model = conditions["both"][0]
    plot_network(both_model, os.path.join(FIG_DIR, "network_graph.png"))
    plot_trajectories(results, os.path.join(FIG_DIR, "trajectories.png"))
    plot_phenotype_bars(pheno_df, os.path.join(FIG_DIR, "phenotype_bars.png"))

    meta = {
        "source": SOURCE,
        "n_sweeps": N_SWEEPS,
        "n_runs": N_RUNS,
        "seed": SEED,
        "rule_nodes": both_model.rule_nodes,
        "input_nodes": both_model.input_nodes,
    }
    with open(os.path.join(OUT_DIR, "run_metadata.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print("\nFigures written to:", FIG_DIR)
    print("Tables written to:", OUT_DIR)


if __name__ == "__main__":
    main()
