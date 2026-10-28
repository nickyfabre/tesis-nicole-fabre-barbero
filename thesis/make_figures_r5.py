#!/usr/bin/env python3
"""New thesis figures (round 5), in Spanish and English.

    python3 thesis/make_figures_r5.py

Draws, from the saved result tables (nothing is recomputed):
  - method figure (five arms on a toy example), with the two panels separated;
  - x-fold deviation of the four sharing strains vs the rest, per arm;
  - estimated vs expected abundance of the three clades that carry the shared signal.

Writes thesis/thesis es/figuras/{metodo,cepas_compartidas,clados}.png
       thesis/figures/{method,sharing_strains,clades}.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
KLEINER = ROOT / "analysis" / "kleiner_benchmark" / "results" / "unfiltered"
OUT = {"es": HERE / "thesis es" / "figuras", "en": HERE / "figures"}
NAMES = {"es": ("metodo", "cepas_compartidas", "clados"), "en": ("method", "sharing_strains", "clades")}

CONDITIONS = ("Constant", "Protein", "Uneven")
ARMS = ("Unique", "Rigid LCA", "Proportional", "EM", "EM + PEP")
ARM_LABEL = {"es": {"Unique": "Únicos", "Rigid LCA": "LCA rígido", "Proportional": "Proporcional",
                    "EM": "EM", "EM + PEP": "EM + PEP"},
             "en": {a: a for a in ARMS}}
INK, QUIET, GRID = "#222222", "#666666", "#cccccc"
SHARING = ("841", "VF", "137", "259")
SHARING_COL = {"841": "#4e79a7", "VF": "#76b7b2", "137": "#e15759", "259": "#f28e2b"}

T = {
    "es": dict(
        comm="Comunidad y evidencia peptídica", where="Dónde pone la intensidad cada estrategia",
        header="péptido   mapea a                 intensidad  PEP",
        rows=[("pep1", "solo A", "60", "0,001"), ("pep2", "solo B", "20", "0,001"),
              ("pep3", "A y B (compartido)", "40", "0,001"), ("pep4", "solo C", "30", "0,001"),
              ("pep5", "solo C", "10", "0,3"), ("pep6", "solo D (ID falsa)", "20", "0,9")],
        foot="Las estrategias 1-4 conservan solo PEP < 0,01 (pep5 y pep6 se descartan).\n"
             "La estrategia 5 conserva todos, cada uno ponderado por 1 - PEP.",
        genus="género\n(A,B)", absent="D\n(ausente)", Genus="Género",
        arms=["1  Únicos", "2  LCA rígido", "3  Proporcional", "4  EM", "5  EM + PEP"],
        notes=["pep3 (compartido) se descarta", "pep3 queda en el nodo del género",
               "pep3 se reparte 3:1, como la evidencia única", "mismo resultado aquí (poco compartido)",
               "C gana 7 por pep5; D recibe 2 por pep6"],
        ylab="desvío en veces respecto de lo agregado\n(1 = perfecto, escala log)",
        rest="resto de las cepas: mediana (línea) y RIC (banda)",
        clade_names={"tax379": "841 + VF\nRhizobium", "tax1280": "137 + 259\nS. aureus",
                     "tax543": "K12 + LT2\nEnterobacteriaceae"},
        expected="esperado", clade_y="abundancia del clado (%)",
        lca="LCA rígido = Proporcional = EM"),
    "en": dict(
        comm="Community and peptide evidence", where="Where each arm puts the intensity",
        header="peptide   maps to                 intensity   PEP",
        rows=[("pep1", "A only", "60", "0.001"), ("pep2", "B only", "20", "0.001"),
              ("pep3", "A and B (shared)", "40", "0.001"), ("pep4", "C only", "30", "0.001"),
              ("pep5", "C only", "10", "0.3"), ("pep6", "D only (false ID)", "20", "0.9")],
        foot="Arms 1-4 keep only PEP < 0.01 (pep5 and pep6 are discarded).\n"
             "Arm 5 keeps all of them, each weighted by 1 - PEP.",
        genus="genus\n(A,B)", absent="D\n(absent)", Genus="Genus",
        arms=["1  Unique", "2  Rigid LCA", "3  Proportional", "4  EM", "5  EM + PEP"],
        notes=["pep3 (shared) is discarded", "pep3 stays at the genus node",
               "pep3 is split 3:1, like the unique evidence", "same result here (little sharing)",
               "C gains 7 from pep5; D receives 2 from pep6"],
        ylab="x-fold deviation from input\n(1 = perfect, log scale)",
        rest="other strains: median (line) and IQR (band)",
        clade_names={"tax379": "841 + VF\nRhizobium", "tax1280": "137 + 259\nS. aureus",
                     "tax543": "K12 + LT2\nEnterobacteriaceae"},
        expected="expected", clade_y="clade abundance (%)",
        lca="Rigid LCA = Proportional = EM"),
}


def fmt(lang):
    def f(x, _pos=None):
        s = f"{x:.3f}".rstrip("0").rstrip(".")
        return s.replace(".", ",") if lang == "es" else s
    return FuncFormatter(f)


def method(lang, path):
    t = T[lang]
    col = {"A": "#4e79a7", "B": "#76b7b2", "Genus": "#b07aa1", "C": "#59a14f", "D": "#e15759"}
    fig = plt.figure(figsize=(12, 6.2))
    ax = fig.add_axes([0.02, 0.05, 0.36, 0.9]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.text(0, 9.6, t["comm"], fontsize=12, weight="bold", color=INK)
    nodes = {"root": (5, 8.4), "Genus": (3, 7.0), "A": (2, 5.6), "B": (4, 5.6), "C": (6.5, 5.6), "D": (8.5, 5.6)}
    for a, b in [("root", "Genus"), ("Genus", "A"), ("Genus", "B"), ("root", "C"), ("root", "D")]:
        (x1, y1), (x2, y2) = nodes[a], nodes[b]
        ax.plot([x1, x2], [y1, y2], color=QUIET, lw=1.2)
    for n, (x, y) in nodes.items():
        ax.add_patch(plt.Circle((x, y), 0.32, color=col.get(n, "white"), ec=QUIET, zorder=3))
        lab = {"Genus": t["genus"], "D": t["absent"], "root": ""}.get(n, n)
        top = n in ("root", "Genus")
        ax.text(x, y + 0.5 if top else y - 0.75, lab, ha="center", va="bottom" if top else "top",
                fontsize=9.5, color=INK)
    ax.text(0.2, 4.0, t["header"], fontsize=9.5, color=QUIET, family="monospace")
    for i, (p, m, inten, e) in enumerate(t["rows"]):
        ax.text(0.2, 3.4 - i * 0.55, f"{p:<9} {m:<23} {inten:>6}     {e}", fontsize=9.5, color=INK,
                family="monospace")
    ax.text(0.2, -0.2, t["foot"], fontsize=9, color=QUIET)

    # vertical separator between the two panels
    fig.add_artist(Line2D([0.405, 0.405], [0.02, 0.98], transform=fig.transFigure, color=GRID, lw=1.5))

    vals = [{"A": 60, "B": 20, "Genus": 0, "C": 30, "D": 0},
            {"A": 60, "B": 20, "Genus": 40, "C": 30, "D": 0},
            {"A": 90, "B": 30, "Genus": 0, "C": 30, "D": 0},
            {"A": 90, "B": 30, "Genus": 0, "C": 30, "D": 0},
            {"A": 90, "B": 30, "Genus": 0, "C": 37, "D": 2}]
    bx = fig.add_axes([0.56, 0.12, 0.43, 0.78])
    for i, (arm, v, note) in enumerate(zip(t["arms"], vals, t["notes"])):
        y = len(vals) - 1 - i; left = 0
        for c in ["A", "B", "Genus", "C", "D"]:
            if v[c]:
                bx.barh(y, v[c], left=left, color=col[c], height=0.55, edgecolor="white")
                big = v[c] >= 14
                name = t["Genus"] if c == "Genus" else c
                bx.text(left + v[c] / 2, y, f"{name}\n{v[c]}" if big else str(v[c]),
                        ha="center", va="center", fontsize=8.5, color="white" if big else INK)
                left += v[c]
        bx.text(-3, y, arm, ha="right", va="center", fontsize=10.5, weight="bold", color=INK)
        bx.text(left + (9 if v["D"] else 3), y, note, ha="left", va="center", fontsize=9, color=QUIET)
    bx.set_xlim(0, 300); bx.set_ylim(-0.6, len(vals) - 0.4); bx.axis("off")
    bx.set_title(t["where"], fontsize=12, weight="bold", color=INK, loc="left", x=-0.3)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def sharing(lang, ps, path):
    t = T[lang]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
    xs = np.arange(len(ARMS))
    for ax, c in zip(axes, CONDITIONS):
        p = ps[(ps.condition == c) & ps.present]
        wide = p.pivot(index="strain", columns="arm", values="xfold")[list(ARMS)]
        rest = wide.drop(index=[s for s in SHARING if s in wide.index])
        q1, med, q3 = rest.quantile(0.25), rest.median(), rest.quantile(0.75)
        ax.fill_between(xs, q1.values, q3.values, color="#bbbbbb", alpha=0.45, lw=0)
        ax.plot(xs, med.values, color="#555555", lw=2.2, marker="o", ms=4)
        for s in SHARING:
            ax.plot(xs, wide.loc[s].values, color=SHARING_COL[s], lw=1.8, marker="o", ms=5, label=s)
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(fmt(lang))
        ax.axhline(1, color="black", lw=0.6, ls=":")
        ax.set_xticks(xs, [ARM_LABEL[lang][a] for a in ARMS], rotation=30, ha="right")
        ax.set_title(c)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel(t["ylab"])
    handles = [Line2D([], [], color=SHARING_COL[s], marker="o", lw=1.8, label=s) for s in SHARING]
    handles.append(Line2D([], [], color="#555555", lw=6, alpha=0.5, label=t["rest"]))
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def clades(lang, ps, nodes, path):
    t = T[lang]
    members = {"tax379": ("841", "VF"), "tax1280": ("137", "259"), "tax543": ("K12", "LT2")}
    show = [("Unique", "#8c8c8c", ARM_LABEL[lang]["Unique"]), ("EM", "#4e79a7", t["lca"])]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    w = 0.26
    for ax, c in zip(axes, CONDITIONS):
        exp = ps[(ps.condition == c) & (ps.arm == "EM")].set_index("strain").expected
        x = np.arange(len(members))
        e = [100 * exp.loc[list(m)].sum() for m in members.values()]
        ax.bar(x - w, e, w, color="white", edgecolor=INK, hatch="//", label=t["expected"])
        for k, (arm, colr, lab) in enumerate(show):
            n = nodes[(nodes.condition == c) & (nodes.arm == arm)].set_index("node").share
            ax.bar(x + k * w, [100 * n.loc[node] for node in members], w, color=colr, label=lab)
        ax.set_xticks(x, [t["clade_names"][n].split("\n")[0] for n in members], fontsize=9)
        for xi, n in zip(x, members):  # taxon names in italics, on a second line
            ax.annotate(t["clade_names"][n].split("\n")[1], (xi, 0), xycoords=("data", "axes fraction"),
                        xytext=(0, -19), textcoords="offset points", ha="center", va="top",
                        fontsize=9, style="italic")
        ax.set_title(c)
        ax.yaxis.set_major_formatter(fmt(lang))
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel(t["clade_y"])
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ps = pd.read_csv(KLEINER / "per_species.tsv", sep="\t", dtype={"strain": str})
    nodes = pd.read_csv(KLEINER / "nodes.tsv", sep="\t", dtype={"node": str})
    for lang in ("es", "en"):
        OUT[lang].mkdir(parents=True, exist_ok=True)
        a, b, c = NAMES[lang]
        method(lang, OUT[lang] / f"{a}.png")
        sharing(lang, ps, OUT[lang] / f"{b}.png")
        clades(lang, ps, nodes, OUT[lang] / f"{c}.png")
    print("done")


if __name__ == "__main__":
    main()
