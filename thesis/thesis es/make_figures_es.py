#!/usr/bin/env python3
"""Figuras 4-6 de la tesis (umbral, abundancia, simulaciones) con los textos en español.
Las Figuras 1-3 las genera thesis/make_figures_r5.py.

Redibuja las figuras desde las tablas de resultados ya guardadas (no recalcula nada),
con el mismo diseño que las figuras en inglés y coma decimal.

    python3 "thesis/thesis es/make_figures_es.py"

Lee:  analysis/kleiner_benchmark/results/unfiltered/{per_species,threshold_sweep}.tsv
      analysis/simulations/results/study{A,B,C}_*_per_community.tsv
Escribe: thesis/thesis es/figuras/figura1.png ... figura5.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
KLEINER = ROOT / "analysis" / "kleiner_benchmark" / "results" / "unfiltered"
SIMS = ROOT / "analysis" / "simulations" / "results"
OUT = HERE / "figuras"

CONDITIONS = ("Constant", "Protein", "Uneven")          # nombres de los diseños: quedan en inglés
ARMS = ("Unique", "Rigid LCA", "Proportional", "EM", "EM + PEP")
ES = {"Unique": "Únicos", "Rigid LCA": "LCA rígido", "Proportional": "Proporcional",
      "EM": "EM", "EM + PEP": "EM + PEP", "EM, no threshold": "EM, sin umbral"}
INK, QUIET, GRID = "#222222", "#666666", "#cccccc"


def coma(x, _pos=None):
    """Número con coma decimal, sin ceros de más."""
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return s.replace(".", ",")


COMA = FuncFormatter(coma)


# ----------------------------------------------------------------- Figura 1: ejemplo simplificado
def figura1(path):
    col = {"A": "#4e79a7", "B": "#76b7b2", "Genus": "#b07aa1", "C": "#59a14f", "D": "#e15759"}
    fig = plt.figure(figsize=(12, 6.2))
    ax = fig.add_axes([0.02, 0.05, 0.36, 0.9]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.text(0, 9.6, "Comunidad y evidencia peptídica", fontsize=12, weight="bold", color=INK)
    nodes = {"root": (5, 8.4), "Genus": (3, 7.0), "A": (2, 5.6), "B": (4, 5.6), "C": (6.5, 5.6), "D": (8.5, 5.6)}
    for a, b in [("root", "Genus"), ("Genus", "A"), ("Genus", "B"), ("root", "C"), ("root", "D")]:
        (x1, y1), (x2, y2) = nodes[a], nodes[b]
        ax.plot([x1, x2], [y1, y2], color=QUIET, lw=1.2)
    for n, (x, y) in nodes.items():
        c = col.get(n, "white")
        ax.add_patch(plt.Circle((x, y), 0.32, color=c if n != "root" else "white", ec=QUIET, zorder=3))
        lab = {"Genus": "género\n(A,B)", "D": "D\n(ausente)", "root": ""}.get(n, n)
        top = n in ("root", "Genus")
        ax.text(x, y + 0.5 if top else y - 0.75, lab, ha="center", va="bottom" if top else "top",
                fontsize=9.5, color=INK)
    rows = [("pep1", "solo A", "60", "0,001"), ("pep2", "solo B", "20", "0,001"),
            ("pep3", "A y B (compartido)", "40", "0,001"), ("pep4", "solo C", "30", "0,001"),
            ("pep5", "solo C", "10", "0,3"), ("pep6", "solo D (ID falsa)", "20", "0,9")]
    ax.text(0.2, 4.0, "péptido   mapea a                 intensidad  PEP", fontsize=9.5, color=QUIET, family="monospace")
    for i, (p, m, inten, e) in enumerate(rows):
        ax.text(0.2, 3.4 - i * 0.55, f"{p:<9} {m:<23} {inten:>6}     {e}", fontsize=9.5, color=INK, family="monospace")
    ax.text(0.2, -0.2, "Las estrategias 1-4 conservan solo PEP < 0,01 (pep5 y pep6 se descartan).\n"
                       "La estrategia 5 conserva todos, cada uno ponderado por 1 - PEP.", fontsize=9, color=QUIET)

    arms = {
        "1  Únicos":       {"A": 60, "B": 20, "Genus": 0, "C": 30, "D": 0},
        "2  LCA rígido":   {"A": 60, "B": 20, "Genus": 40, "C": 30, "D": 0},
        "3  Proporcional": {"A": 90, "B": 30, "Genus": 0, "C": 30, "D": 0},
        "4  EM":           {"A": 90, "B": 30, "Genus": 0, "C": 30, "D": 0},
        "5  EM + PEP":     {"A": 90, "B": 30, "Genus": 0, "C": 37, "D": 2},
    }
    notes = {"1  Únicos": "pep3 (compartido) se descarta", "2  LCA rígido": "pep3 queda en el nodo del género",
             "3  Proporcional": "pep3 se reparte 3:1, como la evidencia única",
             "4  EM": "mismo resultado aquí (poco compartido)",
             "5  EM + PEP": "C gana 7 por pep5; D recibe 2 por pep6"}
    show = {"Genus": "Género"}
    bx = fig.add_axes([0.45, 0.12, 0.53, 0.78])
    for i, (arm, v) in enumerate(arms.items()):
        y = len(arms) - 1 - i; left = 0
        for c in ["A", "B", "Genus", "C", "D"]:
            if v[c]:
                bx.barh(y, v[c], left=left, color=col[c], height=0.55, edgecolor="white")
                big = v[c] >= 14
                bx.text(left + v[c] / 2, y, f"{show.get(c, c)}\n{v[c]}" if big else str(v[c]),
                        ha="center", va="center", fontsize=8.5, color="white" if big else INK)
                left += v[c]
        bx.text(-3, y, arm, ha="right", va="center", fontsize=10.5, weight="bold", color=INK)
        bx.text(left + 3, y, notes[arm], ha="left", va="center", fontsize=9, color=QUIET)
    bx.set_xlim(0, 260); bx.set_ylim(-0.6, len(arms) - 0.4); bx.axis("off")
    bx.set_title("Dónde pone la intensidad cada estrategia", fontsize=12, weight="bold", color=INK, loc="left", x=-0.25)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- Figura 2: error por cepa
def figura2(ps, path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    colours = ["#8c8c8c", "#b07aa1", "#4e79a7", "#59a14f", "#e15759"]
    for ax, c in zip(axes, CONDITIONS):
        p = ps[(ps.condition == c) & ps.present]
        wide = p.pivot(index="strain", columns="arm", values="xfold")[list(ARMS)]
        for _, row in wide.iterrows():
            ax.plot(range(len(ARMS)), row.values, color="#d0d0d0", lw=0.6, zorder=1)
        for j, arm in enumerate(ARMS):
            ax.scatter(np.full(len(wide), j), wide[arm], s=12, color=colours[j], zorder=2)
            ax.hlines(wide[arm].median(), j - 0.3, j + 0.3, color="black", lw=2, zorder=3)
        ax.set_yscale("log")
        ax.set_xticks(range(len(ARMS)), [ES[a] for a in ARMS], rotation=30, ha="right")
        ax.set_title(f"{c} (n = {len(wide)} cepas)")
        ax.axhline(1, color="black", lw=0.5, ls=":")
    axes[0].set_ylabel("desvío en veces respecto de lo agregado\n(1 = perfecto)")
    fig.suptitle("Error por cepa, pareado entre estrategias (barra negra = mediana)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------- Figura 3: barrido de umbrales
def figura3(sw, path):
    sweep = (0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.01)  # 1.01 = sin umbral
    metrics = [("median_xfold", "mediana del\ndesvío en veces"), ("L1", "distancia L1"),
               ("absent_leakage_pct", "señal en cepas\nausentes (%)")]
    fig, axes = plt.subplots(len(metrics), 3, figsize=(12, 8), sharex=True)
    xs = list(range(len(sweep)))
    labels_x = [coma(t) if t <= 1 else "ninguno" for t in sweep]
    for j, c in enumerate(CONDITIONS):
        for i, (col, name) in enumerate(metrics):
            ax = axes[i, j]
            for w, colour, lab in (("none", "#59a14f", "EM, sin ponderar"), ("1 - PEP", "#e15759", "EM + PEP")):
                q = sw[(sw.condition == c) & (sw.weighting == w)].copy()
                q["t"] = q.threshold.replace("none", "1.01").astype(float)
                v = q.sort_values("t")[col].to_numpy()
                ax.plot(xs, v, marker="o", color=colour, label=lab)
            ax.yaxis.set_major_formatter(COMA)
            if i == 0:
                ax.set_title(c)
            if j == 0:
                ax.set_ylabel(name)
            if i == len(metrics) - 1:
                ax.set_xticks(xs, labels_x, rotation=45)
    axes[0, 0].legend(frameon=False)
    fig.supxlabel("umbral de PEP (se conservan los PSM con PEP < umbral)")
    fig.suptitle("¿Depende el resultado del umbral? (línea plana = no)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------- Figura 4: EM vs EM + PEP por abundancia
def figura4(ps, path):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
    for ax, c in zip(axes, CONDITIONS):
        q = ps[(ps.condition == c) & ps.present]
        a = q[q.arm == "EM"].set_index("strain")
        b = q[q.arm == "EM + PEP"].set_index("strain").loc[a.index]
        rare = a.expected < 0.01
        lim = (0.9, max(a.xfold.max(), b.xfold.max()) * 1.6)
        ax.plot(lim, lim, color=GRID, lw=1, zorder=1)
        ax.fill_between(lim, [lim[0]] * 2, lim, color="#e15759", alpha=0.05, zorder=0)
        ax.scatter(a.xfold[~rare], b.xfold[~rare], s=28, color="#4e79a7", label="abundantes (≥ 1 %)", zorder=3)
        ax.scatter(a.xfold[rare], b.xfold[rare], s=28, color="#f28e2b", label="raras (< 1 %)", zorder=3)
        for s in a.index:
            r = a.xfold[s] / b.xfold[s]
            if r > 1.8 or r < 1 / 1.8:
                ax.annotate(s, (a.xfold[s], b.xfold[s]), xytext=(4, -10), textcoords="offset points",
                            fontsize=8.5, color=INK)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
        better = int((b.log10_error < a.log10_error).sum())
        ax.set_title(f"{c}: {better}/{len(a)} cepas mejoran con PEP", fontsize=11)
        ax.set_xlabel("desvío en veces, EM con umbral estricto")
        ax.text(0.97, 0.04, "debajo de la línea = mejor con PEP", transform=ax.transAxes, ha="right",
                fontsize=8.5, color=QUIET)
    axes[0].set_ylabel("desvío en veces, EM + PEP (sin umbral)")
    axes[0].legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------- Figura 5: simulaciones
def figura5(dfa, dfb, dfc, path, n_communities=20):
    colours = {"Unique": "#8c8c8c", "Proportional": "#4e79a7", "EM": "#59a14f",
               "EM + PEP": "#e15759", "EM, no threshold": "#f28e2b"}
    styles = {"EM, no threshold": ":"}
    shown = ["Unique", "Proportional", "EM", "EM + PEP", "EM, no threshold"]
    panels = [
        (dfa, "false_rate", "Identificaciones falsas (%)", "L1", "A  Más errores de identificación", 100),
        (dfa, "false_rate", "Identificaciones falsas (%)", "absent_leakage_pct", "A  Señal en especies ausentes", 100),
        (dfb, "unique_pct", "Péptidos únicos (%)", "L1", "B  Más péptidos compartidos", 1),
        (dfc, "false_rate", "Identificaciones falsas (%)", "rare_median_xfold",
         "C  La confianza depende de la intensidad:\n    especies por debajo del 1 %", 100),
    ]
    ylabels = {"L1": "distancia L1 (menor = mejor)", "absent_leakage_pct": "señal en especies ausentes (%)",
               "rare_median_xfold": "mediana del desvío en veces\n(1 = perfecto)"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, (df, x, xlab, metric, title, scale) in zip(axes.flat, panels):
        if x == "unique_pct":
            agg = (df.groupby(["within_genus_share", "arm"])
                     .agg(x=(x, "mean"), m=(metric, "mean"), s=(metric, "std")).reset_index())
        else:
            agg = df.groupby([x, "arm"]).agg(m=(metric, "mean"), s=(metric, "std")).reset_index()
            agg["x"] = agg[x] * scale
        for arm in shown:
            g = agg[agg.arm == arm].sort_values("x")
            ax.errorbar(g.x, g.m, yerr=g.s / np.sqrt(n_communities), marker="o", ms=4, capsize=2,
                        color=colours[arm], ls=styles.get(arm, "-"), label=ES[arm])
        if x == "unique_pct":
            ax.invert_xaxis()
        ax.xaxis.set_major_formatter(COMA); ax.yaxis.set_major_formatter(COMA)
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylabels[metric])
    axes[0, 0].legend(frameon=False, fontsize=8.5)
    fig.suptitle(f"Comunidades simuladas (30 especies, {n_communities} por punto; barras = error estándar)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    ps = pd.read_csv(KLEINER / "per_species.tsv", sep="\t")
    sw = pd.read_csv(KLEINER / "threshold_sweep.tsv", sep="\t", dtype={"threshold": str})
    r = lambda n: pd.read_csv(SIMS / f"study{n}_per_community.tsv", sep="\t")
    figura3(sw, OUT / "umbral.png")
    figura4(ps, OUT / "abundancia.png")
    figura5(r("A_false_rate"), r("B_sharing"), r("C_intensity_scores"), OUT / "simulaciones.png")
    print(f"Figuras escritas en {OUT}")


if __name__ == "__main__":
    main()
