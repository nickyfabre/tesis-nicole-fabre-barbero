#!/usr/bin/env python3
"""Figure 1 (fusion of the thesis toy example and the poster layout), English and Spanish.

    python3 make_fig1.py <outdir>

Left: the toy community (tree with a genus A+B, an unrelated strain C and an absent strain D)
and its six peptide identifications, drawn as boxes wired to the strains that contain them.
Right: one card per arm, with a mini-tree showing where the shared signal goes and bars with
the intensity each strain receives. Numbers checked with taxoshare.allocate_species.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

INK, QUIET, GRID = "#222222", "#666666", "#cccccc"
COL = {"A": "#4e79a7", "B": "#76b7b2", "Genus": "#b07aa1", "C": "#59a14f", "D": "#e15759"}
SHARED = "#b07aa1"      # part of a bar that comes from the shared peptide pep3
PEP5 = "#a6d49f"        # part of C that comes from the low-confidence pep5
LOWCONF = "#e15759"

T = {
    "en": dict(
        input="Input: community and identifications", output="Output of each arm",
        absent="D\n(absent)", genus="genus", shared="shared", lowconf="low confidence",
        foot=("Number = MS1 intensity.  Dashed boxes = low-confidence identifications (PEP ≥ 0.01):\n"
              "removed by the threshold in arms 1–4, weighted by 1 − PEP in arm 5."),
        arms=["1  Unique", "2  Rigid LCA", "3  Proportional", "4  EM", "5  EM + PEP"],
        desc=["discards pep3",
              "keeps pep3 at the genus node",
              "splits pep3 once, 3:1,\nanchored on unique peptides",
              "starts 50:50 and re-splits\nuntil stable (same here)",
              "uses every identification,\nweighted by 1 − PEP"],
        em_iter="20:20 → … → 30:10",
        legend="bar colour: strain's own peptides;  purple = from pep3;  light green = from pep5;  red = from pep6",
        dec=".",
    ),
    "es": dict(
        input="Entrada: comunidad e identificaciones", output="Resultado de cada estrategia",
        absent="D\n(ausente)", genus="género", shared="compartido", lowconf="baja confianza",
        foot=("Número = intensidad MS1.  Cajas punteadas = identificaciones de baja confianza (PEP ≥ 0,01):\n"
              "el umbral las elimina en las estrategias 1–4; en la 5 se ponderan por 1 − PEP."),
        arms=["1  Únicos", "2  LCA rígido", "3  Proporcional", "4  EM", "5  EM + PEP"],
        desc=["descarta pep3",
              "deja pep3 en el nodo\ndel género",
              "reparte pep3 una vez, 3:1,\nanclado en los únicos",
              "parte de 50:50 y repite\nhasta estabilizarse (igual aquí)",
              "usa todas las identificaciones,\nponderadas por 1 − PEP"],
        em_iter="20:20 → … → 30:10",
        legend="color de la barra: péptidos propios de la cepa;  violeta = de pep3;  verde claro = de pep5;  rojo = de pep6",
        dec=",",
    ),
}

# (peptide, intensity, PEP, strains, x position of the box)
PEPTIDES = [("pep1", 60, 0.001, "A", 1.0), ("pep3", 40, 0.001, "AB", 2.55), ("pep2", 20, 0.001, "B", 4.1),
            ("pep4", 30, 0.001, "C", 5.6), ("pep5", 10, 0.3, "C", 7.1), ("pep6", 20, 0.9, "D", 8.7)]

# what each arm assigns: own (strain's unique peptides), shared (from pep3), extra (pep5 / pep6)
OUT = [
    dict(own={"A": 60, "B": 20, "C": 30}, shared={}, extra={}, genus=0),
    dict(own={"A": 60, "B": 20, "C": 30}, shared={}, extra={}, genus=40),
    dict(own={"A": 60, "B": 20, "C": 30}, shared={"A": 30, "B": 10}, extra={}, genus=0),
    dict(own={"A": 60, "B": 20, "C": 30}, shared={"A": 30, "B": 10}, extra={}, genus=0),
    dict(own={"A": 60, "B": 20, "C": 30}, shared={"A": 30, "B": 10}, extra={"C": 7, "D": 2}, genus=0),
]


def num(x, dec):
    s = f"{x:g}"
    return s.replace(".", ",") if dec == "," else s


def box(ax, x, y, w, h, fc, ec, ls="-", lw=1.4):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.15",
                                fc=fc, ec=ec, ls=ls, lw=lw, zorder=3))


def leaf(ax, x, y, name, r, absent=False, label=True):
    ax.add_patch(Circle((x, y), r, fc="white" if absent else COL[name], ec=COL[name] if absent else "white",
                        ls="--" if absent else "-", lw=1.6, zorder=4))
    if label:
        ax.text(x, y, name, ha="center", va="center", fontsize=10 if r > 0.3 else 8.5, weight="bold",
                color=COL[name] if absent else "white", zorder=5)


def input_panel(ax, t):
    ax.text(0.1, 12.55, t["input"], fontsize=13, weight="bold", color=INK)
    root, gen = (5.0, 11.0), (2.55, 9.9)
    lv = {"A": (1.6, 8.6), "B": (3.5, 8.6), "C": (6.3, 8.6), "D": (8.6, 8.6)}
    for a, b in [(root, gen), (gen, lv["A"]), (gen, lv["B"]), (root, lv["C"]), (root, lv["D"])]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=GRID, lw=2.2, zorder=1)
    ax.add_patch(Circle(root, 0.18, fc=GRID, ec="none", zorder=2))
    ax.add_patch(Circle(gen, 0.22, fc=COL["Genus"], ec="none", zorder=2))
    ax.text(gen[0] - 0.35, gen[1] + 0.1, t["genus"], ha="right", va="center", fontsize=9.5, color=COL["Genus"])
    for n, (x, y) in lv.items():
        leaf(ax, x, y, n, 0.38, absent=(n == "D"))
    ax.text(lv["D"][0] + 0.5, lv["D"][1], t["absent"].split("\n")[1], ha="left", va="center",
            fontsize=9, color=QUIET)
    by = 5.6
    for name, inten, pep, strains, x in PEPTIDES:
        low = pep >= 0.01
        sh = len(strains) > 1
        ec = LOWCONF if low else (COL["Genus"] if sh else "#9aa5b1")
        fc = "#f3eefa" if sh else "white"
        for s in strains:
            sx, sy = lv[s]
            ax.plot([x, sx], [by + 0.75, sy - 0.38], color=ec, lw=1.6, ls="--" if low else "-", zorder=2)
        box(ax, x, by, 1.25, 1.5, fc, ec, ls="--" if low else "-")
        ax.text(x, by + 0.42, name, ha="center", va="center", fontsize=9.5, weight="bold", color=INK, zorder=5)
        ax.text(x, by - 0.02, str(inten), ha="center", va="center", fontsize=11, color=INK, zorder=5)
        ax.text(x, by - 0.45, f"PEP {num(pep, t['dec'])}", ha="center", va="center", fontsize=8, color=QUIET,
                zorder=5)
    ax.text(2.55, by - 1.15, t["shared"], ha="center", va="top", fontsize=9, color=COL["Genus"], weight="bold")
    ax.text(7.9, by - 1.15, t["lowconf"], ha="center", va="top", fontsize=9, color=LOWCONF, weight="bold")
    ax.text(0.1, 2.9, t["foot"], fontsize=8.6, color=QUIET, va="top", linespacing=1.5)


def card(ax, x0, i, t):
    w = 4.0
    ax.add_patch(FancyBboxPatch((x0, 0.55), w, 11.1, boxstyle="round,pad=0,rounding_size=0.2",
                                fc="#f7f8fa" if i < 2 else "#f2f7fd", ec="#cfd8e3", lw=1.3, zorder=0))
    cx = x0 + w / 2
    ax.text(cx, 11.0, t["arms"][i], ha="center", va="center", fontsize=11.5, weight="bold", color=INK)
    ax.text(cx, 10.05, t["desc"][i], ha="center", va="center", fontsize=8.6, color=QUIET, linespacing=1.3)
    # mini-tree
    root, gen = (x0 + 2.25, 8.35), (x0 + 1.15, 7.6)
    lv = {"A": (x0 + 0.65, 6.6), "B": (x0 + 1.65, 6.6), "C": (x0 + 2.7, 6.6), "D": (x0 + 3.5, 6.6)}
    o = OUT[i]
    split = bool(o["shared"])
    for a, b, key in [(root, gen, None), (gen, lv["A"], "A"), (gen, lv["B"], "B"),
                      (root, lv["C"], None), (root, lv["D"], None)]:
        c = SHARED if (split and key) else GRID
        ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=2.4 if c == SHARED else 1.8, zorder=1)
    ax.add_patch(Circle(root, 0.12, fc=GRID, ec="none", zorder=2))
    ax.add_patch(Circle(gen, 0.16, fc=COL["Genus"], ec="none", zorder=2))
    for n, (x, y) in lv.items():
        leaf(ax, x, y, n, 0.24, absent=(n == "D"))
    if i == 0:
        ax.text(gen[0], gen[1] + 0.45, "✗", ha="center", va="center", fontsize=15, color="#b3261e", weight="bold")
    if o["genus"]:
        box(ax, gen[0] - 0.05, gen[1] + 0.62, 0.75, 0.5, COL["Genus"], COL["Genus"])
        ax.text(gen[0] - 0.05, gen[1] + 0.62, str(o["genus"]), ha="center", va="center", fontsize=10,
                weight="bold", color="white", zorder=5)
    if split:
        for s, dx in (("A", -0.42), ("B", 0.42)):
            ax.text((gen[0] + lv[s][0]) / 2 + dx, (gen[1] + lv[s][1]) / 2, str(o["shared"][s]), ha="center",
                    va="center", fontsize=9, weight="bold", color=SHARED)
    if i == 3:
        ax.text(cx, 5.85, t["em_iter"], ha="center", va="center", fontsize=8.6, color=SHARED, weight="bold")
    if i == 4:
        ax.text(lv["C"][0] + 0.02, 5.9, "+7", ha="center", va="center", fontsize=9, weight="bold", color="#3d8b3d")
        ax.text(lv["D"][0], 5.9, "+2", ha="center", va="center", fontsize=9, weight="bold", color=LOWCONF)
    # bars
    base, scale, bw = 1.15, 4.0 / 90, 0.5
    for n, (x, _) in lv.items():
        parts = [(o["own"].get(n, 0), COL[n]), (o["shared"].get(n, 0), SHARED),
                 (o["extra"].get(n, 0), PEP5 if n == "C" else COL["D"])]
        y = base
        for v, c in parts:
            if v:
                ax.add_patch(Rectangle((x - bw / 2, y), bw, v * scale, fc=c, ec="white", lw=0.8, zorder=3))
                y += v * scale
        total = sum(v for v, _ in parts)
        if total:
            ax.text(x, y + 0.12, str(total), ha="center", va="bottom", fontsize=9.5, weight="bold", color=INK)
    ax.plot([x0 + 0.25, x0 + w - 0.25], [base, base], color=GRID, lw=1, zorder=2)


def main(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for lang, name in (("en", "method.png"), ("es", "metodo.png")):
        t = T[lang]
        fig = plt.figure(figsize=(17, 7.2))
        ax = fig.add_axes([0, 0.03, 1, 0.95])
        ax.set_xlim(0, 32.6); ax.set_ylim(0, 13.2); ax.axis("off")
        input_panel(ax, t)
        ax.annotate("", xy=(11.15, 6.4), xytext=(10.35, 6.4),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=2, mutation_scale=18))
        ax.text(11.4 + 2.5 * 4.2 - 0.1, 12.55, t["output"], fontsize=13, weight="bold", color=INK, ha="center")
        for i in range(5):
            card(ax, 11.4 + i * 4.25, i, t)
        ax.text(11.4 + 2.5 * 4.25 - 0.1, 0.12, t["legend"], ha="center", va="center", fontsize=8.6, color=QUIET)
        fig.savefig(outdir / name, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print("written", outdir / name)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
