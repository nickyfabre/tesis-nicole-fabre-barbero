#!/usr/bin/env python3
"""Figure 1 (fusion of the thesis toy example and the poster layout), English and Spanish.

    python3 make_fig1.py <outdir>

Left: the toy community (tree with two strains A+B of one species, an unrelated strain C and an absent strain D)
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
        absent="D\n(absent)", genus="species", shared="shared", lowconf="low confidence",
        foot=("Number = MS1 intensity.  Dashed boxes = low-confidence identifications (PEP ≥ 0.01):\n"
              "removed by the threshold in arms 1–4, weighted by 1 − PEP in arm 5."),
        arms=["1  Unique", "2  Rigid LCA", "3  Proportional", "4  EM", "5  EM + PEP"],
        desc=["discards pep3",
              "keeps pep3 at the species node",
              "splits pep3 once, 3:1,\nanchored on unique peptides",
              "starts 50:50 and re-splits\nuntil stable (same here)",
              "uses every identification,\nweighted by 1 − PEP"],
        em_iter="20:20 → … → 30:10",
        legend="bar colour: strain's own peptides;  purple = from pep3;  light green = from pep5;  red = from pep6",
        dec=".",
    ),
    "es": dict(
        input="Entrada: comunidad e identificaciones", output="Resultado de cada estrategia",
        absent="D\n(ausente)", genus="especie", shared="compartido", lowconf="baja confianza",
        foot=("Número = intensidad MS1.  Cajas punteadas = identificaciones de baja confianza (PEP ≥ 0,01):\n"
              "el umbral las elimina en las estrategias 1–4; en la 5 se ponderan por 1 − PEP."),
        arms=["1  Únicos", "2  LCA rígido", "3  Proporcional", "4  EM", "5  EM + PEP"],
        desc=["descarta pep3",
              "deja pep3 en el nodo\nde la especie",
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
        ax.text(x, y, name, ha="center", va="center", fontsize=15 if r > 0.4 else 11.5, weight="bold",
                color=COL[name] if absent else "white", zorder=5)



SIDE = {
    "en": ["A and B are strains of one species; C is unrelated;",
           "D is absent from the sample.",
           "",
           "Each box is one identification, with its",
           "MS1 intensity and its PEP.",
           "",
           "pep3 is shared by A and B.",
           "",
           "pep5 and pep6 have low confidence (PEP ≥ 0.01):",
           "arms 1–4 remove them; arm 5 weights them",
           "by 1 − PEP. pep6 is a false identification."],
    "es": ["A y B son cepas de una misma especie; C no está",
           "emparentada; D no está en la muestra.",
           "",
           "Cada caja es una identificación, con su",
           "intensidad MS1 y su PEP.",
           "",
           "pep3 es compartido por A y B.",
           "",
           "pep5 y pep6 son de baja confianza (PEP ≥ 0,01):",
           "las estrategias 1–4 los eliminan; la 5 los",
           "pondera por 1 − PEP. pep6 es una identificación falsa."],
}
DESC = {
    "en": ["discards pep3", "keeps pep3 at\nthe species node", "splits pep3 once,\nanchored on unique",
           "starts 50:50, re-splits\nuntil stable", "all identifications,\nweighted by 1 − PEP"],
    "es": ["descarta pep3", "deja pep3 en el\nnodo de la especie", "reparte pep3 una vez,\nanclado en los únicos",
           "parte de 50:50 y repite\nhasta estabilizarse", "todas las identificaciones,\nponderadas por 1 − PEP"],
}
PEPX = {"pep1": 1.3, "pep3": 3.45, "pep2": 5.6, "pep4": 7.6, "pep5": 9.55, "pep6": 11.5}


def top_panel(ax, t, lang):
    ax.text(0.1, 9.55, t["input"], fontsize=17, weight="bold", color=INK, va="center")
    root, gen = (6.6, 8.6), (3.45, 7.55)
    lv = {"A": (2.3, 6.2), "B": (4.6, 6.2), "C": (8.6, 6.2), "D": (11.5, 6.2)}
    for a, b in [(root, gen), (gen, lv["A"]), (gen, lv["B"]), (root, lv["C"]), (root, lv["D"])]:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=GRID, lw=2.6, zorder=1)
    ax.add_patch(Circle(root, 0.2, fc=GRID, ec="none", zorder=2))
    ax.add_patch(Circle(gen, 0.26, fc=COL["Genus"], ec="none", zorder=2))
    ax.text(gen[0] - 0.4, gen[1] + 0.1, t["genus"], ha="right", va="center", fontsize=14, color=COL["Genus"])
    for n, (x, y) in lv.items():
        leaf(ax, x, y, n, 0.45, absent=(n == "D"))
    by = 3.1
    for name, inten, pep, strains, _ in PEPTIDES:
        x = PEPX[name]
        low, sh = pep >= 0.01, len(strains) > 1
        ec = LOWCONF if low else (COL["Genus"] if sh else "#9aa5b1")
        for s in strains:
            sx, sy = lv[s]
            ax.plot([x, sx], [by + 1.05, sy - 0.45], color=ec, lw=2, ls="--" if low else "-", zorder=2)
        box(ax, x, by, 1.7, 2.1, "#f3eefa" if sh else "white", ec, ls="--" if low else "-", lw=1.8)
        ax.text(x, by + 0.6, name, ha="center", va="center", fontsize=14, weight="bold", color=INK, zorder=5)
        ax.text(x, by, str(inten), ha="center", va="center", fontsize=16, color=INK, zorder=5)
        ax.text(x, by - 0.62, f"PEP {num(pep, t['dec'])}", ha="center", va="center", fontsize=11.5,
                color=QUIET, zorder=5)
    ax.text(PEPX["pep3"], by - 1.4, t["shared"], ha="center", va="top", fontsize=13, color=COL["Genus"], weight="bold")
    ax.text((PEPX["pep5"] + PEPX["pep6"]) / 2, by - 1.4, t["lowconf"], ha="center", va="top", fontsize=13,
            color=LOWCONF, weight="bold")
    ax.add_patch(FancyBboxPatch((13.1, 1.3), 10.7, 7.6, boxstyle="round,pad=0,rounding_size=0.25",
                                fc="#f7f8fa", ec="#cfd8e3", lw=1.3, zorder=0))
    for k, line in enumerate(SIDE[lang]):
        ax.text(13.5, 8.35 - k * 0.62, line, fontsize=13, color=INK if line[:3] in ("pep", "Eac", "Cad", "A a", "A y") or k in (1,) else INK,
                va="center")


def card(ax, x0, i, t, lang):
    w = 4.6
    ax.add_patch(FancyBboxPatch((x0, 0.9), w, 10.0, boxstyle="round,pad=0,rounding_size=0.25",
                                fc="#f7f8fa" if i < 2 else "#f2f7fd", ec="#cfd8e3", lw=1.3, zorder=0))
    cx = x0 + w / 2
    ax.text(cx, 10.25, t["arms"][i], ha="center", va="center", fontsize=16, weight="bold", color=INK)
    ax.text(cx, 9.25, DESC[lang][i], ha="center", va="center", fontsize=12, color=QUIET, linespacing=1.25)
    root, gen = (x0 + 2.6, 8.05), (x0 + 1.3, 7.25)
    lv = {"A": (x0 + 0.7, 6.25), "B": (x0 + 1.9, 6.25), "C": (x0 + 3.1, 6.25), "D": (x0 + 4.0, 6.25)}
    o = OUT[i]
    split = bool(o["shared"])
    for a, b, key in [(root, gen, None), (gen, lv["A"], "A"), (gen, lv["B"], "B"),
                      (root, lv["C"], None), (root, lv["D"], None)]:
        c = SHARED if (split and key) else GRID
        ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=3 if c == SHARED else 2.2, zorder=1)
    ax.add_patch(Circle(root, 0.14, fc=GRID, ec="none", zorder=2))
    ax.add_patch(Circle(gen, 0.19, fc=COL["Genus"], ec="none", zorder=2))
    for n, (x, y) in lv.items():
        leaf(ax, x, y, n, 0.31, absent=(n == "D"))
    if i == 0:
        ax.text(gen[0], gen[1] + 0.5, "✗", ha="center", va="center", fontsize=20, color="#b3261e", weight="bold")
    if o["genus"]:
        box(ax, gen[0] - 0.1, gen[1] + 0.7, 0.95, 0.6, COL["Genus"], COL["Genus"])
        ax.text(gen[0] - 0.1, gen[1] + 0.7, str(o["genus"]), ha="center", va="center", fontsize=14,
                weight="bold", color="white", zorder=5)
    if split:
        for s, dx in (("A", -0.45), ("B", 0.45)):
            ax.text((gen[0] + lv[s][0]) / 2 + dx, (gen[1] + lv[s][1]) / 2 + 0.05, str(o["shared"][s]),
                    ha="center", va="center", fontsize=13, weight="bold", color=SHARED)
    if i == 3:
        ax.text(cx, 5.45, t["em_iter"], ha="center", va="center", fontsize=12, color=SHARED, weight="bold")
    if i == 4:
        ax.text(lv["C"][0], 5.45, "+7", ha="center", va="center", fontsize=13, weight="bold", color="#3d8b3d")
        ax.text(lv["D"][0], 5.45, "+2", ha="center", va="center", fontsize=13, weight="bold", color=LOWCONF)
    base, scale, bw = 1.35, 3.0 / 90, 0.62
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
            ax.text(x, y + 0.1, str(total), ha="center", va="bottom", fontsize=13, weight="bold", color=INK)
    ax.plot([x0 + 0.25, x0 + w - 0.25], [base, base], color=GRID, lw=1.2, zorder=2)


def main(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for lang, name in (("en", "method.png"), ("es", "metodo.png")):
        t = T[lang]
        fig = plt.figure(figsize=(12.5, 11.2))
        top = fig.add_axes([0, 0.56, 1, 0.44]); top.set_xlim(0, 24); top.set_ylim(0.6, 10); top.axis("off")
        top_panel(top, t, lang)
        mid = fig.add_axes([0, 0.49, 1, 0.07]); mid.set_xlim(0, 24); mid.set_ylim(0, 1); mid.axis("off")
        mid.annotate("", xy=(12, 0.05), xytext=(12, 0.95),
                     arrowprops=dict(arrowstyle="-|>", color=INK, lw=2.6, mutation_scale=26))
        mid.text(12.5, 0.5, t["output"], ha="left", va="center", fontsize=17, weight="bold", color=INK)
        bot = fig.add_axes([0, 0, 1, 0.49]); bot.set_xlim(0, 24); bot.set_ylim(0, 11); bot.axis("off")
        for i in range(5):
            card(bot, i * 4.85, i, t, lang)
        bot.text(12, 0.35, t["legend"], ha="center", va="center", fontsize=12, color=QUIET)
        fig.savefig(outdir / name, dpi=220, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print("written", outdir / name)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
