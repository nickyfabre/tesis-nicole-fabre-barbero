"""Figure 1 (method schematic) and Figure 4 (EM vs EM + PEP by abundance)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

out = Path(sys.argv[2]); res = Path(sys.argv[1])
INK, QUIET, GRID = "#222222", "#666666", "#cccccc"
COL = {"A": "#4e79a7", "B": "#76b7b2", "Genus": "#b07aa1", "C": "#59a14f", "D": "#e15759"}

# ---------------- Figure 1: toy example of the five arms
fig = plt.figure(figsize=(12, 6.2))
ax = fig.add_axes([0.02, 0.05, 0.36, 0.9]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.text(0, 9.6, "Toy community and peptide evidence", fontsize=12, weight="bold", color=INK)
# tree
nodes = {"root": (5, 8.4), "Genus": (3, 7.0), "A": (2, 5.6), "B": (4, 5.6), "C": (6.5, 5.6), "D": (8.5, 5.6)}
for a, b in [("root", "Genus"), ("Genus", "A"), ("Genus", "B"), ("root", "C"), ("root", "D")]:
    (x1, y1), (x2, y2) = nodes[a], nodes[b]
    ax.plot([x1, x2], [y1, y2], color=QUIET, lw=1.2)
for n, (x, y) in nodes.items():
    c = COL.get(n, "white")
    ax.add_patch(plt.Circle((x, y), 0.32, color=c if n != "root" else "white", ec=QUIET, zorder=3))
    lab = {"Genus": "genus\n(A,B)", "D": "D\n(absent)"}.get(n, n)
    ax.text(x, y - 0.75 if n not in ("root", "Genus") else y + 0.5, lab, ha="center", va="top" if n not in ("root", "Genus") else "bottom",
            fontsize=9.5, color=INK)
rows = [("pep1", "A only", "60", "0.001"), ("pep2", "B only", "20", "0.001"), ("pep3", "A and B (shared)", "40", "0.001"),
        ("pep4", "C only", "30", "0.001"), ("pep5", "C only", "10", "0.3"), ("pep6", "D only (false ID)", "20", "0.9")]
ax.text(0.2, 4.0, "peptide   maps to                 intensity   PEP", fontsize=9.5, color=QUIET, family="monospace")
for i, (p, m, I, e) in enumerate(rows):
    y = 3.4 - i * 0.55
    ax.text(0.2, y, f"{p:<9} {m:<23} {I:>6}     {e}", fontsize=9.5, color=INK, family="monospace")
ax.text(0.2, -0.2, "Arms 1-4 keep only PEP < 0.01 (pep5, pep6 dropped).\nArm 5 keeps all, each weighted by 1 - PEP.",
        fontsize=9, color=QUIET)

# bars per arm
arms = {
    "1  Unique":       {"A": 60, "B": 20, "Genus": 0, "C": 30, "D": 0},
    "2  Rigid LCA":    {"A": 60, "B": 20, "Genus": 40, "C": 30, "D": 0},
    "3  Proportional": {"A": 90, "B": 30, "Genus": 0, "C": 30, "D": 0},
    "4  EM":           {"A": 90, "B": 30, "Genus": 0, "C": 30, "D": 0},
    "5  EM + PEP":     {"A": 90, "B": 30, "Genus": 0, "C": 37, "D": 2},
}
notes = {"1  Unique": "shared pep3 discarded", "2  Rigid LCA": "pep3 parked on the genus node",
         "3  Proportional": "pep3 split 3:1, like the unique evidence", "4  EM": "same answer here (little sharing)",
         "5  EM + PEP": "C gains 7 from pep5; D gets 2 from pep6"}
bx = fig.add_axes([0.45, 0.12, 0.53, 0.78])
cats = ["A", "B", "Genus", "C", "D"]
for i, (arm, v) in enumerate(arms.items()):
    y = len(arms) - 1 - i; left = 0
    for c in cats:
        if v[c]:
            bx.barh(y, v[c], left=left, color=COL[c], height=0.55, edgecolor="white")
            bx.text(left + v[c] / 2, y, f"{c}\n{v[c]}" if v[c] >= 14 else str(v[c]), ha="center", va="center", fontsize=8.5, color="white" if v[c] >= 14 else INK)
            left += v[c]
    bx.text(-3, y, arm, ha="right", va="center", fontsize=10.5, weight="bold", color=INK)
    bx.text(left + 3, y, notes[arm], ha="left", va="center", fontsize=9, color=QUIET)
bx.set_xlim(0, 260); bx.set_ylim(-0.6, len(arms) - 0.4); bx.axis("off")
bx.set_title("Where each arm puts the intensity", fontsize=12, weight="bold", color=INK, loc="left", x=-0.25)
fig.savefig(out / "figure1_method.png", dpi=200, bbox_inches="tight")

# ---------------- Figure 4: EM vs EM + PEP per strain, by abundance
ps = pd.read_csv(res / "per_species.tsv", sep="\t")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
for ax, c in zip(axes, ("Constant", "Protein", "Uneven")):
    q = ps[(ps.condition == c) & ps.present]
    a = q[q.arm == "EM"].set_index("strain"); b = q[q.arm == "EM + PEP"].set_index("strain").loc[a.index]
    rare = a.expected < 0.01
    lim = (0.9, max(a.xfold.max(), b.xfold.max()) * 1.6)
    ax.plot(lim, lim, color=GRID, lw=1, zorder=1)
    ax.fill_between(lim, [lim[0]] * 2, lim, color="#e15759", alpha=0.05, zorder=0)
    ax.scatter(a.xfold[~rare], b.xfold[~rare], s=28, color="#4e79a7", label="abundant (≥ 1 %)", zorder=3)
    ax.scatter(a.xfold[rare], b.xfold[rare], s=28, color="#f28e2b", label="rare (< 1 %)", zorder=3)
    for s in a.index:
        r = a.xfold[s] / b.xfold[s]
        if r > 1.8 or r < 1 / 1.8:
            ax.annotate(s, (a.xfold[s], b.xfold[s]), xytext=(4, -10), textcoords="offset points", fontsize=8.5, color=INK)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    better = int((b.log10_error < a.log10_error).sum())
    ax.set_title(f"{c}: {better}/{len(a)} strains better with PEP", fontsize=11)
    ax.set_xlabel("x-fold error, EM with strict threshold")
    ax.text(0.97, 0.04, "below the line = PEP better", transform=ax.transAxes, ha="right", fontsize=8.5, color=QUIET)
axes[0].set_ylabel("x-fold error, EM + PEP (no threshold)")
axes[0].legend(frameon=False, loc="upper left", fontsize=9)
fig.tight_layout()
fig.savefig(out / "figure4_pep_by_abundance.png", dpi=200)
# ---------------- Figure 3: threshold sweep (redrawn from threshold_sweep.tsv)
sw_path = res / "threshold_sweep.tsv"
if sw_path.exists():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from thesis_benchmark import sweep_figure
    sw = pd.read_csv(sw_path, sep="\t")
    sw["threshold"] = sw["threshold"].astype(str)
    sweep_figure(sw, out / "figure_threshold_sweep.png")
print("figures written")
