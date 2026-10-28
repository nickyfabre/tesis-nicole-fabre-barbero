"""Simulation study: the five arms on communities whose composition is known exactly.

The Kleiner benchmark is one dataset with little peptide sharing and few false
identifications. Simulation lets us vary exactly those two things and see when
each ingredient of the method matters.

Same arms, same rules and same metrics as analysis/kleiner_benchmark:

    Unique        unique PSMs only                       (PSMs with PEP < 0.01)
    Rigid LCA     shared PSMs kept at their LCA          (PEP < 0.01)
    Proportional  one EM step from the unique estimate   (PEP < 0.01)
    EM            EM to convergence                      (PEP < 0.01)
    EM + PEP      EM on all PSMs, weighted by 1 - PEP    (no threshold)
    [EM, no threshold, unweighted]  reference only: what happens if all PSMs are used naively

PEPs are estimated from target-decoy competition (taxoshare.pep), as on the real data.

Three studies:
  A  false-identification rate 0, 2, 5, 10 %          (sharing as in the mock community)
  B  peptide sharing from 83 % to 12 % unique peptides (5 % false identifications)
  D  leaves-only EM (this work) vs all-ranks EM (the MiCId modelling choice), with and
     without PEP weighting, at low and high peptide sharing
  C  as A, but the search score of a correct PSM depends on its intensity, so
     low-abundance species have more low-confidence identifications (as observed in
     the real data). This is the situation where a strict threshold discards real evidence.

The unit of replication is the simulated community: each one is drawn independently,
so arms are compared with a paired Wilcoxon test across communities.

    python run_simulations.py --outdir results           (about 2 minutes)
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))  # repository root, where the taxoshare library lives

from taxoshare.allocate import (  # noqa: E402
    all_rank_candidate_sets, allocate_species, node_totals_at_or_below, nodes_from_species, rigid_lca_nodes)
from taxoshare.pep import estimate_peps  # noqa: E402
from taxoshare.simulate import simulate  # noqa: E402

ARMS = ["Unique", "Rigid LCA", "Proportional", "EM", "EM + PEP", "EM, no threshold"]
MAIN_ARMS = ARMS[:5]
CONTRASTS = [("Unique", "Proportional"), ("Proportional", "EM"), ("EM", "EM + PEP")]
THRESHOLD = 0.01
FLOOR = 1e-6
N_COMMUNITIES = 20
N_PEPTIDES = 8000


def intensity_dependent_scores(sim, rng) -> np.ndarray:
    """Study C: a correct PSM's score rises with its log intensity (mean 2.5 to 7.5)."""
    scores = sim.scores.copy()
    n_t = sim.intensities.size
    correct = ~sim.is_false
    z = np.log(sim.intensities[correct])
    z = (z - z.mean()) / z.std()
    scores[:n_t][correct] = rng.normal(np.clip(5.0 + 1.25 * z, 2.5, 7.5), 1.0)
    return scores


def evaluate(pi: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    present = truth > 0
    err = np.abs(np.log10(np.maximum(pi[present], FLOOR) / truth[present]))
    rare = present & (truth < 0.01)
    err_rare = np.abs(np.log10(np.maximum(pi[rare], FLOOR) / truth[rare])) if rare.any() else np.array([np.nan])
    return dict(median_log_error=float(np.median(err)), median_xfold=float(10 ** np.median(err)),
                rare_median_xfold=float(10 ** np.nanmedian(err_rare)),
                L1=float(np.abs(pi - truth).sum()), absent_leakage_pct=float(pi[~present].sum() * 100))


def run_community(seed: int, *, false_rate: float, share: float, score_model: str) -> list[dict]:
    rng = np.random.default_rng(seed)
    sim = simulate(rng, n_peptides=N_PEPTIDES, false_rate=false_rate,
                   within_genus_share=share, across_genus_share=share / 25)
    scores = intensity_dependent_scores(sim, rng) if score_model == "intensity" else sim.scores
    n_t = sim.intensities.size
    peps = estimate_peps(scores, sim.is_decoy)[:n_t] if sim.is_decoy.any() else np.zeros(n_t)
    strict = peps < THRESHOLD
    sets = np.empty(n_t, dtype=object)
    sets[:] = sim.candidate_sets
    ints = sim.intensities
    sp = sim.species

    est = {
        "Unique": allocate_species(list(sets[strict]), ints[strict], species=sp, method="unique").pi,
        "Proportional": allocate_species(list(sets[strict]), ints[strict], species=sp, method="proportional").pi,
        "EM": allocate_species(list(sets[strict]), ints[strict], species=sp, method="em").pi,
        "EM + PEP": allocate_species(list(sets), ints, species=sp, method="em", peps=peps).pi,
        "EM, no threshold": allocate_species(list(sets), ints, species=sp, method="em").pi,
    }
    lca = rigid_lca_nodes(sim.tree, [n for n, k in zip(sim.lca_nodes, strict) if k], ints[strict])
    est["Rigid LCA"] = np.array([lca.get(s, 0.0) for s in sp])

    info = dict(seed=seed, false_rate=false_rate, within_genus_share=share, score_model=score_model,
                unique_pct=100 * float(np.mean([len(c) == 1 for c in sim.candidate_sets])),
                false_kept_strict_pct=100 * float(sim.is_false[strict].mean()) if strict.any() else 0.0,
                correct_lost_strict_pct=100 * float((~strict & ~sim.is_false).sum() / max(1, (~sim.is_false).sum())))
    return [dict(info, arm=a, **evaluate(est[a], sim.true_abundance)) for a in ARMS]


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj, run = np.empty_like(p), 0.0
    for r, i in enumerate(order):
        run = max(run, (len(p) - r) * p[i])
        adj[i] = min(run, 1.0)
    return adj


def paired_tests(df: pd.DataFrame, level_col: str) -> pd.DataFrame:
    rows = []
    for ref, alt in CONTRASTS:
        for metric in ("L1", "median_log_error"):
            block = []
            for lvl, g in df.groupby(level_col):
                w = g.pivot(index="seed", columns="arm", values=metric)
                a, b = w[ref].to_numpy(), w[alt].to_numpy()
                d = b - a
                p = 1.0 if np.allclose(d, 0) else stats.wilcoxon(b, a).pvalue
                block.append(dict(reference=ref, alternative=alt, metric=metric, level=level_col, value=lvl,
                                  mean_ref=a.mean(), mean_alt=b.mean(), mean_difference=d.mean(),
                                  communities_alt_better=int((d < 0).sum()), n=len(d), wilcoxon_p=p))
            for r, padj in zip(block, holm([r["wilcoxon_p"] for r in block])):
                r["wilcoxon_p_holm"] = padj
                rows.append(r)
    return pd.DataFrame(rows)


def figure(dfa: pd.DataFrame, dfb: pd.DataFrame, dfc: pd.DataFrame, path: Path) -> None:
    """Four panels. Rigid LCA is left out: its leaves receive only unique peptides,
    so its leaf error is far off the scale and hides the differences between the others."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colours = {"Unique": "#8c8c8c", "Proportional": "#4e79a7", "EM": "#59a14f",
               "EM + PEP": "#e15759", "EM, no threshold": "#f28e2b"}
    styles = {"EM, no threshold": ":"}
    shown = [a for a in ARMS if a != "Rigid LCA"]
    panels = [
        (dfa, "false_rate", "False identifications (%)", "L1", "A  More identification errors", 100),
        (dfa, "false_rate", "False identifications (%)", "absent_leakage_pct", "A  Signal on absent species", 100),
        (dfb, "unique_pct", "Unique peptides (%)", "L1", "B  More peptide sharing", 1),
        (dfc, "false_rate", "False identifications (%)", "rare_median_xfold",
         "C  Confidence depends on intensity:\n    species below 1 %", 100),
    ]
    ylabels = {"L1": "L1 distance (lower = better)", "absent_leakage_pct": "absent leakage (%)",
               "rare_median_xfold": "median x-fold error (1 = perfect)"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, (df, x, xlab, metric, title, scale) in zip(axes.flat, panels):
        if x == "unique_pct":
            agg = df.groupby(["within_genus_share", "arm"]).agg(x=(x, "mean"), m=(metric, "mean"), s=(metric, "std")).reset_index()
        else:
            agg = df.groupby([x, "arm"]).agg(m=(metric, "mean"), s=(metric, "std")).reset_index()
            agg["x"] = agg[x] * scale
        for arm in shown:
            g = agg[agg.arm == arm].sort_values("x")
            ax.errorbar(g.x, g.m, yerr=g.s / np.sqrt(N_COMMUNITIES), marker="o", ms=4, capsize=2,
                        color=colours[arm], ls=styles.get(arm, "-"), label=arm)
        if x == "unique_pct":
            ax.invert_xaxis()
        ax.set_title(title, loc="left", fontsize=11)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylabels[metric])
    axes[0, 0].legend(frameon=False, fontsize=8.5)
    fig.suptitle(f"Simulated communities (30 species, {N_COMMUNITIES} per point; bars = standard error)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)


def redraw(outdir: Path) -> None:
    """Redraw the figure from the saved tables, without re-running the simulations."""
    r = lambda n: pd.read_csv(outdir / f"study{n}_per_community.tsv", sep="\t")
    figure(r("A_false_rate"), r("B_sharing"), r("C_intensity_scores"), outdir / "figure5_simulations.png")


# ----------------------------------------------------------------------------- study D

D_MODELS = ["EM, leaves", "EM, all ranks", "EM + PEP, leaves", "EM + PEP, all ranks"]


def run_community_d(seed: int, *, share: float) -> list[dict]:
    """Study D: leaves-only EM (ours) vs all-ranks EM (the MiCId modelling choice).

    All-ranks: a peptide's candidates are every node on the path from its LCA down to
    the species that contain it, so a genus-wide peptide supports the genus itself.
    A node's abundance is its own mass plus everything below it. Same evidence, same
    threshold rules as the other studies.
    """
    rng = np.random.default_rng(seed)
    sim = simulate(rng, n_peptides=N_PEPTIDES, false_rate=0.05,
                   within_genus_share=share, across_genus_share=share / 25)
    tree, sp = sim.tree, sim.species
    n_t = sim.intensities.size
    peps = estimate_peps(sim.scores, sim.is_decoy)[:n_t]
    strict = peps < THRESHOLD
    sets = np.empty(n_t, dtype=object)
    sets[:] = sim.candidate_sets
    ar = np.empty(n_t, dtype=object)
    ar[:] = all_rank_candidate_sets(tree, sim.candidate_sets)
    nodes = tree.nodes()
    ints = sim.intensities

    est = {
        "EM, leaves": nodes_from_species(tree, allocate_species(
            list(sets[strict]), ints[strict], species=sp, method="em").as_dict()),
        "EM, all ranks": node_totals_at_or_below(tree, allocate_species(
            list(ar[strict]), ints[strict], species=nodes, method="em").as_dict()),
        "EM + PEP, leaves": nodes_from_species(tree, allocate_species(
            list(sets), ints, species=sp, method="em", peps=peps).as_dict()),
        "EM + PEP, all ranks": node_totals_at_or_below(tree, allocate_species(
            list(ar), ints, species=nodes, method="em", peps=peps).as_dict()),
    }
    truth = nodes_from_species(tree, dict(zip(sp, sim.true_abundance)))
    leaves = set(sp)
    internal = [n for n in nodes if n not in leaves and n != tree.root and truth[n] > 0]
    unique_pct = 100 * float(np.mean([len(c) == 1 for c in sim.candidate_sets]))
    rows = []
    for name, e in est.items():
        pi = np.array([e.get(s, 0.0) for s in sp])
        m = evaluate(pi, sim.true_abundance)
        node_err = [abs(np.log10(max(e.get(n, 0.0), FLOOR) / truth[n])) for n in internal]
        rows.append(dict(seed=seed, within_genus_share=share, unique_pct=unique_pct, model=name,
                         leaf_median_xfold=m["median_xfold"], leaf_L1=m["L1"],
                         genus_median_xfold=float(10 ** np.median(node_err)),
                         absent_leakage_pct=m["absent_leakage_pct"]))
    return rows


def study_d(outdir: Path) -> pd.DataFrame:
    rows = []
    for i, sh in enumerate((0.025, 0.20)):
        for s in range(N_COMMUNITIES):
            rows += run_community_d(40_000 + 100 * i + s, share=sh)
        print(f"  study D, sharing {sh} done")
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "studyD_all_ranks_per_community.tsv", sep="\t", index=False)
    summ = df.groupby(["within_genus_share", "model"], sort=False).agg(
        unique_pct=("unique_pct", "mean"), leaf_median_xfold=("leaf_median_xfold", "mean"),
        leaf_L1=("leaf_L1", "mean"), genus_median_xfold=("genus_median_xfold", "mean"),
        absent_leakage_pct=("absent_leakage_pct", "mean")).reset_index()
    tests = []
    for sh, g in df.groupby("within_genus_share"):
        for ref, alt in (("EM, leaves", "EM, all ranks"), ("EM + PEP, leaves", "EM + PEP, all ranks"),
                         ("EM, leaves", "EM + PEP, leaves"), ("EM, all ranks", "EM + PEP, all ranks")):
            for metric in ("leaf_L1", "genus_median_xfold"):
                w = g.pivot(index="seed", columns="model", values=metric)
                a, b = w[ref].to_numpy(), w[alt].to_numpy()
                p = 1.0 if np.allclose(a, b) else stats.wilcoxon(b, a).pvalue
                tests.append(dict(within_genus_share=sh, reference=ref, alternative=alt, metric=metric,
                                  mean_ref=a.mean(), mean_alt=b.mean(),
                                  communities_alt_better=int((b < a).sum()), n=len(a), wilcoxon_p=p))
    t = pd.DataFrame(tests)
    t["wilcoxon_p_holm"] = holm(t["wilcoxon_p"].to_numpy())
    summ.to_csv(outdir / "studyD_summary.tsv", sep="\t", index=False)
    t.to_csv(outdir / "studyD_stats.tsv", sep="\t", index=False)
    return summ


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    ap.add_argument("--figure-only", action="store_true", help="redraw the figure from saved tables")
    ap.add_argument("--only-d", action="store_true", help="run only study D (all-ranks comparison)")
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    if a.figure_only:
        redraw(a.outdir)
        return
    if a.only_d:
        print(study_d(a.outdir).round(3).to_string(index=False))
        return
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    false_rates = (0.0, 0.02, 0.05, 0.10)
    shares = (0.025, 0.05, 0.10, 0.20, 0.35)
    rows_a, rows_b, rows_c = [], [], []
    for i, fr in enumerate(false_rates):
        for s in range(N_COMMUNITIES):
            rows_a += run_community(10_000 + 100 * i + s, false_rate=fr, share=0.025, score_model="fixed")
            rows_c += run_community(30_000 + 100 * i + s, false_rate=fr, share=0.025, score_model="intensity")
        print(f"  false-ID rate {fr:.0%} done")
    for i, sh in enumerate(shares):
        for s in range(N_COMMUNITIES):
            rows_b += run_community(20_000 + 100 * i + s, false_rate=0.05, share=sh, score_model="fixed")
        print(f"  sharing {sh} done")

    dfa, dfb, dfc = pd.DataFrame(rows_a), pd.DataFrame(rows_b), pd.DataFrame(rows_c)
    for name, df in (("A_false_rate", dfa), ("B_sharing", dfb), ("C_intensity_scores", dfc)):
        df.to_csv(a.outdir / f"study{name}_per_community.tsv", sep="\t", index=False)

    keys = {"A_false_rate": (dfa, "false_rate"), "B_sharing": (dfb, "within_genus_share"),
            "C_intensity_scores": (dfc, "false_rate")}
    summaries, tests = [], []
    for name, (df, lvl) in keys.items():
        s = df.groupby([lvl, "arm"]).agg(
            unique_pct=("unique_pct", "mean"), L1_mean=("L1", "mean"), L1_sd=("L1", "std"),
            median_xfold=("median_xfold", "mean"), rare_median_xfold=("rare_median_xfold", "mean"),
            absent_leakage_pct=("absent_leakage_pct", "mean"),
            false_kept_strict_pct=("false_kept_strict_pct", "mean"),
            correct_lost_strict_pct=("correct_lost_strict_pct", "mean")).reset_index()
        s.insert(0, "study", name)
        summaries.append(s.rename(columns={lvl: "level"}).assign(level_name=lvl))
        t = paired_tests(df, lvl)
        t.insert(0, "study", name)
        tests.append(t)
    summary = pd.concat(summaries)
    summary["arm"] = pd.Categorical(summary["arm"], ARMS, ordered=True)
    summary.sort_values(["study", "level", "arm"]).to_csv(a.outdir / "summary.tsv", sep="\t", index=False)
    pd.concat(tests).to_csv(a.outdir / "stats.tsv", sep="\t", index=False)
    figure(dfa, dfb, dfc, a.outdir / "figure5_simulations.png")

    study_d(a.outdir)
    pd.set_option("display.width", 220)
    print(summary.pivot_table(index=["study", "level"], columns="arm", values="L1_mean", observed=True).round(3))
    print(f"\nall tables written to {a.outdir}")


if __name__ == "__main__":
    main()
