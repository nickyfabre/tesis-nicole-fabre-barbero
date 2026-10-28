"""Thesis benchmark: five allocation arms on the Kleiner mock community, with statistics.

One script, one run, everything the Results chapter needs.

    python thesis_benchmark.py \
        --peptides  triqler_unfiltered.tsv   (or the standard 1 % FDR triqler_in.tsv)
        --mode      unfiltered               (or filtered)
        --sdrf      fileSDRF_Run_1_2.sdrf.tsv
        --recipe    amounts_species.tsv
        --taxonomy  taxonomy.json            (made by build_taxonomy.py)
        --outdir    results_unfiltered

The five arms (see the Methods):
    1 Unique        only peptides that map to one strain
    2 Rigid LCA     shared peptides kept at their last common ancestor (no splitting)
    3 Proportional  shared peptides split once, by the unique-peptide anchors
    4 EM            the same split iterated to convergence
    5 EM + PEP      EM, every PSM weighted by 1 - PEP, no confidence threshold

Arms 1-4 use only PSMs with PEP < --threshold (default 0.01), unweighted.
Arm 5 uses every target PSM, weighted by 1 - PEP.
In --mode unfiltered, PEP is re-estimated from target-decoy competition
(taxoshare.pep.estimate_peps); in --mode filtered the reported 1 - searchScore is used.

Outputs (tab-separated, in --outdir):
    strain_sharing.tsv  per strain: fraction of its signal that is shared; defines the sharing stratum
    per_species.tsv   one row per arm x condition x strain: expected, measured, log10 error, x-fold
    summary.tsv       one row per arm x condition: median x-fold (IQR), L1, absent leakage, r, rho
    stats.tsv         paired Wilcoxon + Holm + paired bootstrap for the four planned contrasts
    foldchange.tsv    between-condition log2 fold changes (secondary analysis)
    absent.tsv        estimate on each truly absent strain (negative controls)
    nodes.tsv         every tree node, every arm, every condition (supplementary)
    figure_xfold.png  per-strain x-fold deviation, paired across arms
    threshold_sweep.tsv + figure_threshold_sweep.png  (unfiltered mode only) EM at each PEP
                      threshold, with and without PEP weighting
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# Lukas's package lives in the sibling folder "taxoshare by lukas (boss)"
HERE = Path(__file__).resolve().parent
# the taxoshare library sits at the root of this repository (two levels up)
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT))

from taxoshare.allocate import allocate_species, nodes_from_species, rigid_lca_nodes  # noqa: E402
from taxoshare.pep import estimate_peps  # noqa: E402
from taxoshare.realdata import load_recipe, load_taxonomy  # noqa: E402

CONDITIONS = ("Constant", "Protein", "Uneven")
ARMS = ("Unique", "Rigid LCA", "Proportional", "EM", "EM + PEP")
# planned comparisons: (reference, alternative, what it tests, kind)
CONTRASTS = [
    ("Unique", "Proportional", "Does redistributing shared peptides help?", "difference"),
    ("Rigid LCA", "Proportional", "Splitting to strains vs parking at the ancestor", "difference"),
    ("Proportional", "EM", "Does iterating the split change anything?", "equivalence"),
    ("EM", "EM + PEP", "Does weighting by identification confidence help?", "difference"),
]
NI_MARGIN = 0.05  # non-inferiority margin, log10 units (= 12 % on the x-fold scale)
FLOOR = 1e-6  # measured share floor for log errors (0.0001 %)
SHARING_CUTOFF = 0.50  # a strain is a "sharing strain" if >= 50 % of its signal is shared


# ----------------------------------------------------------------------------- loading

def read_design(sdrf: Path) -> pd.DataFrame:
    """Run number (row order in the SDRF, as quantMS numbers them) -> bio / tech replicate."""
    d = pd.read_csv(sdrf, sep="\t")
    out = pd.DataFrame({
        "run": np.arange(1, len(d) + 1),
        "bio_rep": d["characteristics[biological replicate]"].astype(int).to_numpy(),
        "tech_rep": d["comment[technical replicate]"].astype(int).to_numpy(),
        "file": d["comment[data file]"].to_numpy(),
    })
    return out


def load_psms(path: Path, labels: set[str], mode: str) -> pd.DataFrame:
    """All PSM rows with a candidate strain set and a PEP.

    PEP: target-decoy re-estimate over all PSMs (unfiltered) or 1 - searchScore (filtered).
    """
    cols = ["run", "condition", "searchScore", "intensity", "proteins"]
    df = pd.read_csv(path, sep="\t", usecols=cols)
    df["is_decoy"] = df["proteins"].astype(str).str.contains("DECOY", regex=False)

    if mode == "unfiltered":
        df["pep"] = estimate_peps(df["searchScore"].to_numpy(float), df["is_decoy"].to_numpy())
        n_t, n_d = int((~df.is_decoy).sum()), int(df.is_decoy.sum())
        print(f"target-decoy PEP: {n_t:,} targets, {n_d:,} decoys; "
              f"sum of target PEP = {df.loc[~df.is_decoy, 'pep'].sum():,.0f}")
    else:
        df["pep"] = (1.0 - df["searchScore"].astype(float)).clip(0.0, 1.0)

    df = df[~df["is_decoy"]].copy()
    cache: dict[str, frozenset[str]] = {}

    def strains(acc: str) -> frozenset[str]:
        if acc not in cache:
            cache[acc] = frozenset(a.strip().split("_")[0] for a in acc.split(";")
                                   if a.strip().split("_")[0] in labels)
        return cache[acc]

    df["candidates"] = df["proteins"].astype(str).map(strains)
    df = df[df["candidates"].map(len) > 0].copy()
    df["intensity"] = df["intensity"].astype(float)
    df = df[df["intensity"] > 0]
    return df.drop(columns=["proteins", "is_decoy"])


# ----------------------------------------------------------------------------- allocation

def allocate_run(r: pd.DataFrame, labels: list[str], tree, threshold: float) -> dict[str, dict[str, float]]:
    """Node-level shares (every tree node) for each arm, for one run."""
    strict = r[r["pep"] < threshold]
    sets_s, ints_s = list(strict["candidates"]), strict["intensity"].to_numpy()
    out: dict[str, dict[str, float]] = {}

    for arm, method in (("Unique", "unique"), ("Proportional", "proportional"), ("EM", "em")):
        res = allocate_species(sets_s, ints_s, species=labels, method=method)
        out[arm] = nodes_from_species(tree, res.as_dict())

    lca_nodes = [tree.lca(c) for c in strict["candidates"]]
    out["Rigid LCA"] = rigid_lca_nodes(tree, lca_nodes, ints_s)

    res = allocate_species(list(r["candidates"]), r["intensity"].to_numpy(), species=labels,
                           method="em", peps=r["pep"].to_numpy())
    out["EM + PEP"] = nodes_from_species(tree, res.as_dict())
    return out


def run_all(df, design, labels, tree, threshold) -> pd.DataFrame:
    """Long table: arm, condition, bio_rep, tech_rep, node, share."""
    rows = []
    for run, r in df.groupby("run"):
        info = design.set_index("run").loc[int(run)]
        cond = r["condition"].iloc[0]
        per_arm = allocate_run(r, labels, tree, threshold)
        for arm, shares in per_arm.items():
            for node, v in shares.items():
                rows.append((arm, cond, int(info.bio_rep), int(info.tech_rep), int(run), node, v))
        print(f"  run {int(run):>2} ({cond}, bio {info.bio_rep}, tech {info.tech_rep}) done")
    return pd.DataFrame(rows, columns=["arm", "condition", "bio_rep", "tech_rep", "run", "node", "share"])


def collapse_replicates(long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Technical duplicates -> mean per biological replicate -> mean per condition."""
    bio = long.groupby(["arm", "condition", "bio_rep", "node"], as_index=False)["share"].mean()
    cond = bio.groupby(["arm", "condition", "node"], as_index=False)["share"].mean()
    return bio, cond


# ----------------------------------------------------------------------------- metrics

def strain_sharing(df: pd.DataFrame, labels, threshold: float) -> pd.DataFrame:
    """Per strain: share of the intensity of its candidate PSMs that is shared with another strain."""
    d = df[df["pep"] < threshold]
    shared = d["candidates"].map(len) > 1
    rows = []
    for s in labels:
        has = d["candidates"].map(lambda c: s in c)
        tot = d.loc[has, "intensity"].sum()
        sh = d.loc[has & shared, "intensity"].sum()
        rows.append(dict(strain=s, shared_fraction=sh / tot if tot > 0 else np.nan))
    out = pd.DataFrame(rows)
    out["sharing_strain"] = out["shared_fraction"] >= SHARING_CUTOFF
    return out


def per_species_table(cond_est: pd.DataFrame, labels, truth) -> pd.DataFrame:
    rows = []
    for (arm, c), g in cond_est.groupby(["arm", "condition"]):
        m = g.set_index("node")["share"]
        for i, s in enumerate(labels):
            e = truth[c][i]
            meas = float(m.get(s, 0.0))
            err = abs(np.log10(max(meas, FLOOR) / e)) if e > 0 else np.nan
            rows.append(dict(arm=arm, condition=c, strain=s, expected=e, measured=meas,
                             present=e > 0, floored=(e > 0 and meas < FLOOR),
                             log10_error=err, xfold=10 ** err if e > 0 else np.nan))
    return pd.DataFrame(rows)


def summarise(ps: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (arm, c), g in ps.groupby(["arm", "condition"]):
        p = g[g.present]
        ok = p[p.measured > 0]
        rows.append(dict(
            arm=arm, condition=c, n_present=len(p),
            median_xfold=p.xfold.median(), xfold_Q1=p.xfold.quantile(0.25), xfold_Q3=p.xfold.quantile(0.75),
            L1=float((g.measured - g.expected).abs().sum()),
            absent_leakage_pct=float(g.loc[~g.present, "measured"].sum() * 100),
            log_pearson_r=stats.pearsonr(np.log10(ok.measured), np.log10(ok.expected))[0] if len(ok) > 2 else np.nan,
            spearman_rho=stats.spearmanr(p.measured, p.expected)[0],
            n_floored=int(p.floored.sum()),
        ))
    out = pd.DataFrame(rows)
    out["arm"] = pd.Categorical(out["arm"], ARMS, ordered=True)
    return out.sort_values(["condition", "arm"])


# ----------------------------------------------------------------------------- statistics

def holm(pvals: list[float]) -> list[float]:
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * pvals[idx])
        adj[idx] = min(running, 1.0)
    return adj.tolist()


def compare(ps: pd.DataFrame, n_boot: int, seed: int, sharing: pd.DataFrame) -> pd.DataFrame:
    """Paired, per-strain comparison of two arms within each condition.

    Error = |log10(measured/expected)| per present strain. Wilcoxon signed-rank on the
    paired errors; Holm across the three conditions; bootstrap over strains of the
    difference (alternative - reference) in median error and in log Pearson r,
    both arms evaluated on the same resample.
    """
    rng = np.random.default_rng(seed)
    is_sharing = sharing.set_index("strain")["sharing_strain"]
    rows = []
    for ref, alt, question, kind in CONTRASTS:
      for stratum in ("all strains", "sharing strains", "non-sharing strains"):
        block = []
        for c in CONDITIONS:
            a = ps[(ps.arm == ref) & (ps.condition == c) & ps.present].set_index("strain")
            b = ps[(ps.arm == alt) & (ps.condition == c) & ps.present].set_index("strain")
            a = a.loc[b.index]
            if stratum != "all strains":
                keep = is_sharing.loc[b.index].to_numpy() == (stratum == "sharing strains")
                a, b = a[keep], b[keep]
            ea, eb = a.log10_error.to_numpy(), b.log10_error.to_numpy()
            la = np.log10(np.maximum(a.measured.to_numpy(), FLOOR))
            lb = np.log10(np.maximum(b.measured.to_numpy(), FLOOR))
            le = np.log10(a.expected.to_numpy())
            diff = eb - ea
            if len(diff) < 6:  # too few strains for a meaningful signed-rank test
                p = np.nan
            elif np.allclose(diff, 0):
                p = 1.0
            else:
                p = stats.wilcoxon(eb, ea, zero_method="wilcox", alternative="two-sided").pvalue
            n = len(ea)
            idx = rng.integers(0, n, size=(n_boot, n))
            d_med = np.median(eb[idx], axis=1) - np.median(ea[idx], axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                d_r = np.array([np.corrcoef(lb[i], le[i])[0, 1] - np.corrcoef(la[i], le[i])[0, 1] for i in idx])
            lo, hi = np.nanpercentile(d_med, [2.5, 97.5])
            rlo, rhi = np.nanpercentile(d_r, [2.5, 97.5])
            block.append(dict(
                reference=ref, alternative=alt, question=question, kind=kind, stratum=stratum,
                condition=c, n_strains=n,
                median_fold_improvement=float(np.median(10 ** (ea - eb))),
                median_error_ref=np.median(ea), median_error_alt=np.median(eb),
                median_xfold_ref=10 ** np.median(ea), median_xfold_alt=10 ** np.median(eb),
                strains_alt_better=int((diff < 0).sum()), strains_alt_worse=int((diff > 0).sum()),
                wilcoxon_p=p,
                delta_median_error=np.median(eb) - np.median(ea), delta_median_error_CI95=f"[{lo:.4f}, {hi:.4f}]",
                delta_log_r=np.corrcoef(lb, le)[0, 1] - np.corrcoef(la, le)[0, 1],
                delta_log_r_CI95=f"[{rlo:.4f}, {rhi:.4f}]",
                equivalence_bound_error=max(abs(lo), abs(hi)), equivalence_bound_r=max(abs(rlo), abs(rhi)),
                non_inferior=bool(hi < NI_MARGIN),  # alternative not worse than reference by more than the margin
            ))
        ps_ = [r["wilcoxon_p"] for r in block]
        adj = holm(ps_) if not any(np.isnan(ps_)) else [np.nan] * len(ps_)
        for row, padj in zip(block, adj):
            row["wilcoxon_p_holm"] = padj
            rows.append(row)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- secondary analyses

def fold_changes(cond_est, labels, truth) -> pd.DataFrame:
    """log2 fold change per strain for each pair of designs, observed vs expected."""
    rows = []
    idx = {s: i for i, s in enumerate(labels)}
    for arm in ARMS:
        m = cond_est[cond_est.arm == arm].pivot(index="node", columns="condition", values="share")
        for c1, c2 in (("Constant", "Protein"), ("Constant", "Uneven"), ("Protein", "Uneven")):
            exp, obs, absent = [], [], []
            n_both_absent = 0
            for s in labels:
                e1, e2 = truth[c1][idx[s]], truth[c2][idx[s]]
                if e1 == 0 and e2 == 0:  # absent in both designs: no fold change to recover
                    n_both_absent += 1
                    continue
                o1, o2 = max(m.loc[s, c1], FLOOR), max(m.loc[s, c2], FLOOR)
                obs.append(np.log2(o2 / o1))
                if e1 > 0 and e2 > 0:
                    exp.append(np.log2(e2 / e1)); absent.append(False)
                else:  # absent in one design: expected direction only, ranked at the extreme
                    exp.append(np.inf if e1 == 0 else -np.inf); absent.append(True)
            exp, obs, absent = np.array(exp), np.array(obs), np.array(absent)
            both = ~absent
            rank_exp = stats.rankdata(np.where(np.isinf(exp), np.sign(exp) * 1e9, exp))
            rows.append(dict(
                arm=arm, contrast=f"{c2} vs {c1}", n_both_present=int(both.sum()), n_absent_in_one=int(absent.sum()),
                n_absent_in_both_skipped=n_both_absent,
                pearson_r_present=stats.pearsonr(obs[both], exp[both])[0],
                spearman_rho_all=stats.spearmanr(obs, rank_exp)[0],
                median_abs_log2_error_present=float(np.median(np.abs(obs[both] - exp[both]))),
                absent_direction_correct=(f"{int((np.sign(obs[absent]) == np.sign(exp[absent])).sum())}/{int(absent.sum())}"
                                          if absent.any() else "n/a"),
            ))
    return pd.DataFrame(rows)


def absent_controls(bio_est, labels, truth) -> pd.DataFrame:
    rows = []
    idx = {s: i for i, s in enumerate(labels)}
    for c in ("Constant", "Protein"):
        for s in labels:
            if truth[c][idx[s]] > 0:
                continue
            g = bio_est[(bio_est.condition == c) & (bio_est.node == s)]
            for arm in ARMS:
                v = g[g.arm == arm]["share"] * 100
                rows.append(dict(condition=c, strain=s, arm=arm, mean_pct=v.mean(), min_pct=v.min(), max_pct=v.max()))
    return pd.DataFrame(rows)


SWEEP = (0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.01)  # 1.01 = no threshold


def threshold_sweep(df, design, labels, truth) -> pd.DataFrame:
    """EM at each PEP threshold, unweighted vs PEP-weighted.

    A strict cutoff forces a choice of threshold; this shows how much the answer
    depends on that choice with and without weighting.
    """
    rows = []
    reps = design.set_index("run")["bio_rep"]
    for thr in SWEEP:
        for weighted in (False, True):
            est = []
            for run, r in df.groupby("run"):
                k = r[r["pep"] < thr]
                res = allocate_species(list(k["candidates"]), k["intensity"].to_numpy(), species=labels,
                                       method="em", peps=k["pep"].to_numpy() if weighted else None)
                est.append((r["condition"].iloc[0], reps.loc[int(run)], res.pi))
            for c in CONDITIONS:
                by_bio = {}
                for cc, b, pi in est:
                    if cc == c:
                        by_bio.setdefault(b, []).append(pi)
                m = np.mean([np.mean(v, axis=0) for v in by_bio.values()], axis=0)
                e = truth[c]
                present = e > 0
                err = np.abs(np.log10(np.maximum(m[present], FLOOR) / e[present]))
                rows.append(dict(threshold="none" if thr > 1 else thr, weighting="1 - PEP" if weighted else "none",
                                 condition=c, median_xfold=10 ** np.median(err), L1=float(np.abs(m - e).sum()),
                                 absent_leakage_pct=float(m[~present].sum() * 100)))
        print(f"  sweep threshold {thr} done")
    return pd.DataFrame(rows)


def sweep_figure(sw: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = [("median_xfold", "median x-fold deviation"), ("L1", "L1 distance"),
               ("absent_leakage_pct", "absent leakage (%)")]
    fig, axes = plt.subplots(len(metrics), 3, figsize=(12, 8), sharex=True)
    xs = list(range(len(SWEEP)))
    labels_x = [str(t) if t <= 1 else "none" for t in SWEEP]
    for j, c in enumerate(CONDITIONS):
        for i, (col, name) in enumerate(metrics):
            ax = axes[i, j]
            for w, colour, lab in (("none", "#59a14f", "EM, unweighted"), ("1 - PEP", "#e15759", "EM + PEP")):
                v = sw[(sw.condition == c) & (sw.weighting == w)][col].to_numpy()
                ax.plot(xs, v, marker="o", color=colour, label=lab)
            if i == 0:
                ax.set_title(c)
            if j == 0:
                ax.set_ylabel(name)
            if i == len(metrics) - 1:
                ax.set_xticks(xs, labels_x, rotation=45)
                ax.set_xlabel("PEP threshold (PSMs kept if PEP < threshold)")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Does the result depend on the threshold? (flat line = no)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)


def figure(ps: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
        ax.set_xticks(range(len(ARMS)), ARMS, rotation=30, ha="right")
        ax.set_title(f"{c} (n = {len(wide)} strains)")
        ax.axhline(1, color="black", lw=0.5, ls=":")
    axes[0].set_ylabel("x-fold deviation from input (1 = perfect)")
    fig.suptitle("Per-strain error, paired across arms (black bar = median)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)


# ----------------------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peptides", required=True, type=Path)
    ap.add_argument("--mode", choices=["filtered", "unfiltered"], required=True)
    ap.add_argument("--sdrf", required=True, type=Path)
    ap.add_argument("--recipe", required=True, type=Path)
    ap.add_argument("--taxonomy", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--threshold", type=float, default=0.01, help="PEP cutoff for arms 1-4")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261030)
    a = ap.parse_args()
    warnings.filterwarnings("ignore", category=RuntimeWarning)  # NaN correlations in degenerate bootstrap resamples
    a.outdir.mkdir(parents=True, exist_ok=True)

    labels, truth = load_recipe(a.recipe)
    tree = load_taxonomy(a.taxonomy, labels)
    design = read_design(a.sdrf)
    df = load_psms(a.peptides, set(labels), a.mode)

    shared = df["candidates"].map(len) > 1
    strict = df["pep"] < a.threshold
    print(f"target PSMs: {len(df):,} | below PEP {a.threshold}: {strict.sum():,} | "
          f"shared: {shared.mean():.2%} of rows, {df.intensity[shared].sum() / df.intensity.sum():.2%} of intensity")

    sharing = strain_sharing(df, labels, a.threshold)
    sharing.to_csv(a.outdir / "strain_sharing.tsv", sep="\t", index=False)
    print("sharing strains (>= {:.0%} of signal shared): {}".format(
        SHARING_CUTOFF, ", ".join(sharing.loc[sharing.sharing_strain, "strain"])))

    long = run_all(df, design, labels, tree, a.threshold)
    bio, cond = collapse_replicates(long)
    cond.to_csv(a.outdir / "nodes.tsv", sep="\t", index=False)

    ps = per_species_table(cond, labels, truth)
    ps.to_csv(a.outdir / "per_species.tsv", sep="\t", index=False)
    summ = summarise(ps)
    summ.to_csv(a.outdir / "summary.tsv", sep="\t", index=False)
    st = compare(ps, a.n_boot, a.seed, sharing)
    st.to_csv(a.outdir / "stats.tsv", sep="\t", index=False)
    fold_changes(cond, labels, truth).to_csv(a.outdir / "foldchange.tsv", sep="\t", index=False)
    absent_controls(bio, labels, truth).to_csv(a.outdir / "absent.tsv", sep="\t", index=False)
    figure(ps, a.outdir / "figure_xfold.png")
    if a.mode == "unfiltered":
        print("threshold sweep (EM unweighted vs EM + PEP)...")
        sw = threshold_sweep(df, design, labels, truth)
        sw.to_csv(a.outdir / "threshold_sweep.tsv", sep="\t", index=False)
        try:
            sweep_figure(sw, a.outdir / "figure_threshold_sweep.png")
        except OSError as err:
            print(f"\n[!] could not save figure_threshold_sweep.png ({err}).\n"
                  "    All tables were written. The image is probably open in another program:\n"
                  "    close it and run make_result_figures.py to redraw the figures.")

    pd.set_option("display.width", 200)
    print("\n=== headline (leaves only) ===")
    print(summ[["condition", "arm", "median_xfold", "xfold_Q1", "xfold_Q3", "L1", "absent_leakage_pct",
                "log_pearson_r", "spearman_rho"]].round(4).to_string(index=False))
    print("\n=== planned comparisons ===")
    print(st[["reference", "alternative", "stratum", "condition", "n_strains", "median_xfold_ref", "median_xfold_alt",
              "median_fold_improvement", "strains_alt_better", "strains_alt_worse", "wilcoxon_p_holm", "delta_median_error_CI95",
              "delta_log_r_CI95"]].round(4).to_string(index=False))
    print(f"\nall tables written to {a.outdir}")


if __name__ == "__main__":
    main()
