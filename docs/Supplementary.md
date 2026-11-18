# Supplementary material {.unnumbered}

## I. Code and reproducibility {.unnumbered}

The code, the SDRF file and the result tables are available in the GitHub repository <https://github.com/nickyfabre/tesis-nicole-fabre-barbero>. Paths are relative to the root of the repository.

| Folder or file | Content |
|----------------|---------|
| `taxoshare/` | The library: taxonomy tree and LCA (`tree.py`), the five allocation arms, including EM (`allocate.py`), PEP estimation by target--decoy competition (`pep.py`), metrics, data loaders and the simulator (`simulate.py`) |
| `tests/` | Unit tests of the library |
| `analysis/kleiner_benchmark/` | Analysis of the mock communities: tree building, the five arms, statistics and figures |
| `analysis/simulations/` | Simulation studies A--D |
| `files/` | Input data: SDRF file, mock-community recipe and reference database (the large quantMS tables are not stored in the repository) |

The results are reproduced with the following steps:

1.  Build the species tree from NCBI: `analysis/kleiner_benchmark/build_taxonomy.py` (writes `taxonomy.json`; manual name assignments in Table S1).
2.  Main analysis on the unfiltered quantMS table: `analysis/kleiner_benchmark/thesis_benchmark.py --mode unfiltered`.
3.  Sensitivity analysis on the standard 1 % FDR table: `thesis_benchmark.py --mode filtered`.
4.  Figures: Figures 1–3, `thesis/make_figures_r5.py`; the others, `analysis/kleiner_benchmark/make_result_figures.py` (Spanish versions: `thesis/thesis es/make_figures_es.py`).
5.  Simulations: `analysis/simulations/run_simulations.py` (about 10 minutes; fixed seeds).

## II. quantMS parameters {.unnumbered}

quantMS 1.7.0 was run with `--add_decoys true` and `--add_triqler_output true`. For the unfiltered table, all FDR control points were set to 1.0: `--fdr_threshold`, `--run_fdr_cutoff`, `--psm_level_fdr_cutoff` and `--protein_level_fdr_cutoff`.

## III. Supplementary tables {.unnumbered}

The supplementary tables are tab-separated files in the repository (Supplementary material I). Kleiner tables come from the unfiltered table; the same tables for the standard table (FDR ≤ 1 %) are in `analysis/kleiner_benchmark/results/filtered/`.

| Table | Content | File |
|-------|---------|------|
| S1 | Strains whose name could not be matched to NCBI automatically, with the manual assignment and its source | `analysis/kleiner_benchmark/strain_names.tsv` |
| S2 | Fraction of each strain's intensity that is shared with other strains, and classification as sharing strain | `analysis/kleiner_benchmark/results/unfiltered/strain_sharing.tsv` |
| S3 | Headline metrics per arm and design: median x-fold deviation and IQR, L1 distance, absent leakage, log Pearson *r*, Spearman ρ, number of imputed estimates | `analysis/kleiner_benchmark/results/unfiltered/summary.tsv` |
| S4 | Expected and estimated share, log error and x-fold deviation for every strain, arm and design | `analysis/kleiner_benchmark/results/unfiltered/per_species.tsv` |
| S5 | Planned comparisons between arms: strata, medians, strains improved and worsened, Wilcoxon *p*, bootstrap confidence intervals, equivalence and non-inferiority | `analysis/kleiner_benchmark/results/unfiltered/stats.tsv` |
| S6 | Estimated share at every node of the species tree, per arm and design | `analysis/kleiner_benchmark/results/unfiltered/nodes.tsv` |
| S7 | Threshold analysis: EM at seven PEP thresholds, with and without 1 − PEP weighting | `analysis/kleiner_benchmark/results/unfiltered/threshold_sweep.tsv` |
| S8 | Signal on the strains absent from Constant and Protein, per arm (mean and range over biological replicates) | `analysis/kleiner_benchmark/results/unfiltered/absent.tsv` |
| S9 | Fold-change recovery between designs, per arm | `analysis/kleiner_benchmark/results/unfiltered/foldchange.tsv` |
| S10 | Simulation studies A--C: mean metrics per arm and setting | `analysis/simulations/results/summary.tsv` |
| S11 | Simulation studies A--C: paired tests between arms | `analysis/simulations/results/stats.tsv` |
| S12 | Simulation study D: leaves-only versus all-ranks model | `analysis/simulations/results/studyD_summary.tsv` |
