# Kleiner mock-community benchmark (thesis analysis)

Five ways of turning peptide intensities into strain abundances — Unique, Rigid LCA, Proportional, EM, EM + PEP — compared against the known composition of the Kleiner et al. (2017) mock communities, with paired per-strain statistics.

## Files

| File | What it does |
|---|---|
| `build_taxonomy.py` | Builds the species tree from NCBI (ete3) -> `taxonomy.json` |
| `strain_names.tsv` | Manual NCBI assignments for names that cannot be matched automatically (only F0) |
| `thesis_benchmark.py` | The five strategies, all metrics and statistics |
| `make_result_figures.py` | Result figures in English |
| `results/unfiltered/` | Main results (unfiltered table) |
| `results/filtered/` | Sensitivity analysis on the standard 1 % FDR table |

## Run (PowerShell, from this folder)

```
python build_taxonomy.py --recipe ..\..\files\amounts_species.tsv --out taxonomy.json
python thesis_benchmark.py --mode unfiltered --peptides ..\..\files\unfiltered\fileSDRF_Run_1_2.sdrf_openms_design_triqler_in.tsv --sdrf ..\..\files\fileSDRF_Run_1_2.sdrf.tsv --recipe ..\..\files\amounts_species.tsv --taxonomy taxonomy.json --outdir results\unfiltered
python make_result_figures.py results\unfiltered results\unfiltered
python thesis_benchmark.py --mode filtered --peptides ..\..\files\fileSDRF_Run_1_2.sdrf_openms_design_triqler_in.tsv --sdrf ..\..\files\fileSDRF_Run_1_2.sdrf.tsv --recipe ..\..\files\amounts_species.tsv --taxonomy taxonomy.json --outdir results\filtered
```

Close any open figure before running (Windows cannot overwrite an open image).
The bootstrap seed is fixed, so a re-run reproduces every number exactly.

What to look at first: `results/unfiltered/figure1_method.png`, then `summary.tsv` and `stats.tsv`.
