# tesis-nicole-fabre

Code and results for the Licenciatura thesis (Universidad de Buenos Aires, Facultad de Ciencias Exactas y Naturales):

**Asignación filogenética proporcional de péptidos compartidos para la estimación de biomasa microbiana en metaproteómica**
*(Proportional phylogenetic allocation of shared peptides for microbial biomass estimation in metaproteomics)*

Author: Nicole Fabre Barbero. Director: Dr. Lukas Käll (SciLifeLab, KTH). Co-director: Dr. Germán Rosano.

Peptide intensities from quantMS are allocated over a species tree. Five strategies are compared: discarding shared peptides (Unique), keeping them at the lowest common ancestor (Rigid LCA), splitting them among candidate strains once (Proportional) or iteratively (EM), and weighting every identification by 1 − PEP instead of using a confidence threshold (EM + PEP). The strategies are evaluated on the mock communities of Kleiner et al. (2017, PRIDE PXD006118) and on simulated communities.

## Contents

| Folder or file | Content |
|---|---|
| `taxoshare/` | The library: species tree and LCA (`tree.py`), the five allocation strategies including EM (`allocate.py`), PEP estimation by target–decoy competition (`pep.py`), metrics, data loading and the simulator (`simulate.py`) |
| `tests/` | Unit tests of the library |
| `analysis/kleiner_benchmark/` | Mock-community analysis: tree building, the five strategies, statistics, figures and result tables (`results/unfiltered/`, `results/filtered/`) |
| `analysis/simulations/` | Simulation studies A–D and their result tables (`results/`) |
| `files/` | Inputs: SDRF file, mock-community recipe (`amounts_species.tsv`) and reference database (`Mock_Comm_RefDB_V3.fasta`) |
| `thesis/` | Scripts that draw the thesis figures, and the figures (`figures/` in English, `thesis es/figuras/` in Spanish) |

The supplementary tables of the thesis (S1–S12) are the `.tsv` files listed in Supplementary material III.

## Installation

```
pip install -r requirements.txt
pytest tests
```

## Reproducing the results

The quantMS peptide tables (about 2.96 million PSMs, more than 200 MB) are not stored here. They are obtained by processing PRIDE PXD006118 (runs 1 and 2, 24 RAW files) with quantMS 1.7.0, using `files/fileSDRF_Run_1_2.sdrf.tsv` and `files/Mock_Comm_RefDB_V3.fasta`, with `--add_decoys true` and `--add_triqler_output true`. For the unfiltered table, all FDR control points are set to 1.0 (`--fdr_threshold`, `--run_fdr_cutoff`, `--psm_level_fdr_cutoff`, `--protein_level_fdr_cutoff`). Place the resulting `*_triqler_in.tsv` files in `files/unfiltered/` (unfiltered search) and `files/` (standard 1 % FDR search).

From `analysis/kleiner_benchmark/`:

1. Species tree from NCBI: `python build_taxonomy.py --recipe ../../files/amounts_species.tsv --out taxonomy.json`
2. Main analysis (unfiltered table):
   `python thesis_benchmark.py --mode unfiltered --peptides ../../files/unfiltered/fileSDRF_Run_1_2.sdrf_openms_design_triqler_in.tsv --sdrf ../../files/fileSDRF_Run_1_2.sdrf.tsv --recipe ../../files/amounts_species.tsv --taxonomy taxonomy.json --outdir results/unfiltered`
3. Sensitivity analysis (standard 1 % FDR table): the same command with `--mode filtered`, `--peptides ../../files/fileSDRF_Run_1_2.sdrf_openms_design_triqler_in.tsv` and `--outdir results/filtered`.
4. Figures: `python make_result_figures.py results/unfiltered results/unfiltered`; thesis figures: `python thesis/make_figures_r5.py` and `python "thesis/thesis es/make_figures_es.py"` from the repository root.

Simulations (about 10 minutes): `python analysis/simulations/run_simulations.py`.

The bootstrap and the simulations use fixed random seeds, so every number in the thesis is reproduced exactly.

## Acknowledgements

Computations were run on the Arrhenius cluster, provided by the National Academic Infrastructure for Supercomputing in Sweden (NAISS), funded by the Swedish Research Council.

## License

Apache License 2.0 (see `LICENSE`).
