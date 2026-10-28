# Simulation study (thesis analysis)

The five arms of the Kleiner benchmark, with the same rules and metrics, on simulated communities whose composition is known exactly.

| Study | What varies |
|---|---|
| A | false identifications 0–10 % |
| B | peptide sharing, 83 % to 12 % unique peptides |
| C | as A, with identification confidence depending on intensity (as in the real data) |
| D | leaves-only EM vs all-ranks EM (the MiCId modelling choice), with and without PEP |

Run from this folder (about 10 minutes, fixed seeds):

```
python run_simulations.py --outdir results
python run_simulations.py --outdir results --figure-only     # only redraw the figure
python run_simulations.py --outdir results --only-d          # only study D
```

Outputs in `results/`: `summary.tsv` (means per setting and arm), `stats.tsv` (paired Wilcoxon across communities, Holm-corrected), one table per study with every community, and `figure5_simulations.png`.
