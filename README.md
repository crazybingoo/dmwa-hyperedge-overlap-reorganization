# DMWA role–reciprocity: revised code and aggregate source data

Release 0.3.0 accompanies the November 2026 revision: main Figs. 2–8 and Supplementary Figs. 1–10.

- `python/`: PLV/COC selection, DMWA weights, source/sink gaps, exact motifs, hypergraph controls, held-out model comparison and independent simulations.
- `matlab/`: compact DMWA construction reference implementation.
- `data/`: reviewed group summaries and simulation results. `data/catalog.csv` maps every dataset to its analysis unit; filenames identify the associated figure.
- `r_figures/render_aggregates.R`: executable R plotting implementation for the public aggregate results.
- `tests/verify_core.py`: synthetic parity, structural identities, motif census and null-preservation checks.
- `docs/methods.md`: formulas and analysis conventions, including the main/SI standardization distinction.

## Run

From the repository root, with Python 3.10+ and R 4.3+:

```sh
python -m pip install -r requirements.txt
python tests/verify_core.py
Rscript r_figures/render_aggregates.R
```

The renderer uses base R and writes PDF/PNG aggregate panels to `outputs/aggregate_figures`. No private inputs are needed. These panels expose the released numerical summaries; the manuscript's individual-point layers and clinical traces are excluded.

Optional full independent simulations:

```sh
python python/simulate.py --pilot
python python/simulate.py
python python/kuramoto_count_null.py
```

Python core functions can also be imported with `python/` on the module search path. Input matrices and synthetic examples are described in `docs/methods.md` and the tests.

## Data availability

Only group-level human-derived results and simulation summaries are public. Raw/processed SEEG and participant-level measurements are not included. See `docs/privacy_and_data_boundary.md`. CSV values preserve the source analysis precision; manuscript displays use rounded values.

The latest release supersedes the earlier figure numbering and source workbook. Git history retains the prior public release. Cite the associated manuscript and this repository using `CITATION.cff`.
