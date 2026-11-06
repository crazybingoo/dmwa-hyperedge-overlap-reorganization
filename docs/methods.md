# Core implementation and conventions

`primary_pipeline.py` computes Hilbert phases, pairwise PLV and triplet circular omega complexity (COC). COC is one minus the normalized eigenvalue entropy of the 3 x 3 circular-correlation matrix. Each score family uses its density-curve elbow over thresholds 0:0.0025:1; the selected thresholds are held across the analyzed windows.

`dmwa.py` takes a binary hyperedge-by-member incidence matrix containing pairs and triplets. For overlap c_ij, order k_i and overlap degree d_i, W_ij = (c_ij/k_i) d_j/(d_i+d_j) for equal orders, and W_ij = c_ij/k_i for different orders. The diagonal is zero. S = sum(min(W_ij,W_ji))/sum(max(W_ij,W_ji)) over unordered pairs. Undefined denominators produce NaN. Roles use q_i = (out_i-in_i)/(out_i+in_i+epsilon).

Primary role gaps use uncorrected Gini values, retain zero loads, and are undefined without positive mass. `D_G_finite_size` is the sensitivity result multiplied by m/(m-1). Entropy and participation-ratio gaps use sink minus source; the Gini gap uses source minus sink.

The exact induced 13-triad census is available as a direct enumeration and a bitset implementation. `kernels.py` computes the three reported motif families at the primary dyadic ratio cutoff 0.85. Bit masks support at most 64 members; the incidence-matrix implementation does not have that encoding limit. The participation-preserving null kernels use signed int64 masks, four chains and a replicate count divisible by four. Their dense membership lookup uses 2^n entries and is intended for the small member sets used here. The synthetic test demonstrates the required inputs.

Main Fig. 5 uses population SDs (ddof=0) over 70 participant-stage observations. Supplementary Fig. 1 uses sample SDs (ddof=1) and keeps those scales fixed during bootstrap resampling. The combined ictal displacement is 2.27313813062039 in the main analysis and 2.25684302353242 in the supplementary construction analysis. These are separately specified analysis conventions.

`held_out_comparison.py` contains the fixed logistic comparison with scaling learned within each training fold. The public data report only aggregate Brier scores. Clinical analyses require authorized private inputs; the public aggregate renderer uses only the included CSV files.

`simulate.py` reproduces independent Kuramoto, common-input and degree-constrained hypergraph experiments. Its output folder is created under `outputs/simulation`. The supplied source tables are the completed experiment summaries; running the full simulation ensemble is computationally intensive.
