# Private Input Schema

The public repository does not include the inputs described below. They are required only for authorized reproduction from restricted clinical data and must remain outside the Git repository.

Set `PAPER1_PROJECT_ROOT` to the authorized private workspace. Expected local folders include:

```text
data_tables/
stats/
figures_data/
simulation/
nature_fig/
notes/
```

## Case metadata

`data_tables/paper1_case_metadata.csv` is private and must not be committed. Required fields are `case_id`, optional `patient_id`, `file_path`, `node_region_path`, and the pre-ictal/ictal/post-ictal window boundaries. `file_path` points to a restricted SEEG `.mat` file containing channels × samples data; `node_region_path` points to private SOZ/PZ/NIZ annotations.

## Threshold inputs

`data_tables/paper2_threshold_elbow_case_summary.csv` supplies one seizure-specific PLV threshold. `data_tables/paper2_threshold_elbow_window_density.csv` supplies window-level density and mean-PLV fields. Both remain private because they retain internal case identifiers before public relabelling.

## Simulation inputs

`simulation/competing_mechanism_raw.csv` contains the abstract model grid used by the Fig. 7 stability analysis. Required feature columns are `S`, `source_sink_gini_gap`, `reciprocal_fraction` and `dense_fraction`.

## Public output boundary

Only the de-identified, derived values copied into `Source_Data.xlsx` are public. Raw recordings, clinical metadata, internal identifiers, linkage keys, absolute paths and intermediate patient/case files must never be added to this repository.
