# Privacy and public-data boundary

This repository contains analysis code, simulation outputs and aggregate statistical source data only.

The public release excludes all row-level human-participant information, including raw or processed recordings, participant and seizure labels, individual measurements, clinical annotations, electrode labels, recording dates, internal case codes, linkage keys and local file paths. Supplementary Table 1 and the clinical signal underlying Fig. 1 are not distributed here.

`Source_Data.xlsx` contains only group-level estimates, uncertainty intervals, test summaries, model summaries and simulation results. Counts such as `n_patients` describe the sample size of an aggregate estimate and do not identify any participant.

Scripts that require restricted inputs use the `PAPER1_PROJECT_ROOT` environment variable. Those inputs must remain in an access-controlled institutional workspace and must never be committed.

The public aggregate Fig. 6 renderer intentionally omits participant-level points. It reconstructs the aggregate layers from `Source_Data.xlsx` without accessing restricted inputs.
