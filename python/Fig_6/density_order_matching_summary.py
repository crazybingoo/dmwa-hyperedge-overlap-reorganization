"""Create aggregate Fig. 6a summaries from restricted matched-window results.

The input remains in an authorized private workspace. This script writes only
group-level estimates and never exports identifiers or row-level observations.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(os.environ["PAPER1_PROJECT_ROOT"])
INPUT = ROOT / "null_models" / "hyperedge_order_matched_preictal_case_diffs.csv"
OUTPUT = ROOT / "stats" / "results_logic_density_order_matched_patient.csv"
METRICS = ["S", "source_sink_gini_gap", "sink_minus_source_entropy_gap",
           "sink_minus_source_pr_gap", "reciprocal_fraction", "dense_fraction"]
SEED = 20260713


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def bootstrap_mean_ci(values: np.ndarray, n_boot: int = 20_000) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    values = np.asarray(values, dtype=float)
    draw = rng.integers(0, len(values), size=(n_boot, len(values)))
    return tuple(np.quantile(values[draw].mean(axis=1), [0.025, 0.975]))


def main() -> None:
    data = pd.read_csv(INPUT)
    required = {"patient_id", "target_stage", "metric", "target_minus_matched_pre"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    patient = (
        data[data["metric"].isin(METRICS)]
        .groupby(["patient_id", "target_stage", "metric"], as_index=False)
        ["target_minus_matched_pre"].mean()
    )
    rows = []
    for (stage, metric), block in patient.groupby(["target_stage", "metric"]):
        values = block["target_minus_matched_pre"].to_numpy(float)
        lo, hi = bootstrap_mean_ci(values)
        rows.append({"target_stage": stage, "metric": metric, "n_patients": len(values),
                     "mean_target_minus_matched_pre": float(values.mean()),
                     "bootstrap_ci_low": float(lo), "bootstrap_ci_high": float(hi),
                     "wilcoxon_p": float(wilcoxon(values, alternative="two-sided", zero_method="wilcox").pvalue)})
    out = pd.DataFrame(rows)
    out["bh_fdr_p_24_contrasts"] = bh_adjust(out["wilcoxon_p"].to_numpy())
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, index=False)


if __name__ == "__main__":
    main()
