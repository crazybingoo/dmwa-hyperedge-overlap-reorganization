"""Patient-level analysis for elbow-threshold S ablation.

Phase 14 recomputes the Original, Coverage-only and DegreeBias-only DMWA
variants from seizure-specific threshold-density elbow thresholds. This script
aggregates those case-level outputs to patient-level summaries for Fig. 6.
"""

from __future__ import annotations

import itertools
from pathlib import Path
import os

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("PAPER1_PROJECT_ROOT", ".")).resolve()
STATS_DIR = ROOT / "stats"
NOTES_DIR = ROOT / "notes"

INPUT = STATS_DIR / "phase14_elbow_s_ablation_stage_summary.csv"
PHASE9_WINDOWS = STATS_DIR / "phase9_elbow_threshold_window_metrics.csv"

PATIENT_STAGE_OUT = STATS_DIR / "phase14_elbow_s_ablation_patient_stage.csv"
SUMMARY_OUT = STATS_DIR / "phase14_elbow_s_ablation_patient_summary.csv"
CONTRAST_OUT = STATS_DIR / "phase14_elbow_s_ablation_method_contrasts.csv"
INTEGRITY_OUT = STATS_DIR / "phase14_elbow_s_ablation_integrity_check.csv"
NOTE_OUT = NOTES_DIR / "phase14_elbow_s_ablation_patient_analysis.md"

PHASES = ["early", "mid", "late", "post-ictal"]
METHODS = [
    ("Original", "Delta_orig"),
    ("Coverage only", "Delta_cov"),
    ("DegreeBias only", "Delta_deg"),
]
CONTRASTS = [
    ("Original - Coverage only", "Delta_orig", "Delta_cov"),
    ("Original - DegreeBias only", "Delta_orig", "Delta_deg"),
    ("Coverage only - DegreeBias only", "Delta_cov", "Delta_deg"),
]


def exact_signflip_p(values: np.ndarray) -> float:
    """Two-sided exact sign-flip p value for patient-level paired effects."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = values.size
    if n == 0:
        return float("nan")
    observed = abs(float(np.mean(values)))
    count = 0
    total = 2**n
    for signs in itertools.product((1.0, -1.0), repeat=n):
        signed_mean = abs(float(np.mean(values * np.asarray(signs))))
        if signed_mean >= observed - 1e-12:
            count += 1
    return count / total


def bootstrap_ci(
    values: np.ndarray, n_boot: int = 10000, seed: int = 20260617
) -> tuple[float, float]:
    """Patient-level bootstrap CI for the mean."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, values.size, size=(n_boot, values.size))
    means = values[idx].mean(axis=1)
    return tuple(np.quantile(means, [0.025, 0.975]))


def holm_adjust(p_values: list[float]) -> list[float]:
    """Holm-adjust p values while preserving input order."""
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full_like(p, np.nan, dtype=float)
    finite_idx = np.where(np.isfinite(p))[0]
    if finite_idx.size == 0:
        return adjusted.tolist()
    order = finite_idx[np.argsort(p[finite_idx])]
    m = order.size
    running = 0.0
    for rank, idx in enumerate(order):
        candidate = (m - rank) * p[idx]
        running = max(running, candidate)
        adjusted[idx] = min(running, 1.0)
    return adjusted.tolist()


def summarize_effects(patient_stage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for phase in PHASES:
        phase_df = patient_stage.loc[patient_stage["phase5"] == phase]
        for method, column in METHODS:
            values = phase_df[column].to_numpy(dtype=float)
            ci_low, ci_high = bootstrap_ci(values)
            rows.append(
                {
                    "phase": phase,
                    "method": method,
                    "n_patients": int(np.isfinite(values).sum()),
                    "mean_delta_S": float(np.nanmean(values)),
                    "median_delta_S": float(np.nanmedian(values)),
                    "ci_low": float(ci_low),
                    "ci_high": float(ci_high),
                    "exact_signflip_p": exact_signflip_p(values),
                }
            )
    out = pd.DataFrame(rows)
    out["holm_p_all_tests"] = holm_adjust(out["exact_signflip_p"].tolist())
    return out


def summarize_contrasts(patient_stage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for phase in PHASES:
        phase_df = patient_stage.loc[patient_stage["phase5"] == phase]
        for label, left, right in CONTRASTS:
            values = (phase_df[left] - phase_df[right]).to_numpy(dtype=float)
            ci_low, ci_high = bootstrap_ci(values)
            rows.append(
                {
                    "phase": phase,
                    "contrast": label,
                    "n_patients": int(np.isfinite(values).sum()),
                    "mean_difference": float(np.nanmean(values)),
                    "median_difference": float(np.nanmedian(values)),
                    "ci_low": float(ci_low),
                    "ci_high": float(ci_high),
                    "exact_signflip_p": exact_signflip_p(values),
                }
            )
    out = pd.DataFrame(rows)
    out["holm_p_all_tests"] = holm_adjust(out["exact_signflip_p"].tolist())
    return out


def check_phase9_alignment(stage_df: pd.DataFrame) -> pd.DataFrame:
    """Confirm Phase14 Original S reproduces Phase9 elbow-threshold S."""
    phase9 = pd.read_csv(PHASE9_WINDOWS)
    phase9 = phase9.loc[phase9["phase5"].isin(["pre-ictal", *PHASES])].copy()
    phase9_stage = (
        phase9.groupby(["case_id", "patient_id", "phase5"], as_index=False)["S"]
        .mean()
        .rename(columns={"S": "phase9_S_mean"})
    )
    merged = stage_df.merge(
        phase9_stage, on=["case_id", "patient_id", "phase5"], how="outer"
    )
    merged["S_orig_minus_phase9"] = merged["S_orig_mean"] - merged["phase9_S_mean"]
    return merged[
        ["case_id", "patient_id", "phase5", "S_orig_mean", "phase9_S_mean", "S_orig_minus_phase9"]
    ]


def markdown_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False, floatfmt=".6f")
    except ImportError:
        return "```\n" + df.to_string(index=False) + "\n```"


def write_note(
    stage_df: pd.DataFrame,
    patient_stage: pd.DataFrame,
    summary: pd.DataFrame,
    contrasts: pd.DataFrame,
    integrity: pd.DataFrame,
) -> None:
    max_abs_diff = float(integrity["S_orig_minus_phase9"].abs().max())
    lines = [
        "# Phase 14 Elbow-Threshold S Ablation Patient Analysis",
        "",
        "Date: 2026-06-17",
        "",
        "## Purpose",
        "",
        "Replaces the legacy Phase 8 S ablation inputs for Fig. 6 with patient-level summaries recomputed from seizure-specific threshold-density elbow thresholds.",
        "",
        "## Inputs",
        "",
        f"- Phase 14 stage table: `{INPUT}`",
        f"- Phase 9 elbow-threshold window table for integrity check: `{PHASE9_WINDOWS}`",
        f"- Case-stage rows: {len(stage_df)}",
        f"- Unique cases/seizures: {stage_df['case_id'].nunique()}",
        f"- Unique patients after case aggregation: {patient_stage['patient_id'].nunique()}",
        "",
        "Multiple seizures from the same patient were averaged within each stage before statistical testing.",
        "",
        "## Integrity Check",
        "",
        f"- Maximum absolute difference between Phase14 Original S and Phase9 stage-mean S: {max_abs_diff:.3e}.",
        "",
        "## Patient-Level Effects",
        "",
        markdown_table(summary),
        "",
        "## Method Contrasts",
        "",
        markdown_table(contrasts),
        "",
        "## Interpretation Boundary",
        "",
        "- This analysis tests whether the S increase survives when Coverage or DegreeBias is isolated under the same elbow-threshold hypergraph construction.",
        "- It should be framed as component sensitivity, not as proof that the full DMWA formula outperforms every reduced variant for S alone.",
        "- Fig. 6 still requires the simulation target to be tied to elbow-threshold empirical features before the entire figure can be described as fully elbow-threshold based.",
        "",
        "## Outputs",
        "",
        f"- `{PATIENT_STAGE_OUT}`",
        f"- `{SUMMARY_OUT}`",
        f"- `{CONTRAST_OUT}`",
        f"- `{INTEGRITY_OUT}`",
    ]
    NOTE_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    if not PHASE9_WINDOWS.exists():
        raise FileNotFoundError(PHASE9_WINDOWS)

    STATS_DIR.mkdir(parents=True, exist_ok=True)
    NOTES_DIR.mkdir(parents=True, exist_ok=True)

    stage_df = pd.read_csv(INPUT)
    required = {
        "case_id",
        "patient_id",
        "phase5",
        "S_orig_mean",
        "Delta_orig",
        "Delta_cov",
        "Delta_deg",
    }
    missing = required - set(stage_df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    integrity = check_phase9_alignment(stage_df)
    max_abs_diff = integrity["S_orig_minus_phase9"].abs().max()
    if not np.isfinite(max_abs_diff) or max_abs_diff > 1e-10:
        raise ValueError(
            "Phase14 Original S does not reproduce Phase9 stage means; "
            f"max absolute difference = {max_abs_diff}"
        )

    target = stage_df.loc[stage_df["phase5"].isin(PHASES)].copy()
    patient_stage = (
        target.groupby(["patient_id", "phase5"], as_index=False)[
            ["Delta_orig", "Delta_cov", "Delta_deg"]
        ]
        .mean()
        .sort_values(["patient_id", "phase5"])
    )

    summary = summarize_effects(patient_stage)
    contrasts = summarize_contrasts(patient_stage)

    patient_stage.to_csv(PATIENT_STAGE_OUT, index=False)
    summary.to_csv(SUMMARY_OUT, index=False)
    contrasts.to_csv(CONTRAST_OUT, index=False)
    integrity.to_csv(INTEGRITY_OUT, index=False)
    write_note(stage_df, patient_stage, summary, contrasts, integrity)

    print(f"Wrote {PATIENT_STAGE_OUT}")
    print(f"Wrote {SUMMARY_OUT}")
    print(f"Wrote {CONTRAST_OUT}")
    print(f"Wrote {INTEGRITY_OUT}")
    print(f"Wrote {NOTE_OUT}")
    print(f"Patients: {patient_stage['patient_id'].nunique()}")
    print(f"Maximum absolute Phase14-vs-Phase9 S difference: {max_abs_diff:.3e}")


if __name__ == "__main__":
    main()

