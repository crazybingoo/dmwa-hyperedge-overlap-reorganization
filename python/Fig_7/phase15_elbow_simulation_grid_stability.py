"""Elbow-threshold empirical target for the simulation grid.

This phase keeps the existing abstract simulation raw grid fixed, but rebuilds
the empirical target vector from elbow-threshold DMWA outputs:

- Phase 9 for S and source-sink Gini gap.
- Phase 13 for reciprocal and dense motif fractions.

The target is constructed from the common set of seizures with finite values
for all four components at pre-ictal, early, mid and late ictal stages. Mean
case-level target-stage minus pre-ictal differences are then averaged over
early, mid and late ictal stages.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get("PAPER1_PROJECT_ROOT", Path(__file__).resolve().parents[2])
).resolve()
CODE_DIR = Path(__file__).resolve().parent
SIM_DIR = PROJECT_ROOT / "simulation"
STATS_DIR = PROJECT_ROOT / "stats"
NOTES_DIR = PROJECT_ROOT / "notes"

sys.path.insert(0, str(CODE_DIR))
import phase10_simulation_grid_stability as phase10  # noqa: E402


PHASE9_WINDOWS = STATS_DIR / "phase9_elbow_threshold_window_metrics.csv"
PHASE13_STAGE = STATS_DIR / "phase13_elbow_motif_stage_summary.csv"
RAW_PATH = SIM_DIR / "competing_mechanism_raw.csv"

EMP_STAGE_DIFF_OUT = STATS_DIR / "phase15_elbow_empirical_stage_differences.csv"
EMP_VECTOR_OUT = STATS_DIR / "phase15_elbow_empirical_vector.csv"
OBSERVED_OUT = STATS_DIR / "phase15_simulation_grid_observed_fit.csv"
BOOT_SUMMARY_OUT = STATS_DIR / "phase15_simulation_grid_bootstrap_summary.csv"
RANKING_OUT = STATS_DIR / "phase15_simulation_grid_ranking_stability.csv"
BOOT_BEST_OUT = STATS_DIR / "phase15_simulation_grid_bootstrap_best_by_model.csv"
FIT_MIRROR_OUT = SIM_DIR / "competing_mechanism_fit_elbow.csv"
NOTE_OUT = NOTES_DIR / "phase15_elbow_simulation_grid_stability.md"

FEATURES = [
    "S",
    "source_sink_gini_gap",
    "reciprocal_fraction",
    "dense_fraction",
]
STAGES = ["early", "mid", "late", "post-ictal"]
ICTAL_STAGES = ["early", "mid", "late"]


def build_elbow_case_stage_table() -> pd.DataFrame:
    phase9 = pd.read_csv(PHASE9_WINDOWS)
    phase9 = phase9.loc[phase9["phase5"].isin(["pre-ictal", *STAGES])].copy()
    role_stage = (
        phase9.groupby(["case_id", "patient_id", "phase5"], as_index=False)[
            ["S", "source_sink_gini_gap"]
        ]
        .mean()
    )

    motif = pd.read_csv(PHASE13_STAGE)
    motif_stage = motif[
        [
            "case_id",
            "patient_id",
            "phase5",
            "reciprocal_fraction_mean",
            "dense_fraction_mean",
        ]
    ].rename(
        columns={
            "reciprocal_fraction_mean": "reciprocal_fraction",
            "dense_fraction_mean": "dense_fraction",
        }
    )

    out = role_stage.merge(
        motif_stage, on=["case_id", "patient_id", "phase5"], how="outer"
    )
    return out.sort_values(["case_id", "phase5"]).reset_index(drop=True)


def summarize_case_stage_diffs(case_stage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in FEATURES:
        for target in STAGES:
            diffs = []
            for case_id, case_df in case_stage.groupby("case_id", sort=False):
                pre = case_df.loc[case_df["phase5"] == "pre-ictal", metric]
                tar = case_df.loc[case_df["phase5"] == target, metric]
                if pre.empty or tar.empty:
                    continue
                a = float(pre.iloc[0])
                b = float(tar.iloc[0])
                if np.isfinite(a) and np.isfinite(b):
                    diffs.append(b - a)
            values = np.asarray(diffs, dtype=float)
            rows.append(
                {
                    "metric": metric,
                    "target_stage": target,
                    "n_cases": int(values.size),
                    "mean_diff": float(np.nanmean(values)) if values.size else np.nan,
                    "median_diff": float(np.nanmedian(values)) if values.size else np.nan,
                    "positive_fraction": float(np.mean(values > 0)) if values.size else np.nan,
                }
            )
    return pd.DataFrame(rows)


def restrict_to_complete_target_cases(
    case_stage: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    required_stages = ["pre-ictal", *ICTAL_STAGES]
    complete_case_ids = []
    for case_id, case_df in case_stage.groupby("case_id", sort=False):
        required = case_df[case_df["phase5"].isin(required_stages)].copy()
        if set(required["phase5"]) != set(required_stages):
            continue
        required = required.set_index("phase5").loc[required_stages, FEATURES]
        if np.isfinite(required.to_numpy(dtype=float)).all():
            complete_case_ids.append(str(case_id))
    restricted = case_stage[
        case_stage["case_id"].astype(str).isin(complete_case_ids)
    ].copy()
    return restricted, complete_case_ids


def empirical_vector_from_stage_diffs(stage_diffs: pd.DataFrame) -> np.ndarray:
    values = []
    for feature in FEATURES:
        idx = (stage_diffs["metric"] == feature) & (
            stage_diffs["target_stage"].isin(ICTAL_STAGES)
        )
        values.append(float(stage_diffs.loc[idx, "mean_diff"].mean()))
    return np.asarray(values, dtype=float)


def make_fit_mirror(observed: pd.DataFrame, emp_vec: np.ndarray) -> pd.DataFrame:
    best_obs = (
        observed.sort_values(
            ["scaled_euclidean_distance", "cosine_similarity"],
            ascending=[True, False],
        )
        .groupby("model", as_index=False)
        .first()
        .sort_values("scaled_euclidean_distance")
    )
    best_obs = best_obs.rename(columns={"beta": "best_beta"})
    best_obs.insert(2, "empirical_S", emp_vec[0])
    best_obs.insert(3, "empirical_source_sink_gini_gap", emp_vec[1])
    best_obs.insert(4, "empirical_reciprocal_fraction", emp_vec[2])
    best_obs.insert(5, "empirical_dense_fraction", emp_vec[3])
    return best_obs


def dataframe_to_markdown(frame: pd.DataFrame) -> str:
    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        cells = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                cells.append(f"{float(value):.6f}")
            else:
                cells.append(str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_report(
    stage_diffs: pd.DataFrame,
    emp_vec: np.ndarray,
    observed: pd.DataFrame,
    rank_summary: pd.DataFrame,
    n_complete_cases: int,
) -> None:
    best_obs = (
        observed.sort_values(
            ["scaled_euclidean_distance", "cosine_similarity"],
            ascending=[True, False],
        )
        .groupby("model", as_index=False)
        .first()
        .sort_values("scaled_euclidean_distance")
    )
    top_model = best_obs.iloc[0]
    top_rank = rank_summary.sort_values("top1_rate", ascending=False).iloc[0]
    hybrid_rank = rank_summary[
        rank_summary["model"] == "role_reciprocity_hybrid"
    ].iloc[0]

    lines = [
        "# Phase 15 elbow-threshold simulation grid stability",
        "",
        "Date: 2026-06-17",
        "",
        "## Purpose",
        "",
        "Rebuilds the empirical simulation target from elbow-threshold DMWA outputs while reusing the existing abstract simulation raw grid.",
        "",
        "## Inputs",
        "",
        f"- Phase 9 elbow-threshold window metrics: `{PHASE9_WINDOWS}`",
        f"- Phase 13 elbow-threshold motif stage metrics: `{PHASE13_STAGE}`",
        f"- Simulation raw grid: `{RAW_PATH}`",
        "",
        "## Empirical Target Definition",
        "",
        f"The target vector uses the common set of {n_complete_cases} seizures with finite values for all four components at pre-ictal, early, mid and late ictal stages. Case-level target-stage minus pre-ictal differences are averaged across early, mid and late ictal stages.",
        "",
        f"- `S`: {emp_vec[0]:.6f}",
        f"- `source_sink_gini_gap`: {emp_vec[1]:.6f}",
        f"- `reciprocal_fraction`: {emp_vec[2]:.6f}",
        f"- `dense_fraction`: {emp_vec[3]:.6f}",
        "",
        "Stage-level inputs:",
        "",
        dataframe_to_markdown(stage_diffs),
        "",
        "## Observed Best-Beta Ranking",
        "",
    ]
    for _, row in best_obs.iterrows():
        lines.append(
            f"- `{row['model']}`: beta = {row['beta']:.1f}, distance = {row['scaled_euclidean_distance']:.3f}, cosine = {row['cosine_similarity']:.3f}, sign agreement = {row['sign_agreement']:.2f}."
        )
    lines.extend(
        [
            "",
            "## Bootstrap Ranking Stability",
            "",
        ]
    )
    for _, row in rank_summary.sort_values("top1_rate", ascending=False).iterrows():
        lines.append(
            f"- `{row['model']}`: top-1 rate = {row['top1_rate']:.3f}, top-2 rate = {row['top2_rate']:.3f}, median rank = {row['median_rank']:.1f}, modal beta = {row['best_beta_mode']:.1f}."
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            f"- Observed best model under the elbow-threshold target: `{top_model['model']}`.",
            f"- Bootstrap top-1 leader under the elbow-threshold target: `{top_rank['model']}`.",
            f"- Hybrid model top-1 rate: {hybrid_rank['top1_rate']:.3f}; top-2 rate: {hybrid_rank['top2_rate']:.3f}.",
            "- If the leading model differs from the legacy Phase10 result, the manuscript should state that the simulation is a stress test of abstract mechanisms rather than a direct fit to the empirical network.",
            "",
            "## Outputs",
            "",
            f"- `{EMP_STAGE_DIFF_OUT}`",
            f"- `{EMP_VECTOR_OUT}`",
            f"- `{OBSERVED_OUT}`",
            f"- `{BOOT_SUMMARY_OUT}`",
            f"- `{RANKING_OUT}`",
            f"- `{BOOT_BEST_OUT}`",
            f"- `{FIT_MIRROR_OUT}`",
        ]
    )
    NOTE_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    NOTES_DIR.mkdir(parents=True, exist_ok=True)
    SIM_DIR.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(RAW_PATH)
    raw["model"] = raw["model"].astype(str)

    case_stage_all = build_elbow_case_stage_table()
    case_stage, complete_case_ids = restrict_to_complete_target_cases(case_stage_all)
    stage_diffs = summarize_case_stage_diffs(case_stage)
    emp_vec = empirical_vector_from_stage_diffs(stage_diffs)

    observed = phase10.compute_observed_grid(raw, emp_vec)
    boot_grid, boot_records = phase10.bootstrap_grid(raw, emp_vec)
    rank_summary, best_by_model = phase10.ranking_stability(boot_records)
    fit_mirror = make_fit_mirror(observed, emp_vec)

    stage_diffs.to_csv(EMP_STAGE_DIFF_OUT, index=False)
    pd.DataFrame(
        [{"feature": feature, "empirical_effect": value} for feature, value in zip(FEATURES, emp_vec)]
    ).to_csv(EMP_VECTOR_OUT, index=False)
    observed.to_csv(OBSERVED_OUT, index=False)
    boot_grid.to_csv(BOOT_SUMMARY_OUT, index=False)
    rank_summary.to_csv(RANKING_OUT, index=False)
    best_by_model.to_csv(BOOT_BEST_OUT, index=False)
    fit_mirror.to_csv(FIT_MIRROR_OUT, index=False)
    write_report(
        stage_diffs,
        emp_vec,
        observed,
        rank_summary,
        n_complete_cases=len(complete_case_ids),
    )

    best_obs = (
        observed.sort_values(
            ["scaled_euclidean_distance", "cosine_similarity"],
            ascending=[True, False],
        )
        .groupby("model", as_index=False)
        .first()
        .sort_values("scaled_euclidean_distance")
    )
    print("Phase 15 elbow simulation grid stability complete.")
    print(f"Complete seizures used for all target components: {len(complete_case_ids)}")
    print("Empirical vector [S, source_sink_gini_gap, reciprocal_fraction, dense_fraction]:")
    print(emp_vec)
    print("\nObserved best-beta ranking:")
    print(
        best_obs[
            [
                "model",
                "beta",
                "scaled_euclidean_distance",
                "cosine_similarity",
                "sign_agreement",
                "dense_fraction_ratio",
            ]
        ].to_string(index=False)
    )
    print("\nBootstrap ranking stability:")
    print(
        rank_summary.sort_values("top1_rate", ascending=False)[
            [
                "model",
                "top1_rate",
                "top2_rate",
                "median_rank",
                "best_beta_mode",
                "selected_distance_median",
                "selected_dense_ratio_median",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
