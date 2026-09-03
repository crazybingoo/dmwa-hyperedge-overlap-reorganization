"""Analyze elbow-threshold Paper 1 comparator metrics.

This script uses the same threshold rule now documented for Paper 2:
one seizure-specific PLV threshold selected from the threshold-density elbow,
then applied to all windows from that seizure.

The analysis asks whether DMWA direction/source-sink variables add a
patient-aware state-space signal beyond pairwise PLV density, hyperedge
density/order descriptors and an undirected hyperedge line-graph comparator.
"""

from __future__ import annotations

import itertools
from pathlib import Path
import os

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(os.environ.get("PAPER1_PROJECT_ROOT", ".")).resolve()
STATS_DIR = PROJECT_ROOT / "stats"
NOTES_DIR = PROJECT_ROOT / "notes"

WINDOW_METRICS = STATS_DIR / "phase9_elbow_threshold_window_metrics.csv"
STAGE_ORDER = ["pre-ictal", "early", "mid", "late", "post-ictal"]
TARGET_STAGES = ["early", "mid", "late", "post-ictal"]

ORDER_EMPTY_ZERO_COLS = [
    "mean_order",
    "frac_order2",
    "frac_order3",
    "order_entropy_norm",
]

FEATURE_SETS = {
    "DMWA_direction_role": [
        "S",
        "source_sink_gini_gap",
        "sink_minus_source_entropy_gap",
        "sink_minus_source_pr_gap",
    ],
    "pairwise_plv_density": [
        "pair_density_elbow",
        "mean_plv",
    ],
    "hyperedge_density_order": [
        "n_hyperedges",
        "mean_order",
        "frac_order2",
        "frac_order3",
        "order_entropy_norm",
    ],
    "undirected_linegraph": [
        "edge_density_undirected",
        "mean_weight_nonzero",
        "weight_sum_norm",
        "degree_gini",
        "degree_entropy_norm",
        "degree_pr_norm",
        "two_path_triple_density",
    ],
    "undirected_plus_order": [
        "edge_density_undirected",
        "mean_weight_nonzero",
        "weight_sum_norm",
        "degree_gini",
        "degree_entropy_norm",
        "degree_pr_norm",
        "two_path_triple_density",
        "n_hyperedges",
        "mean_order",
        "frac_order2",
        "frac_order3",
        "order_entropy_norm",
    ],
}


def exact_signflip_pvalue(diffs: np.ndarray, observed_norm: float) -> tuple[float, int]:
    n = diffs.shape[0]
    if n > 20:
        rng = np.random.default_rng(20260617)
        signs = rng.choice([-1.0, 1.0], size=(200_000, n))
    else:
        signs = np.asarray(list(itertools.product([-1.0, 1.0], repeat=n)), dtype=float)
    null_vecs = signs @ diffs / n
    null_norms = np.linalg.norm(null_vecs, axis=1)
    p_value = (np.sum(null_norms >= observed_norm - 1e-12) + 1.0) / (len(null_norms) + 1.0)
    return float(p_value), int(len(null_norms))


def bootstrap_distance_ci(diffs: np.ndarray, n_boot: int = 20_000) -> tuple[float, float, float]:
    rng = np.random.default_rng(20260617)
    n = diffs.shape[0]
    idx = rng.integers(0, n, size=(n_boot, n))
    boot_mean = diffs[idx].mean(axis=1)
    boot_norm = np.linalg.norm(boot_mean, axis=1)
    return (
        float(np.mean(boot_norm)),
        float(np.quantile(boot_norm, 0.025)),
        float(np.quantile(boot_norm, 0.975)),
    )


def usable_features(df: pd.DataFrame, features: list[str]) -> list[str]:
    usable = []
    for feature in features:
        values = pd.to_numeric(df[feature], errors="coerce")
        if values.notna().all() and float(values.std(ddof=0)) > 0:
            usable.append(feature)
    if not usable:
        raise ValueError(f"No usable nonzero-variance features in {features}")
    return usable


def zscore(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    out = df.copy()
    mu = out[features].mean(axis=0)
    sigma = out[features].std(axis=0, ddof=0)
    if (sigma == 0).any():
        zero_cols = sigma.index[sigma == 0].tolist()
        raise ValueError(f"Cannot z-score zero-variance features: {zero_cols}")
    out[features] = (out[features] - mu) / sigma
    return out


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    window = pd.read_csv(WINDOW_METRICS)
    all_features = sorted({feature for features in FEATURE_SETS.values() for feature in features})

    for col in all_features:
        window[col] = pd.to_numeric(window[col], errors="coerce")
    window[ORDER_EMPTY_ZERO_COLS] = window[ORDER_EMPTY_ZERO_COLS].fillna(0.0)

    keep = ["case_id", "patient_id", "phase5"] + all_features
    case_stage = window[keep].groupby(["case_id", "patient_id", "phase5"], as_index=False)[all_features].mean()
    patient_stage = case_stage.groupby(["patient_id", "phase5"], as_index=False)[all_features].mean()
    patient_stage = patient_stage.rename(columns={"phase5": "stage"})
    patient_stage = patient_stage[patient_stage["stage"].isin(STAGE_ORDER)].copy()
    return window, patient_stage


def paired_distance(df: pd.DataFrame, features: list[str], target_stage: str) -> tuple[dict, pd.DataFrame]:
    features = usable_features(df, features)
    zdf = zscore(df[["patient_id", "stage"] + features], features)
    pre = zdf[zdf["stage"] == "pre-ictal"].set_index("patient_id")
    if target_stage == "ictal-combined":
        target = (
            zdf[zdf["stage"].isin(["early", "mid", "late"])]
            .groupby("patient_id", as_index=True)[features]
            .mean()
        )
    else:
        target = zdf[zdf["stage"] == target_stage].set_index("patient_id")[features]
    patients = sorted(set(pre.index) & set(target.index))
    diffs = target.loc[patients].to_numpy(dtype=float) - pre.loc[patients, features].to_numpy(dtype=float)
    mean_vec = diffs.mean(axis=0)
    distance = float(np.linalg.norm(mean_vec))
    p_value, n_perm = exact_signflip_pvalue(diffs, distance)
    boot_mean, boot_low, boot_high = bootstrap_distance_ci(diffs)

    feature_rows = []
    for feature, mean_diff in zip(features, mean_vec):
        feature_rows.append(
            {
                "target_stage": target_stage,
                "feature": feature,
                "mean_diff_z": float(mean_diff),
                "direction": "increase" if mean_diff > 0 else "decrease" if mean_diff < 0 else "no_change",
            }
        )

    return (
        {
            "target_stage": target_stage,
            "n_patients": len(patients),
            "n_features": len(features),
            "features_used": ";".join(features),
            "centroid_distance_z": distance,
            "bootstrap_mean_distance_z": boot_mean,
            "bootstrap_ci_low_z": boot_low,
            "bootstrap_ci_high_z": boot_high,
            "exact_signflip_p": p_value,
            "n_signflip_permutations": n_perm,
        },
        pd.DataFrame(feature_rows),
    )


def lopo_nearest_centroid(df: pd.DataFrame, features: list[str], target_stage: str) -> dict:
    features = usable_features(df, features)
    if target_stage == "ictal-combined":
        temp = df[df["stage"].isin(["pre-ictal", "early", "mid", "late"])].copy()
        temp["label"] = np.where(temp["stage"] == "pre-ictal", "pre-ictal", "ictal")
        temp = temp.groupby(["patient_id", "label"], as_index=False)[features].mean()
        positive_label = "ictal"
    else:
        temp = df[df["stage"].isin(["pre-ictal", target_stage])].copy()
        temp["label"] = temp["stage"]
        positive_label = target_stage

    patients = sorted(temp["patient_id"].unique())
    y_true: list[str] = []
    y_pred: list[str] = []
    for held_out in patients:
        train = temp[temp["patient_id"] != held_out]
        test = temp[temp["patient_id"] == held_out]
        mu = train[features].mean(axis=0)
        sigma = train[features].std(axis=0, ddof=0).replace(0, np.nan)
        train_z = train.copy()
        test_z = test.copy()
        train_z[features] = (train[features] - mu) / sigma
        test_z[features] = (test[features] - mu) / sigma
        centroids = train_z.groupby("label")[features].mean()
        for _, row in test_z.iterrows():
            distances = {
                label: float(np.linalg.norm(row[features].to_numpy(dtype=float) - centroids.loc[label].to_numpy(dtype=float)))
                for label in centroids.index
            }
            y_true.append(str(row["label"]))
            y_pred.append(min(distances, key=distances.get))

    y_true_arr = np.asarray(y_true)
    y_pred_arr = np.asarray(y_pred)
    labels = ["pre-ictal", positive_label]
    recalls = []
    for label in labels:
        mask = y_true_arr == label
        recalls.append(float(np.mean(y_pred_arr[mask] == label)) if np.any(mask) else np.nan)
    return {
        "n_rows": len(y_true),
        "accuracy": float(np.mean(y_true_arr == y_pred_arr)),
        "balanced_accuracy": float(np.nanmean(recalls)),
        "preictal_recall": recalls[0],
        "target_recall": recalls[1],
    }


def write_report(window: pd.DataFrame, comparison: pd.DataFrame, feature_shifts: pd.DataFrame) -> None:
    out = NOTES_DIR / "phase9_elbow_threshold_comparator_analysis.md"
    combined = comparison[comparison["target_stage"] == "ictal-combined"].sort_values(
        "centroid_distance_z", ascending=False
    )
    density_check = np.nanmax(np.abs(window["pair_density_elbow"] - window["pair_density_from_paper2"]))
    mean_plv_check = np.nanmax(np.abs(window["mean_plv"] - window["mean_plv_from_paper2"]))
    lines = [
        "# Phase 9 elbow-threshold comparator analysis",
        "",
        "Date: 2026-06-17",
        "",
        "Purpose:",
        "",
        "This analysis updates Paper 1 comparator evidence to use Paper 2's threshold-density elbow rule rather than the older per-window 55% PLV quantile threshold.",
        "",
        "Threshold rule:",
        "",
        "- One PLV threshold was selected per seizure from the mean threshold-density curve elbow.",
        "- The same seizure-specific threshold was applied to every 3 s window from that seizure.",
        "- The threshold source is expected in `<PAPER1_PROJECT_ROOT>/data_tables`.",
        "",
        "Extraction QA:",
        "",
        f"- Windows analyzed: {len(window)}.",
        f"- Seizures: {window['case_id'].nunique()}; patients: {window['patient_id'].nunique()}.",
        f"- Threshold range: {window['selected_threshold'].min():.4f} to {window['selected_threshold'].max():.4f}.",
        f"- Max absolute pair-density difference versus Paper 2 window table: {density_check:.6g}.",
        f"- Max absolute mean-PLV difference versus Paper 2 window table: {mean_plv_check:.6g}.",
        "",
        "Combined ictal versus pre-ictal comparator summary:",
        "",
    ]
    for _, row in combined.iterrows():
        lines.append(
            f"- `{row['feature_set']}`: distance = {row['centroid_distance_z']:.3f}, "
            f"95% CI = {row['bootstrap_ci_low_z']:.3f} to {row['bootstrap_ci_high_z']:.3f}, "
            f"p = {row['exact_signflip_p']:.5f}, LOPO balanced accuracy = {row['balanced_accuracy']:.3f}."
        )

    dmwa = combined[combined["feature_set"] == "DMWA_direction_role"].iloc[0]
    undirected = combined[combined["feature_set"] == "undirected_linegraph"].iloc[0]
    pairwise = combined[combined["feature_set"] == "pairwise_plv_density"].iloc[0]

    lines.extend(
        [
            "",
            "Manuscript implication:",
            "",
            f"- The elbow-threshold DMWA direction/role vector shows a combined ictal shift of {dmwa['centroid_distance_z']:.3f} z with exact sign-flip p = {dmwa['exact_signflip_p']:.5f}.",
            f"- The undirected line-graph comparator shows distance = {undirected['centroid_distance_z']:.3f}; pairwise PLV/density shows distance = {pairwise['centroid_distance_z']:.3f}.",
            "- This supports a conservative claim that DMWA source/sink variables capture a directed hyperedge-level organization not identical to pairwise density or an undirected line graph.",
            "- If a comparator performs similarly or better, the manuscript should not claim superiority. It should instead claim nonredundant decomposition and interpret DMWA as the framework that makes source/sink polarity measurable.",
            "",
            "Dominant DMWA direction/role shifts for combined ictal:",
            "",
        ]
    )
    dmwa_shifts = feature_shifts[
        (feature_shifts["feature_set"] == "DMWA_direction_role")
        & (feature_shifts["target_stage"] == "ictal-combined")
    ].copy()
    dmwa_shifts["abs_shift"] = dmwa_shifts["mean_diff_z"].abs()
    dmwa_shifts = dmwa_shifts.sort_values("abs_shift", ascending=False)
    for _, row in dmwa_shifts.iterrows():
        lines.append(f"- `{row['feature']}`: {row['mean_diff_z']:.3f} z ({row['direction']}).")

    lines.extend(
        [
            "",
            "Output files:",
            "",
            "- `<PAPER1_PROJECT_ROOT>/stats/phase9_elbow_threshold_comparator_summary.csv`",
            "- `<PAPER1_PROJECT_ROOT>/stats/phase9_elbow_threshold_feature_shifts.csv`",
            "- `<PAPER1_PROJECT_ROOT>/stats/phase9_elbow_threshold_patient_stage_features.csv`",
        ]
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    window, patient_stage = load_tables()

    rows = []
    feature_rows = []
    for feature_set, features in FEATURE_SETS.items():
        for stage in TARGET_STAGES + ["ictal-combined"]:
            row, shifts = paired_distance(patient_stage, features, stage)
            row.update(lopo_nearest_centroid(patient_stage, features, stage))
            row["feature_set"] = feature_set
            shifts["feature_set"] = feature_set
            rows.append(row)
            feature_rows.append(shifts)

    comparison = pd.DataFrame(rows)
    feature_shifts = pd.concat(feature_rows, ignore_index=True)
    cols = [
        "feature_set",
        "target_stage",
        "n_patients",
        "n_features",
        "features_used",
        "centroid_distance_z",
        "bootstrap_mean_distance_z",
        "bootstrap_ci_low_z",
        "bootstrap_ci_high_z",
        "exact_signflip_p",
        "n_signflip_permutations",
        "n_rows",
        "accuracy",
        "balanced_accuracy",
        "preictal_recall",
        "target_recall",
    ]
    comparison = comparison[cols]

    STATS_DIR.mkdir(parents=True, exist_ok=True)
    NOTES_DIR.mkdir(parents=True, exist_ok=True)
    patient_stage.to_csv(STATS_DIR / "phase9_elbow_threshold_patient_stage_features.csv", index=False)
    comparison.to_csv(STATS_DIR / "phase9_elbow_threshold_comparator_summary.csv", index=False)
    feature_shifts.to_csv(STATS_DIR / "phase9_elbow_threshold_feature_shifts.csv", index=False)
    write_report(window, comparison, feature_shifts)

    print("Phase 9 elbow-threshold comparator analysis complete.")
    print(
        comparison[comparison["target_stage"] == "ictal-combined"][
            [
                "feature_set",
                "centroid_distance_z",
                "bootstrap_ci_low_z",
                "bootstrap_ci_high_z",
                "exact_signflip_p",
                "balanced_accuracy",
            ]
        ].sort_values("centroid_distance_z", ascending=False).to_string(index=False)
    )


if __name__ == "__main__":
    main()

