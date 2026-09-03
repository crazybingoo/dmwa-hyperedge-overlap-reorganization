"""Rebuild the Paper 1 joint state-space from threshold-consistent inputs.

Phase9 provides elbow-threshold DMWA direction/role and density/order metrics.
Phase13 provides elbow-threshold reciprocal and dense motif fractions. This
module merges both sources at the patient-stage level before any z-scoring.
"""

from __future__ import annotations

import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get("PAPER1_PROJECT_ROOT", Path(__file__).resolve().parents[2])
).resolve()
STATS_DIR = PROJECT_ROOT / "stats"
DMWA_PATH = STATS_DIR / "phase9_elbow_threshold_patient_stage_features.csv"
MOTIF_PATH = STATS_DIR / "phase13_elbow_motif_stage_summary.csv"

STAGE_ORDER = ["pre-ictal", "early", "mid", "late", "post-ictal"]
TARGET_STAGES = ["early", "mid", "late", "post-ictal"]
ICTAL_STAGES = ["early", "mid", "late"]

STATE_FEATURES = [
    "S",
    "source_sink_gini_gap",
    "sink_minus_source_entropy_gap",
    "sink_minus_source_pr_gap",
    "reciprocal_fraction_mean",
    "dense_fraction_mean",
]

ROLE_FEATURES = [
    "source_sink_gini_gap",
    "sink_minus_source_entropy_gap",
    "sink_minus_source_pr_gap",
]

BLOCK_BALANCED_FEATURES = [
    "S",
    "role_composite",
    "reciprocal_fraction_mean",
    "dense_fraction_mean",
]

FEATURE_SETS = {
    "DMWA_full": STATE_FEATURES,
    "S_only": ["S"],
    "role_only": ROLE_FEATURES,
    "motif_only": ["reciprocal_fraction_mean", "dense_fraction_mean"],
    "density_order_only": [
        "n_hyperedges",
        "mean_order",
        "frac_order2",
        "frac_order3",
        "order_entropy_norm",
    ],
}


def normalize_stage(value: object) -> str:
    return str(value).replace("_", "-")


def load_patient_stage_table(dmwa_path: Path, motif_path: Path) -> pd.DataFrame:
    dmwa = pd.read_csv(dmwa_path)
    motif = pd.read_csv(motif_path)

    if "stage" not in dmwa.columns or "phase5" not in motif.columns:
        raise ValueError("Expected Phase9 'stage' and Phase13 'phase5' columns.")

    dmwa = dmwa.copy()
    motif = motif.copy()
    dmwa["stage"] = dmwa["stage"].map(normalize_stage)
    motif["phase5"] = motif["phase5"].map(normalize_stage)

    required_dmwa = {"patient_id", "stage"}.union(
        {feature for features in FEATURE_SETS.values() for feature in features}
        - {"reciprocal_fraction_mean", "dense_fraction_mean"}
    )
    required_motif = {
        "patient_id",
        "phase5",
        "reciprocal_fraction_mean",
        "dense_fraction_mean",
    }
    missing_dmwa = sorted(required_dmwa - set(dmwa.columns))
    missing_motif = sorted(required_motif - set(motif.columns))
    if missing_dmwa or missing_motif:
        raise ValueError(
            f"Missing required columns: Phase9={missing_dmwa}, Phase13={missing_motif}"
        )

    if dmwa.duplicated(["patient_id", "stage"]).any():
        raise ValueError("Phase9 input has duplicate patient-stage rows.")

    motif_patient_stage = (
        motif.groupby(["patient_id", "phase5"], as_index=False)[
            ["reciprocal_fraction_mean", "dense_fraction_mean"]
        ]
        .mean()
        .rename(columns={"phase5": "stage"})
    )

    merged = dmwa.merge(
        motif_patient_stage,
        on=["patient_id", "stage"],
        how="inner",
        validate="one_to_one",
    )
    merged = merged[merged["stage"].isin(STAGE_ORDER)].copy()

    all_features = sorted({feature for features in FEATURE_SETS.values() for feature in features})
    merged[all_features] = merged[all_features].apply(pd.to_numeric, errors="coerce")
    if merged[all_features].isna().any().any():
        missing = merged.loc[
            merged[all_features].isna().any(axis=1), ["patient_id", "stage"] + all_features
        ]
        raise ValueError(f"Incomplete patient-stage features:\n{missing.to_string(index=False)}")

    expected_stages = set(STAGE_ORDER)
    stage_sets = merged.groupby("patient_id")["stage"].agg(lambda values: set(values))
    incomplete = stage_sets[stage_sets != expected_stages]
    if not incomplete.empty:
        raise ValueError(f"Patients without all five stages: {incomplete.to_dict()}")

    stage_rank = {stage: index for index, stage in enumerate(STAGE_ORDER)}
    merged["_stage_rank"] = merged["stage"].map(stage_rank)
    merged = merged.sort_values(["patient_id", "_stage_rank"]).drop(columns="_stage_rank")
    return merged.reset_index(drop=True)


def zscore_table(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    out = df.copy()
    mu = out[features].mean(axis=0)
    sigma = out[features].std(axis=0, ddof=0)
    if (sigma == 0).any():
        zero_variance = sigma.index[sigma == 0].tolist()
        raise ValueError(f"Cannot z-score zero-variance features: {zero_variance}")
    out[features] = (out[features] - mu) / sigma
    return out


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment with monotonicity enforcement."""
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def exact_signflip_pvalue(
    diffs: np.ndarray, observed_norm: float
) -> tuple[float, int, int, str]:
    """Return a strict exact p value when every sign assignment is enumerated.

    The Monte Carlo +1 correction is used only when the permutation space is
    sampled rather than enumerated. Applying that correction to the complete
    2^n sign-flip distribution would make an exact p value unnecessarily
    conservative.
    """
    n_patients = diffs.shape[0]
    if n_patients > 20:
        rng = np.random.default_rng(20260617)
        signs = rng.choice([-1.0, 1.0], size=(200_000, n_patients))
        mode = "monte_carlo"
    else:
        signs = np.asarray(
            list(itertools.product([-1.0, 1.0], repeat=n_patients)), dtype=float
        )
        mode = "exact_enumeration"
    null_vectors = signs @ diffs / n_patients
    null_norms = np.linalg.norm(null_vectors, axis=1)
    extreme_count = int(np.sum(null_norms >= observed_norm - 1e-12))
    if mode == "exact_enumeration":
        p_value = extreme_count / len(null_norms)
    else:
        p_value = (extreme_count + 1.0) / (len(null_norms) + 1.0)
    return float(p_value), int(len(null_norms)), extreme_count, mode


def patient_stage_array(
    raw_df: pd.DataFrame, features: list[str]
) -> tuple[list[str], np.ndarray]:
    patients = sorted(raw_df["patient_id"].unique())
    indexed = raw_df.set_index(["patient_id", "stage"])
    values = np.empty((len(patients), len(STAGE_ORDER), len(features)), dtype=float)
    for patient_index, patient in enumerate(patients):
        values[patient_index] = indexed.loc[
            [(patient, stage) for stage in STAGE_ORDER], features
        ].to_numpy(dtype=float)
    return patients, values


def standardize_feature_array(values: np.ndarray) -> np.ndarray:
    """Z-score features across the patient and stage axes.

    Supports either one dataset with shape (patient, stage, feature) or a batch
    of bootstrap datasets with shape (bootstrap, patient, stage, feature).
    """
    array = np.asarray(values, dtype=float)
    patient_axis = array.ndim - 3
    stage_axis = array.ndim - 2
    axes = (patient_axis, stage_axis)
    mean = array.mean(axis=axes, keepdims=True)
    sd = array.std(axis=axes, ddof=0, keepdims=True)
    if np.any(sd == 0):
        raise ValueError("Cannot z-score a zero-variance bootstrap feature.")
    return (array - mean) / sd


def block_balanced_coordinates(values: np.ndarray) -> np.ndarray:
    """Collapse the three correlated role metrics into one equal-weight block."""
    standardized = standardize_feature_array(values)
    role_mean = standardized[..., 1:4].mean(axis=-1)
    patient_axis = role_mean.ndim - 2
    stage_axis = role_mean.ndim - 1
    axes = (patient_axis, stage_axis)
    role_mean = (
        role_mean - role_mean.mean(axis=axes, keepdims=True)
    ) / role_mean.std(axis=axes, ddof=0, keepdims=True)
    return np.stack(
        [
            standardized[..., 0],
            role_mean,
            standardized[..., 4],
            standardized[..., 5],
        ],
        axis=-1,
    )


def paired_differences_from_coordinates(
    coordinates: np.ndarray, target_stage: str
) -> np.ndarray:
    pre = coordinates[..., 0, :]
    if target_stage == "ictal-combined":
        target = coordinates[..., 1:4, :].mean(axis=-2)
    else:
        target = coordinates[..., STAGE_ORDER.index(target_stage), :]
    return target - pre


def bootstrap_distance_ci_reestimated(
    raw_df: pd.DataFrame,
    features: list[str],
    target_stage: str,
    n_boot: int,
    *,
    block_balanced: bool = False,
    chunk_size: int = 250,
) -> tuple[float, float, float]:
    """Patient-bootstrap CI with standardization re-estimated per resample."""
    _, values = patient_stage_array(raw_df, features)
    rng = np.random.default_rng(20260617)
    n_patients = values.shape[0]
    indices = rng.integers(0, n_patients, size=(n_boot, n_patients))
    bootstrap_norms = np.empty(n_boot, dtype=float)
    for start in range(0, n_boot, chunk_size):
        stop = min(start + chunk_size, n_boot)
        samples = values[indices[start:stop]]
        if block_balanced:
            coordinates = block_balanced_coordinates(samples)
        else:
            coordinates = standardize_feature_array(samples)
        diffs = paired_differences_from_coordinates(coordinates, target_stage)
        bootstrap_norms[start:stop] = np.linalg.norm(diffs.mean(axis=-2), axis=-1)
    return (
        float(np.mean(bootstrap_norms)),
        float(np.quantile(bootstrap_norms, 0.025)),
        float(np.quantile(bootstrap_norms, 0.975)),
    )


def bootstrap_mean_ci(diffs: np.ndarray, n_boot: int) -> tuple[float, float]:
    """Return a signed patient-bootstrap CI for a one-dimensional mean shift."""
    values = np.asarray(diffs, dtype=float).reshape(-1)
    rng = np.random.default_rng(20260617)
    indices = rng.integers(0, len(values), size=(n_boot, len(values)))
    bootstrap_means = values[indices].mean(axis=1)
    return (
        float(np.quantile(bootstrap_means, 0.025)),
        float(np.quantile(bootstrap_means, 0.975)),
    )


def paired_matrices(
    df: pd.DataFrame, features: list[str], target_stage: str
) -> tuple[list[str], np.ndarray]:
    pre = df[df["stage"] == "pre-ictal"].set_index("patient_id")
    if target_stage == "ictal-combined":
        target = (
            df[df["stage"].isin(ICTAL_STAGES)]
            .groupby("patient_id", as_index=True)[features]
            .mean()
        )
    else:
        target = df[df["stage"] == target_stage].set_index("patient_id")[features]
    patients = sorted(set(pre.index) & set(target.index))
    diffs = target.loc[patients].to_numpy(dtype=float) - pre.loc[
        patients, features
    ].to_numpy(dtype=float)
    return patients, diffs


def paired_result(
    zdf: pd.DataFrame,
    raw_df: pd.DataFrame,
    features: list[str],
    target_stage: str,
    n_boot: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    patients, diffs = paired_matrices(zdf, features, target_stage)
    _, raw_diffs = paired_matrices(raw_df, features, target_stage)
    mean_vector = diffs.mean(axis=0)
    mean_raw_vector = raw_diffs.mean(axis=0)
    distance = float(np.linalg.norm(mean_vector))
    (
        p_value,
        n_permutations,
        extreme_count,
        signflip_mode,
    ) = exact_signflip_pvalue(diffs, distance)
    bootstrap_mean, bootstrap_low, bootstrap_high = (
        bootstrap_distance_ci_reestimated(
            raw_df,
            features,
            target_stage,
            n_boot,
        )
    )

    summary = {
        "target_stage": target_stage,
        "n_patients": len(patients),
        "n_features": len(features),
        "features_used": ";".join(features),
        "centroid_distance_z": distance,
        "bootstrap_mean_distance_z": bootstrap_mean,
        "bootstrap_ci_low_z": bootstrap_low,
        "bootstrap_ci_high_z": bootstrap_high,
        "exact_signflip_p": p_value,
        "n_signflip_permutations": n_permutations,
        "signflip_extreme_count": extreme_count,
        "signflip_mode": signflip_mode,
        "bootstrap_scaling": "reestimated_within_resample",
    }
    feature_rows = pd.DataFrame(
        {
            "target_stage": target_stage,
            "feature": features,
            "mean_diff_z": mean_vector,
            "mean_diff_raw": mean_raw_vector,
            "direction": np.where(
                mean_vector > 0,
                "increase",
                np.where(mean_vector < 0, "decrease", "no_change"),
            ),
        }
    )
    return summary, feature_rows


def trajectory_table(zdf: pd.DataFrame) -> pd.DataFrame:
    centroids = zdf.groupby("stage", as_index=True)[STATE_FEATURES].mean().loc[STAGE_ORDER]
    rows = []
    cumulative = 0.0
    for start, end in zip(STAGE_ORDER[:-1], STAGE_ORDER[1:]):
        step = centroids.loc[end].to_numpy() - centroids.loc[start].to_numpy()
        distance = float(np.linalg.norm(step))
        cumulative += distance
        rows.append(
            {
                "transition": f"{start}_to_{end}",
                "from_stage": start,
                "to_stage": end,
                "step_distance_z": distance,
                "cumulative_distance_z": cumulative,
            }
        )
    return pd.DataFrame(rows)


def lopo_nearest_centroid(
    raw_df: pd.DataFrame, features: list[str], target_stage: str
) -> dict[str, float | int]:
    if target_stage == "ictal-combined":
        temp = raw_df[raw_df["stage"].isin(["pre-ictal"] + ICTAL_STAGES)].copy()
        temp["label"] = np.where(temp["stage"] == "pre-ictal", "pre-ictal", "ictal")
        temp = temp.groupby(["patient_id", "label"], as_index=False)[features].mean()
        positive_label = "ictal"
    else:
        temp = raw_df[raw_df["stage"].isin(["pre-ictal", target_stage])].copy()
        temp["label"] = temp["stage"]
        positive_label = target_stage

    patients = sorted(temp["patient_id"].unique())
    y_true: list[str] = []
    y_pred: list[str] = []
    for held_out in patients:
        train = temp[temp["patient_id"] != held_out].copy()
        test = temp[temp["patient_id"] == held_out].copy()
        mu = train[features].mean(axis=0)
        sigma = train[features].std(axis=0, ddof=0)
        if (sigma == 0).any():
            raise ValueError(f"Zero variance in LOPO training fold for patient {held_out}")
        train[features] = (train[features] - mu) / sigma
        test[features] = (test[features] - mu) / sigma
        centroids = train.groupby("label")[features].mean()
        for _, row in test.iterrows():
            distances = {
                label: float(
                    np.linalg.norm(
                        row[features].to_numpy(dtype=float)
                        - centroids.loc[label].to_numpy(dtype=float)
                    )
                )
                for label in centroids.index
            }
            y_true.append(str(row["label"]))
            y_pred.append(min(distances, key=distances.get))

    y_true_array = np.asarray(y_true)
    y_pred_array = np.asarray(y_pred)
    labels = ["pre-ictal", positive_label]
    recalls = []
    for label in labels:
        mask = y_true_array == label
        recalls.append(float(np.mean(y_pred_array[mask] == label)))
    return {
        "n_rows": len(y_true),
        "accuracy": float(np.mean(y_true_array == y_pred_array)),
        "balanced_accuracy": float(np.mean(recalls)),
        "preictal_recall": recalls[0],
        "target_recall": recalls[1],
    }


def feature_set_comparison(raw_df: pd.DataFrame, n_boot: int) -> pd.DataFrame:
    rows = []
    for feature_set, features in FEATURE_SETS.items():
        zdf = zscore_table(raw_df[["patient_id", "stage"] + features], features)
        summary, _ = paired_result(
            zdf, raw_df, features, "ictal-combined", n_boot=n_boot
        )
        summary.update(lopo_nearest_centroid(raw_df, features, "ictal-combined"))
        summary["feature_set"] = feature_set
        rows.append(summary)
    columns = [
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
        "signflip_extreme_count",
        "signflip_mode",
        "bootstrap_scaling",
        "n_rows",
        "accuracy",
        "balanced_accuracy",
        "preictal_recall",
        "target_recall",
    ]
    return pd.DataFrame(rows)[columns]


def consistency_table(raw_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pre = raw_df[raw_df["stage"] == "pre-ictal"].set_index("patient_id")
    for feature in STATE_FEATURES:
        stage_differences = []
        for stage in ICTAL_STAGES:
            target = raw_df[raw_df["stage"] == stage].set_index("patient_id")
            patients = sorted(set(pre.index) & set(target.index))
            stage_differences.append(
                float((target.loc[patients, feature] - pre.loc[patients, feature]).mean())
            )
        _, combined_diffs = paired_matrices(raw_df, [feature], "ictal-combined")
        combined_difference = float(combined_diffs.mean())
        mean_stage_difference = float(np.mean(stage_differences))
        rows.append(
            {
                "feature": feature,
                "early_difference_raw": stage_differences[0],
                "mid_difference_raw": stage_differences[1],
                "late_difference_raw": stage_differences[2],
                "mean_stage_difference_raw": mean_stage_difference,
                "combined_patient_difference_raw": combined_difference,
                "absolute_consistency_error": abs(
                    combined_difference - mean_stage_difference
                ),
                "post_ictal_weight": 0.0,
                "combined_stage_definition": "early;mid;late",
            }
        )
    result = pd.DataFrame(rows)
    if (result["absolute_consistency_error"] > 1e-12).any():
        raise AssertionError("Combined ictal aggregation does not match early-mid-late mean.")
    return result


def component_inference_table(
    zdf: pd.DataFrame, raw_df: pd.DataFrame, n_boot: int
) -> pd.DataFrame:
    """Audit each state coordinate with paired bootstrap and exact sign flips."""
    rows = []
    for target_stage in TARGET_STAGES + ["ictal-combined"]:
        patients, z_diffs = paired_matrices(zdf, STATE_FEATURES, target_stage)
        _, raw_diffs = paired_matrices(raw_df, STATE_FEATURES, target_stage)
        for index, feature in enumerate(STATE_FEATURES):
            z_values = z_diffs[:, index]
            raw_values = raw_diffs[:, index]
            z_low, z_high = bootstrap_mean_ci(z_values, n_boot)
            raw_low, raw_high = bootstrap_mean_ci(raw_values, n_boot)
            z_mean = float(z_values.mean())
            raw_mean = float(raw_values.mean())
            (
                p_value,
                n_permutations,
                extreme_count,
                signflip_mode,
            ) = exact_signflip_pvalue(z_values.reshape(-1, 1), abs(z_mean))
            rows.append(
                {
                    "target_stage": target_stage,
                    "feature": feature,
                    "n_patients": len(patients),
                    "mean_diff_z": z_mean,
                    "bootstrap_ci_low_z": z_low,
                    "bootstrap_ci_high_z": z_high,
                    "mean_diff_raw": raw_mean,
                    "bootstrap_ci_low_raw": raw_low,
                    "bootstrap_ci_high_raw": raw_high,
                    "exact_signflip_p": p_value,
                    "n_signflip_permutations": n_permutations,
                    "signflip_extreme_count": extreme_count,
                    "signflip_mode": signflip_mode,
                }
            )
    return pd.DataFrame(rows)


def add_stage_family_adjustment(table: pd.DataFrame) -> pd.DataFrame:
    """Treat combined ictal as primary and four stage tests as one BH family."""
    result = table.copy()
    result["inference_role"] = np.where(
        result["target_stage"] == "ictal-combined",
        "primary_combined_ictal",
        "secondary_stage_specific",
    )
    result["secondary_bh_adjusted_p"] = np.nan
    stage_mask = result["target_stage"].isin(TARGET_STAGES)
    result.loc[stage_mask, "secondary_bh_adjusted_p"] = bh_adjust(
        result.loc[stage_mask, "exact_signflip_p"].to_numpy(dtype=float)
    )
    result["multiple_testing_family"] = np.where(
        stage_mask,
        "four_stage_specific_centroid_tests",
        "not_applicable_primary_combined_ictal",
    )
    return result


def block_balanced_analysis(
    patient_stage: pd.DataFrame, n_boot: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    patients, values = patient_stage_array(patient_stage, STATE_FEATURES)
    coordinates = block_balanced_coordinates(values)
    summaries = []
    for target_stage in TARGET_STAGES + ["ictal-combined"]:
        diffs = paired_differences_from_coordinates(coordinates, target_stage)
        mean_vector = diffs.mean(axis=0)
        distance = float(np.linalg.norm(mean_vector))
        (
            p_value,
            n_permutations,
            extreme_count,
            signflip_mode,
        ) = exact_signflip_pvalue(diffs, distance)
        bootstrap_mean, bootstrap_low, bootstrap_high = (
            bootstrap_distance_ci_reestimated(
                patient_stage,
                STATE_FEATURES,
                target_stage,
                n_boot,
                block_balanced=True,
            )
        )
        summaries.append(
            {
                "target_stage": target_stage,
                "n_patients": len(patients),
                "n_features": len(BLOCK_BALANCED_FEATURES),
                "features_used": ";".join(BLOCK_BALANCED_FEATURES),
                "centroid_distance_z": distance,
                "bootstrap_mean_distance_z": bootstrap_mean,
                "bootstrap_ci_low_z": bootstrap_low,
                "bootstrap_ci_high_z": bootstrap_high,
                "exact_signflip_p": p_value,
                "n_signflip_permutations": n_permutations,
                "signflip_extreme_count": extreme_count,
                "signflip_mode": signflip_mode,
                "bootstrap_scaling": (
                    "six_features_and_role_composite_reestimated_within_resample"
                ),
            }
        )

    coordinate_rows = []
    for patient_index, patient in enumerate(patients):
        for stage_index, stage in enumerate(STAGE_ORDER):
            row = {"patient_id": patient, "stage": stage}
            row.update(
                dict(
                    zip(
                        BLOCK_BALANCED_FEATURES,
                        coordinates[patient_index, stage_index],
                    )
                )
            )
            coordinate_rows.append(row)
    return add_stage_family_adjustment(pd.DataFrame(summaries)), pd.DataFrame(
        coordinate_rows
    )


def correlation_audit_table(patient_stage: pd.DataFrame) -> pd.DataFrame:
    """Document role-feature redundancy in observations and paired changes."""
    _, combined_diffs = paired_matrices(
        zscore_table(
            patient_stage[["patient_id", "stage"] + STATE_FEATURES],
            STATE_FEATURES,
        ),
        STATE_FEATURES,
        "ictal-combined",
    )
    scopes = {
        "all_patient_stage_observations": patient_stage[STATE_FEATURES].to_numpy(
            dtype=float
        ),
        "ictal_combined_paired_differences": combined_diffs,
    }
    rows = []
    role_indices = [STATE_FEATURES.index(feature) for feature in ROLE_FEATURES]
    for scope, values in scopes.items():
        role_values = values[:, role_indices]
        correlation = np.corrcoef(role_values, rowvar=False)
        condition_number = float(np.linalg.cond(correlation))
        for first, second in itertools.combinations(range(len(ROLE_FEATURES)), 2):
            rows.append(
                {
                    "scope": scope,
                    "feature_1": ROLE_FEATURES[first],
                    "feature_2": ROLE_FEATURES[second],
                    "pearson_r": float(correlation[first, second]),
                    "role_correlation_matrix_condition_number": condition_number,
                    "n_observations": len(role_values),
                }
            )
    return pd.DataFrame(rows)


def analyze_state_space(
    patient_stage: pd.DataFrame, n_boot: int = 20_000
) -> dict[str, pd.DataFrame]:
    state_z = zscore_table(
        patient_stage[["patient_id", "stage"] + STATE_FEATURES], STATE_FEATURES
    )
    summaries = []
    feature_rows = []
    for stage in TARGET_STAGES + ["ictal-combined"]:
        summary, shifts = paired_result(
            state_z, patient_stage, STATE_FEATURES, stage, n_boot=n_boot
        )
        summaries.append(summary)
        feature_rows.append(shifts)

    patient_stage_output = patient_stage.copy()
    for feature in STATE_FEATURES:
        patient_stage_output[f"{feature}_z"] = state_z[feature].to_numpy()

    block_tests, block_coordinates = block_balanced_analysis(
        patient_stage, n_boot=n_boot
    )
    return {
        "patient_stage_features": patient_stage_output,
        "centroid_tests": add_stage_family_adjustment(pd.DataFrame(summaries)),
        "feature_shifts": pd.concat(feature_rows, ignore_index=True),
        "trajectory": trajectory_table(state_z),
        "feature_set_comparison": feature_set_comparison(patient_stage, n_boot=n_boot),
        "component_inference": component_inference_table(
            state_z, patient_stage, n_boot=n_boot
        ),
        "consistency": consistency_table(patient_stage),
        "block_balanced_centroid_tests": block_tests,
        "block_balanced_patient_stage_coordinates": block_coordinates,
        "role_feature_correlation_audit": correlation_audit_table(patient_stage),
    }


OUTPUT_FILENAMES = {
    "patient_stage_features": "phase17_elbow_joint_state_space_patient_stage_features.csv",
    "centroid_tests": "phase17_elbow_joint_state_space_centroid_tests.csv",
    "feature_shifts": "phase17_elbow_joint_state_space_feature_shifts.csv",
    "trajectory": "phase17_elbow_joint_state_space_trajectory.csv",
    "feature_set_comparison": "phase17_elbow_joint_state_space_feature_set_comparison.csv",
    "component_inference": "phase17_elbow_joint_state_space_component_inference.csv",
    "consistency": "phase17_elbow_joint_state_space_consistency.csv",
    "block_balanced_centroid_tests": (
        "phase17_elbow_joint_state_space_block_balanced_centroid_tests.csv"
    ),
    "block_balanced_patient_stage_coordinates": (
        "phase17_elbow_joint_state_space_block_balanced_patient_stage_coordinates.csv"
    ),
    "role_feature_correlation_audit": (
        "phase17_elbow_joint_state_space_role_feature_correlation_audit.csv"
    ),
}


def write_outputs(outputs: dict[str, pd.DataFrame], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in OUTPUT_FILENAMES.items():
        outputs[key].to_csv(output_dir / filename, index=False)


def main() -> None:
    patient_stage = load_patient_stage_table(DMWA_PATH, MOTIF_PATH)
    outputs = analyze_state_space(patient_stage, n_boot=20_000)
    write_outputs(outputs, STATS_DIR)

    combined = outputs["centroid_tests"].query(
        "target_stage == 'ictal-combined'"
    ).iloc[0]
    dense = outputs["feature_shifts"].query(
        "target_stage == 'ictal-combined' and feature == 'dense_fraction_mean'"
    ).iloc[0]
    dense_inference = outputs["component_inference"].query(
        "target_stage == 'ictal-combined' and feature == 'dense_fraction_mean'"
    ).iloc[0]
    print("Phase17 elbow-threshold joint state-space analysis complete.")
    print(
        f"Combined distance={combined['centroid_distance_z']:.6f}, "
        f"95% CI={combined['bootstrap_ci_low_z']:.6f} to "
        f"{combined['bootstrap_ci_high_z']:.6f}, "
        f"p={combined['exact_signflip_p']:.8f}"
    )
    print(
        f"Dense combined shift={dense['mean_diff_z']:.6f} z "
        f"({dense['mean_diff_raw']:.9f} raw), "
        f"95% CI={dense_inference['bootstrap_ci_low_z']:.6f} to "
        f"{dense_inference['bootstrap_ci_high_z']:.6f}, "
        f"p={dense_inference['exact_signflip_p']:.8f}"
    )


if __name__ == "__main__":
    main()
