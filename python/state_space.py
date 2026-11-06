"""Primary standardization, exact sign flips and bootstrap implementation."""
from __future__ import annotations
import itertools
import numpy as np
import pandas as pd
STAGE_ORDER=["pre-ictal","early","mid","late","post-ictal"]

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
