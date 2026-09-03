"""Aggregate and audit the corrected Phase 11 W-level null experiment.

The MATLAB runner writes versioned, per-case case-stage summaries under
``stats/w_level_null_v2/<tag>/raw``.  This script keeps the inferential unit at
the patient level, performs the Monte Carlo analysis, and writes all derived
files back to the same versioned directory.  It never reads from or writes to
``Manuscript.docx``.

Typical production invocation::

    python phase11_w_level_null_v2_analysis.py --tag production_B5000_YYYYMMDD_HHMMSS

Small preflight runs can override ``--expected-iterations`` and
``--bootstrap-reps``.  Production defaults are B=5000 and 20,000 bootstrap
replicates.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


DEFAULT_PROJECT_ROOT = Path(
    os.environ.get("PAPER1_PROJECT_ROOT", Path(__file__).resolve().parents[2])
).resolve()
DEFAULT_EXPECTED_ITERATIONS = 5_000
DEFAULT_BOOTSTRAP_REPS = 20_000
DEFAULT_BOOTSTRAP_SEED = 20260715
DEFAULT_EXPECTED_MASTER_SEED = 20260715

PHASE9_RELATIVE_PATH = Path("stats") / "phase9_elbow_threshold_window_metrics.csv"
V2_RELATIVE_ROOT = Path("stats") / "w_level_null_v2"

STAGES = ["pre-ictal", "early", "mid", "late", "post-ictal"]
TARGET_STAGES = ["early", "mid", "late", "post-ictal"]
METRICS = [
    "S",
    "source_sink_gini_gap",
    "sink_minus_source_entropy_gap",
    "sink_minus_source_pr_gap",
]
NULL_MODELS = ["weight_shuffle", "orientation_pair_swap"]

EMPIRICAL_REQUIRED = ["case_id", "patient_id", "phase5", "n_windows", *METRICS]
NULL_REQUIRED = ["null_model", "iter", *EMPIRICAL_REQUIRED]

FLOAT_TOLERANCE = 1e-10
ORIENTATION_S_TOLERANCE = 1e-10
CONVERGENCE_MIN_FINAL_B = 1_000
CONVERGENCE_ALPHA = 0.05
CONVERGENCE_MAX_MEAN_DRIFT_Z = 3.0
CONVERGENCE_MAX_P_DRIFT_Z = 3.0
CONVERGENCE_MAX_INTERVAL_ENDPOINT_FRACTION = 0.10

SUMMARY_FILENAME = "phase11_w_level_null_v2_summary.csv"
FIGURE_SOURCE_FILENAME = "phase11_w_level_null_v2_figure_source.csv"
REPORT_FILENAME = "phase11_w_level_null_v2_analysis_report.md"


@dataclass(frozen=True)
class AnalysisPaths:
    project_root: Path
    tag: str

    @property
    def tag_dir(self) -> Path:
        return self.project_root / V2_RELATIVE_ROOT / self.tag

    @property
    def raw_dir(self) -> Path:
        candidate = self.tag_dir / "raw"
        return candidate if candidate.is_dir() else self.tag_dir

    @property
    def manifest_dir(self) -> Path:
        candidate = self.tag_dir / "manifest"
        return candidate if candidate.is_dir() else self.tag_dir

    @property
    def phase9(self) -> Path:
        return self.project_root / PHASE9_RELATIVE_PATH


def _stable_seed(base_seed: int, *parts: object) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return (int(base_seed) + int.from_bytes(digest[:8], "little")) % (2**63 - 1)


def _require_columns(df: pd.DataFrame, required: Sequence[str], label: str) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def _normalise_stage(value: object) -> str:
    text = str(value).strip().lower().replace("_", "-").replace(" ", "-")
    aliases = {
        "preictal": "pre-ictal",
        "pre-ictal": "pre-ictal",
        "early": "early",
        "mid": "mid",
        "middle": "mid",
        "late": "late",
        "post": "post-ictal",
        "postictal": "post-ictal",
        "post-ictal": "post-ictal",
    }
    if text not in aliases:
        raise ValueError(f"Unrecognised phase5/stage value: {value!r}")
    return aliases[text]


def _discover_files(directory: Path, pattern: str) -> list[Path]:
    files = sorted(directory.glob(pattern))
    if not files and directory.name != "raw":
        files = sorted(directory.rglob(pattern))
    return files


def _read_csvs(files: Sequence[Path], label: str) -> pd.DataFrame:
    if not files:
        raise FileNotFoundError(f"No {label} files were found")
    frames: list[pd.DataFrame] = []
    for path in files:
        frame = pd.read_csv(path)
        frame["_source_file"] = path.name
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _drop_identical_key_duplicates(
    df: pd.DataFrame,
    key_columns: Sequence[str],
    value_columns: Sequence[str],
    label: str,
) -> tuple[pd.DataFrame, int]:
    duplicate_mask = df.duplicated(list(key_columns), keep=False)
    if not duplicate_mask.any():
        return df, 0

    duplicates = df.loc[duplicate_mask, [*key_columns, *value_columns]].copy()
    conflicts: list[tuple[object, ...]] = []
    for key, group in duplicates.groupby(list(key_columns), sort=False, dropna=False):
        if group[list(value_columns)].drop_duplicates().shape[0] != 1:
            conflicts.append(key if isinstance(key, tuple) else (key,))
            if len(conflicts) >= 5:
                break
    if conflicts:
        raise ValueError(f"{label} contains conflicting duplicate keys; examples: {conflicts}")

    before = len(df)
    clean = df.drop_duplicates(list(key_columns), keep="first").copy()
    return clean, before - len(clean)


def load_inputs(paths: AnalysisPaths) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    empirical_files = _discover_files(paths.raw_dir, "empirical_case*.csv")
    null_files = _discover_files(paths.raw_dir, "null_case*.csv")
    manifest_files = _discover_files(paths.manifest_dir, "seed_manifest_case*.csv")

    empirical = _read_csvs(empirical_files, "empirical_case*.csv")
    null = _read_csvs(null_files, "null_case*.csv")
    manifests = _read_csvs(manifest_files, "seed_manifest_case*.csv") if manifest_files else pd.DataFrame()

    _require_columns(empirical, EMPIRICAL_REQUIRED, "empirical inputs")
    _require_columns(null, NULL_REQUIRED, "null inputs")

    empirical["phase5"] = empirical["phase5"].map(_normalise_stage)
    null["phase5"] = null["phase5"].map(_normalise_stage)
    null["null_model"] = null["null_model"].astype(str).str.strip()

    null["iter"] = pd.to_numeric(null["iter"], errors="raise").astype(np.int64)
    for column in ["n_windows", *METRICS]:
        empirical[column] = pd.to_numeric(empirical[column], errors="raise")
        null[column] = pd.to_numeric(null[column], errors="raise")

    for label, frame in [("empirical", empirical), ("null", null)]:
        values = frame[["n_windows", *METRICS]].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            bad = np.argwhere(~np.isfinite(values))[0]
            raise ValueError(f"{label} inputs contain a non-finite numeric value at row {bad[0]}")
        if (frame["n_windows"] <= 0).any():
            raise ValueError(f"{label} inputs contain n_windows <= 0")

    empirical, empirical_duplicate_rows = _drop_identical_key_duplicates(
        empirical,
        ["case_id", "patient_id", "phase5"],
        ["n_windows", *METRICS],
        "empirical inputs",
    )
    null, null_duplicate_rows = _drop_identical_key_duplicates(
        null,
        ["null_model", "iter", "case_id", "patient_id", "phase5"],
        ["n_windows", *METRICS],
        "null inputs",
    )

    unknown_models = sorted(set(null["null_model"]) - set(NULL_MODELS))
    if unknown_models:
        raise ValueError(f"Unexpected null_model values: {unknown_models}")

    load_audit: dict[str, object] = {
        "empirical_files": len(empirical_files),
        "null_files": len(null_files),
        "manifest_files": len(manifest_files),
        "empirical_identical_duplicate_rows_dropped": empirical_duplicate_rows,
        "null_identical_duplicate_rows_dropped": null_duplicate_rows,
    }
    empirical = empirical.sort_values(["case_id", "phase5"]).reset_index(drop=True)
    null = null.sort_values(["null_model", "iter", "case_id", "phase5"]).reset_index(drop=True)
    return empirical, null, manifests, load_audit


def audit_coverage(
    empirical: pd.DataFrame,
    null: pd.DataFrame,
    expected_iterations: int,
) -> pd.DataFrame:
    empirical_keys = empirical[
        ["case_id", "patient_id", "phase5", "n_windows"]
    ].drop_duplicates()
    key_columns = ["case_id", "patient_id", "phase5"]
    rows: list[dict[str, object]] = []
    for model in NULL_MODELS:
        model_df = null[null["null_model"] == model]
        iterations = np.sort(model_df["iter"].unique())
        expected_ids = np.arange(1, expected_iterations + 1, dtype=np.int64)
        ids_exact = np.array_equal(iterations, expected_ids)

        matched = model_df.merge(
            empirical_keys.rename(columns={"n_windows": "n_windows_empirical"}).assign(
                _expected_key=1
            ),
            on=key_columns,
            how="left",
        )
        unknown_key_rows = int(matched["_expected_key"].isna().sum())
        matched_rows = int(matched["_expected_key"].notna().sum())
        expected_rows = int(expected_iterations * len(empirical_keys))
        missing_key_rows = max(0, expected_rows - matched_rows)

        per_iteration = model_df.groupby("iter", as_index=False).size()
        complete_iteration_rows = bool(
            len(per_iteration) == expected_iterations
            and (per_iteration["size"] == len(empirical_keys)).all()
        )
        known_rows = matched["_expected_key"].notna()
        n_windows_match = bool(
            known_rows.any()
            and np.array_equal(
                matched.loc[known_rows, "n_windows"].to_numpy(dtype=float),
                matched.loc[known_rows, "n_windows_empirical"].to_numpy(dtype=float),
            )
        )
        rows.append(
            {
                "null_model": model,
                "expected_iterations": expected_iterations,
                "observed_iterations": int(len(iterations)),
                "iteration_min": int(iterations.min()) if len(iterations) else np.nan,
                "iteration_max": int(iterations.max()) if len(iterations) else np.nan,
                "iteration_ids_exact_1_to_B": ids_exact,
                "expected_case_stage_rows": expected_rows,
                "observed_case_stage_rows": int(len(model_df)),
                "unknown_case_stage_key_rows": unknown_key_rows,
                "missing_case_stage_key_rows": missing_key_rows,
                "every_iteration_has_complete_case_stage_rows": complete_iteration_rows,
                "null_n_windows_match_empirical_by_key": n_windows_match,
                "coverage_pass": bool(
                    ids_exact
                    and unknown_key_rows == 0
                    and missing_key_rows == 0
                    and complete_iteration_rows
                    and n_windows_match
                ),
            }
        )
    return pd.DataFrame(rows)


def aggregate_patient_stages(empirical: pd.DataFrame, null: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    empirical_stage = (
        empirical.groupby(["patient_id", "phase5"], as_index=False)
        .agg(
            n_seizures=("case_id", "nunique"),
            **{metric: (metric, "mean") for metric in METRICS},
        )
        .rename(columns={"phase5": "stage"})
    )
    null_stage = (
        null.groupby(["null_model", "iter", "patient_id", "phase5"], as_index=False)
        .agg(
            n_seizures=("case_id", "nunique"),
            **{metric: (metric, "mean") for metric in METRICS},
        )
        .rename(columns={"phase5": "stage"})
    )
    return empirical_stage, null_stage


def empirical_patient_differences(empirical_stage: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for metric in METRICS:
        pivot = empirical_stage.pivot(index="patient_id", columns="stage", values=metric)
        if "pre-ictal" not in pivot.columns:
            raise ValueError("Empirical patient-stage data have no pre-ictal baseline")
        for target in TARGET_STAGES:
            if target not in pivot.columns:
                raise ValueError(f"Empirical patient-stage data have no {target} stage")
            diff = (pivot[target] - pivot["pre-ictal"]).dropna()
            rows.append(
                pd.DataFrame(
                    {
                        "patient_id": diff.index.astype(str),
                        "metric": metric,
                        "target_stage": target,
                        "empirical_diff": diff.to_numpy(dtype=float),
                    }
                )
            )
    return pd.concat(rows, ignore_index=True)


def _patient_null_matrix(
    null_stage: pd.DataFrame,
    model: str,
    metric: str,
    target: str,
    patients: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    subset = null_stage[null_stage["null_model"] == model]
    pivot = subset.pivot(index=["iter", "patient_id"], columns="stage", values=metric)
    if "pre-ictal" not in pivot.columns or target not in pivot.columns:
        raise ValueError(f"Null patient-stage data lack baseline or {target}: {model}, {metric}")
    diff = (pivot[target] - pivot["pre-ictal"]).rename("null_diff").reset_index()
    matrix = diff.pivot(index="iter", columns="patient_id", values="null_diff").sort_index()
    missing_patients = sorted(set(patients) - set(matrix.columns.astype(str)))
    if missing_patients:
        raise ValueError(f"Null matrix lacks patients for {model}/{metric}/{target}: {missing_patients}")
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.loc[:, list(patients)]
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Null patient-difference matrix is incomplete for {model}/{metric}/{target}")
    return matrix.index.to_numpy(dtype=np.int64), values


def _bootstrap_mean_ci(values: np.ndarray, reps: int, rng: np.random.Generator) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.nan, np.nan
    output = np.empty(reps, dtype=float)
    chunk_size = 2_000
    for start in range(0, reps, chunk_size):
        stop = min(start + chunk_size, reps)
        indices = rng.integers(0, values.size, size=(stop - start, values.size))
        output[start:stop] = values[indices].mean(axis=1)
    low, high = np.quantile(output, [0.025, 0.975])
    return float(low), float(high)


def _combined_patient_null_bootstrap_ci(
    empirical_values: np.ndarray,
    null_matrix: np.ndarray,
    reps: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Combine a patient cluster bootstrap with finite-null-mean uncertainty.

    Patients are resampled non-parametrically.  Conditional on each patient
    resample, finite-B uncertainty of the null mean is added with a Gaussian
    multiplier bootstrap using the empirical between-iteration covariance.
    At production B>=1000 this is the computationally stable large-B limit of
    explicitly resampling all iteration rows in every patient replicate, while
    retaining cross-patient covariance among null effects.
    """

    empirical_values = np.asarray(empirical_values, dtype=float)
    null_matrix = np.asarray(null_matrix, dtype=float)
    n_iterations, n_patients = null_matrix.shape
    if empirical_values.shape != (n_patients,):
        raise ValueError("Empirical patient vector and null patient matrix are misaligned")

    output = np.empty(reps, dtype=float)
    patient_levels = np.arange(n_patients)
    null_patient_mean = null_matrix.mean(axis=0)
    null_covariance = np.atleast_2d(np.cov(null_matrix, rowvar=False, ddof=1))
    residual_patient_mean = empirical_values - null_patient_mean
    chunk_size = 2_000
    for start in range(0, reps, chunk_size):
        stop = min(start + chunk_size, reps)
        width = stop - start
        sampled_patients = rng.integers(0, n_patients, size=(width, n_patients))
        weights = (sampled_patients[:, :, None] == patient_levels[None, None, :]).sum(axis=1)
        weights = weights.astype(float) / float(n_patients)  # width x patients
        patient_boot = weights @ residual_patient_mean
        iteration_variance = np.einsum(
            "bi,ij,bj->b", weights, null_covariance, weights, optimize=True
        )
        iteration_standard_error = np.sqrt(
            np.clip(iteration_variance, 0.0, None) / float(n_iterations)
        )
        output[start:stop] = patient_boot + rng.normal(size=width) * iteration_standard_error

    low, high = np.quantile(output, [0.025, 0.975])
    return float(low), float(high)


def _monte_carlo_p(null_values: np.ndarray, observed: float) -> dict[str, float]:
    values = np.asarray(null_values, dtype=float)
    b = len(values)
    upper = (float(np.count_nonzero(values >= observed)) + 1.0) / (b + 1.0)
    lower = (float(np.count_nonzero(values <= observed)) + 1.0) / (b + 1.0)
    smaller_tail = min(upper, lower)
    two_sided = min(1.0, 2.0 * smaller_tail)
    denominator = float(b + 1)
    return {
        "upper_tail_p": upper,
        "lower_tail_p": lower,
        "two_sided_p": two_sided,
        "upper_tail_mc_se": math.sqrt(upper * (1.0 - upper) / denominator),
        "lower_tail_mc_se": math.sqrt(lower * (1.0 - lower) / denominator),
        # The reported two-sided value is 2*min(one-sided tails), so its MC
        # standard error follows the selected one-sided tail rather than a
        # separate binomial estimate based on the clipped two-sided value.
        "two_sided_mc_se": min(
            1.0,
            2.0 * math.sqrt(smaller_tail * (1.0 - smaller_tail) / denominator),
        ),
    }


def _bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy(), kind="mergesort")
    ordered = valid.to_numpy()[order]
    m = len(ordered)
    adjusted = ordered * m / np.arange(1, m + 1, dtype=float)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    ordered_indices = valid.index.to_numpy()[order]
    result.loc[ordered_indices] = adjusted
    return result


def _convergence_checkpoints(b: int) -> list[int]:
    candidates = [10, 20, 50, 100, 250, 500, 1_000, 2_000, 3_000, 4_000, 5_000]
    if b >= 20:
        candidates.append(b // 2)
    checkpoints = [value for value in candidates if value <= b]
    if b not in checkpoints:
        checkpoints.append(b)
    return sorted(set(checkpoints))


def build_statistics(
    empirical_stage: pd.DataFrame,
    null_stage: pd.DataFrame,
    iteration_hashes: pd.DataFrame,
    bootstrap_reps: int,
    bootstrap_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    empirical_diffs = empirical_patient_differences(empirical_stage)
    summary_rows: list[dict[str, object]] = []
    null_effect_rows: list[pd.DataFrame] = []
    patient_effect_rows: list[pd.DataFrame] = []
    convergence_rows: list[dict[str, object]] = []

    unique_by_model = iteration_hashes.groupby("null_model")["signature_sha256"].nunique().to_dict()

    for metric in METRICS:
        for target in TARGET_STAGES:
            emp = empirical_diffs[
                (empirical_diffs["metric"] == metric)
                & (empirical_diffs["target_stage"] == target)
            ].sort_values("patient_id")
            patients = emp["patient_id"].astype(str).tolist()
            empirical_values = emp["empirical_diff"].to_numpy(dtype=float)
            empirical_effect = float(empirical_values.mean())
            empirical_ci_rng = np.random.default_rng(
                _stable_seed(bootstrap_seed, "empirical", metric, target)
            )
            empirical_ci_low, empirical_ci_high = _bootstrap_mean_ci(
                empirical_values, bootstrap_reps, empirical_ci_rng
            )

            for model in NULL_MODELS:
                iteration_ids, null_matrix = _patient_null_matrix(
                    null_stage, model, metric, target, patients
                )
                null_values = null_matrix.mean(axis=1)
                null_mean = float(null_values.mean())
                null_ci_low, null_ci_high = np.quantile(null_values, [0.025, 0.975])
                empirical_minus_null = empirical_effect - null_mean
                status = (
                    "invariant_by_construction"
                    if model == "orientation_pair_swap" and metric == "S"
                    else "valid"
                )
                if status == "invariant_by_construction":
                    # Enforce the exact scientific semantics of this control:
                    # pairwise orientation swapping preserves S by definition.
                    null_mean = empirical_effect
                    null_ci_low = empirical_effect
                    null_ci_high = empirical_effect
                    empirical_minus_null = 0.0
                    combined_ci_low = 0.0
                    combined_ci_high = 0.0
                    p_fields = {
                        key: np.nan
                        for key in _monte_carlo_p(null_values, empirical_effect)
                    }
                else:
                    combined_rng = np.random.default_rng(
                        _stable_seed(bootstrap_seed, "combined", model, metric, target)
                    )
                    combined_ci_low, combined_ci_high = _combined_patient_null_bootstrap_ci(
                        empirical_values, null_matrix, bootstrap_reps, combined_rng
                    )
                    p_fields = _monte_carlo_p(null_values, empirical_effect)

                row: dict[str, object] = {
                    "null_model": model,
                    "metric": metric,
                    "target_stage": target,
                    "inferential_status": status,
                    "n_patients": len(patients),
                    "n_null_iterations": len(iteration_ids),
                    "n_unique_null_iterations": int(unique_by_model.get(model, 0)),
                    "empirical_effect": empirical_effect,
                    "empirical_patient_bootstrap_ci_low": empirical_ci_low,
                    "empirical_patient_bootstrap_ci_high": empirical_ci_high,
                    "null_mean": null_mean,
                    "null_mean_mc_se": float(
                        np.std(null_values, ddof=1) / math.sqrt(len(null_values))
                    ),
                    "null_interval_low": float(null_ci_low),
                    "null_interval_high": float(null_ci_high),
                    "empirical_minus_null": empirical_minus_null,
                    "combined_bootstrap_ci_low": combined_ci_low,
                    "combined_bootstrap_ci_high": combined_ci_high,
                    **p_fields,
                    "bh_adjusted_p": np.nan,
                    "n_patient_bootstrap": bootstrap_reps,
                    "bootstrap_seed": bootstrap_seed,
                    "empirical_bootstrap_method": "patient_cluster_percentile",
                    "combined_bootstrap_method": "patient_cluster_plus_gaussian_multiplier_finite_null_mean",
                }
                summary_rows.append(row)

                null_effect_rows.append(
                    pd.DataFrame(
                        {
                            "null_model": model,
                            "iter": iteration_ids,
                            "metric": metric,
                            "target_stage": target,
                            "null_effect": null_values,
                            "empirical_effect": empirical_effect,
                        }
                    )
                )

                patient_null_means = null_matrix.mean(axis=0)
                if status == "invariant_by_construction":
                    patient_null_means = empirical_values.copy()
                patient_effect_rows.append(
                    pd.DataFrame(
                        {
                            "null_model": model,
                            "patient_id": patients,
                            "metric": metric,
                            "target_stage": target,
                            "empirical_diff": empirical_values,
                            "patient_null_mean_diff": patient_null_means,
                            "patient_empirical_minus_null": empirical_values - patient_null_means,
                            "inferential_status": status,
                        }
                    )
                )

                for checkpoint in _convergence_checkpoints(len(null_values)):
                    prefix = null_values[:checkpoint]
                    prefix_mean = float(prefix.mean())
                    prefix_ci = np.quantile(prefix, [0.025, 0.975])
                    prefix_p = _monte_carlo_p(prefix, empirical_effect)
                    if status != "valid":
                        prefix_p = {key: np.nan for key in prefix_p}
                    convergence_rows.append(
                        {
                            "null_model": model,
                            "metric": metric,
                            "target_stage": target,
                            "inferential_status": status,
                            "checkpoint_B": checkpoint,
                            "null_mean": prefix_mean,
                            "null_mean_mc_se": float(
                                np.std(prefix, ddof=1) / math.sqrt(checkpoint)
                            ) if checkpoint > 1 else np.nan,
                            "null_interval_low": float(prefix_ci[0]),
                            "null_interval_high": float(prefix_ci[1]),
                            "empirical_minus_null": empirical_effect - prefix_mean,
                            "bh_adjusted_p": np.nan,
                            **prefix_p,
                        }
                    )

    summary = pd.DataFrame(summary_rows)
    valid_mask = summary["inferential_status"].eq("valid")
    summary.loc[valid_mask, "bh_adjusted_p"] = _bh_adjust(
        summary.loc[valid_mask, "two_sided_p"]
    )
    summary = summary.sort_values(["null_model", "metric", "target_stage"]).reset_index(drop=True)

    null_effects = pd.concat(null_effect_rows, ignore_index=True)
    patient_effects = pd.concat(patient_effect_rows, ignore_index=True)
    convergence = pd.DataFrame(convergence_rows)
    valid_convergence = convergence["inferential_status"].eq("valid")
    for checkpoint, indices in convergence.loc[valid_convergence].groupby(
        "checkpoint_B"
    ).groups.items():
        convergence.loc[indices, "bh_adjusted_p"] = _bh_adjust(
            convergence.loc[indices, "two_sided_p"]
        )
    return summary, null_effects, patient_effects, convergence


def build_convergence_acceptance(
    convergence: pd.DataFrame,
    final_b: int,
) -> pd.DataFrame:
    """Audit prefix-B stability without treating an integrity PASS as convergence.

    The reference prefix is the first half of the final sequence (B/2).  The
    gate requires stable effect direction and 0.05 decisions, drift in the null
    mean and raw Monte Carlo p within three final-run Monte Carlo standard
    errors, central-null-interval endpoint drift within 10% of the final
    interval width, and final raw/BH values more than one raw-p Monte Carlo
    standard error from 0.05.  The thresholds are deliberately explicit so a
    marginal B=1000 result triggers extension rather than silent acceptance.
    """

    reference_b = final_b // 2
    valid = convergence[convergence["inferential_status"].eq("valid")].copy()
    keys = ["null_model", "metric", "target_stage"]
    final_rows = valid[valid["checkpoint_B"].eq(final_b)].set_index(keys)
    reference_rows = valid[valid["checkpoint_B"].eq(reference_b)].set_index(keys)
    expected_valid_tests = len(NULL_MODELS) * len(METRICS) * len(TARGET_STAGES) - len(TARGET_STAGES)

    if (
        final_b < CONVERGENCE_MIN_FINAL_B
        or len(final_rows) != expected_valid_tests
        or len(reference_rows) != expected_valid_tests
        or set(final_rows.index) != set(reference_rows.index)
    ):
        return pd.DataFrame(
            [
                {
                    "null_model": "__overall__",
                    "metric": "__coverage__",
                    "target_stage": "__coverage__",
                    "reference_B": reference_b,
                    "final_B": final_b,
                    "convergence_pass": False,
                    "extension_recommended": True,
                    "note": (
                        f"Convergence requires B >= {CONVERGENCE_MIN_FINAL_B} and exact "
                        f"coverage of {expected_valid_tests} valid tests at B/2 and B."
                    ),
                }
            ]
        )

    rows: list[dict[str, object]] = []
    two_sided_floor = 2.0 / float(final_b + 1)
    for key in sorted(final_rows.index):
        final = final_rows.loc[key]
        reference = reference_rows.loc[key]
        mean_drift = abs(float(reference["null_mean"]) - float(final["null_mean"]))
        mean_se = float(final["null_mean_mc_se"])
        mean_drift_z = mean_drift / mean_se if mean_se > 0 else (0.0 if mean_drift == 0 else np.inf)

        p_drift = abs(float(reference["two_sided_p"]) - float(final["two_sided_p"]))
        p_se = float(final["two_sided_mc_se"])
        p_drift_scale = max(p_se, two_sided_floor)
        p_drift_z = p_drift / p_drift_scale

        low_drift = abs(float(reference["null_interval_low"]) - float(final["null_interval_low"]))
        high_drift = abs(float(reference["null_interval_high"]) - float(final["null_interval_high"]))
        final_interval_width = max(
            abs(float(final["null_interval_high"]) - float(final["null_interval_low"])),
            np.finfo(float).eps,
        )
        endpoint_drift_fraction = max(low_drift, high_drift) / final_interval_width

        direction_stable = bool(
            np.sign(float(reference["empirical_minus_null"]))
            == np.sign(float(final["empirical_minus_null"]))
        )
        raw_decision_stable = bool(
            (float(reference["two_sided_p"]) < CONVERGENCE_ALPHA)
            == (float(final["two_sided_p"]) < CONVERGENCE_ALPHA)
        )
        bh_decision_stable = bool(
            (float(reference["bh_adjusted_p"]) < CONVERGENCE_ALPHA)
            == (float(final["bh_adjusted_p"]) < CONVERGENCE_ALPHA)
        )
        raw_boundary_clear = abs(float(final["two_sided_p"]) - CONVERGENCE_ALPHA) > p_se
        bh_boundary_clear = abs(float(final["bh_adjusted_p"]) - CONVERGENCE_ALPHA) > p_se

        passed = bool(
            direction_stable
            and raw_decision_stable
            and bh_decision_stable
            and mean_drift_z <= CONVERGENCE_MAX_MEAN_DRIFT_Z
            and p_drift_z <= CONVERGENCE_MAX_P_DRIFT_Z
            and endpoint_drift_fraction <= CONVERGENCE_MAX_INTERVAL_ENDPOINT_FRACTION
            and raw_boundary_clear
            and bh_boundary_clear
        )
        rows.append(
            {
                "null_model": key[0],
                "metric": key[1],
                "target_stage": key[2],
                "reference_B": reference_b,
                "final_B": final_b,
                "reference_null_mean": float(reference["null_mean"]),
                "final_null_mean": float(final["null_mean"]),
                "null_mean_abs_drift": mean_drift,
                "final_null_mean_mc_se": mean_se,
                "null_mean_drift_z": mean_drift_z,
                "reference_null_interval_low": float(reference["null_interval_low"]),
                "final_null_interval_low": float(final["null_interval_low"]),
                "reference_null_interval_high": float(reference["null_interval_high"]),
                "final_null_interval_high": float(final["null_interval_high"]),
                "max_interval_endpoint_drift_fraction": endpoint_drift_fraction,
                "reference_two_sided_p": float(reference["two_sided_p"]),
                "final_two_sided_p": float(final["two_sided_p"]),
                "final_two_sided_mc_se": p_se,
                "two_sided_p_abs_drift": p_drift,
                "two_sided_p_drift_z": p_drift_z,
                "reference_bh_adjusted_p": float(reference["bh_adjusted_p"]),
                "final_bh_adjusted_p": float(final["bh_adjusted_p"]),
                "effect_direction_stable": direction_stable,
                "raw_0.05_decision_stable": raw_decision_stable,
                "bh_0.05_decision_stable": bh_decision_stable,
                "raw_final_clear_of_0.05_by_one_mc_se": raw_boundary_clear,
                "bh_final_clear_of_0.05_by_one_raw_p_mc_se": bh_boundary_clear,
                "convergence_pass": passed,
                "extension_recommended": not passed,
                "note": "B/2-to-B prefix stability audit; intervals are central null intervals",
            }
        )
    return pd.DataFrame(rows)


def build_iteration_hash_audit(null: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    signature_columns = ["case_id", "patient_id", "phase5", "n_windows", *METRICS]
    records: list[dict[str, object]] = []
    for (model, iteration), group in null.groupby(["null_model", "iter"], sort=True):
        canonical = group.sort_values(["case_id", "patient_id", "phase5"])[signature_columns]
        row_hashes = pd.util.hash_pandas_object(canonical, index=False, categorize=False)
        digest = hashlib.sha256(row_hashes.to_numpy(dtype=np.uint64).tobytes()).hexdigest()
        records.append(
            {
                "null_model": model,
                "iter": int(iteration),
                "signature_sha256": digest,
                "n_case_stage_rows": len(canonical),
            }
        )
    hashes = pd.DataFrame(records)
    hashes["first_iter_with_signature"] = hashes.groupby(
        ["null_model", "signature_sha256"]
    )["iter"].transform("min")
    hashes["is_duplicate_signature"] = hashes["iter"] != hashes["first_iter_with_signature"]

    rows: list[dict[str, object]] = []
    for model in NULL_MODELS:
        subset = hashes[hashes["null_model"] == model]
        n_observed = int(subset["iter"].nunique())
        n_unique = int(subset["signature_sha256"].nunique())
        rows.append(
            {
                "null_model": model,
                "n_observed_iterations": n_observed,
                "n_unique_iteration_signatures": n_unique,
                "n_duplicate_iteration_signatures": n_observed - n_unique,
                "all_iteration_signatures_unique": n_observed == n_unique,
            }
        )
    return hashes, pd.DataFrame(rows)


def audit_phase9_empirical(empirical: pd.DataFrame, phase9_path: Path) -> pd.DataFrame:
    if not phase9_path.is_file():
        return pd.DataFrame(
            [
                {
                    "metric": "__coverage__",
                    "n_empirical_case_stage_rows": len(empirical),
                    "n_phase9_case_stage_rows": 0,
                    "n_matched_rows": 0,
                    "max_abs_diff": np.nan,
                    "mean_abs_diff": np.nan,
                    "rmse": np.nan,
                    "empirical_rows_covered_by_phase9": False,
                    "numerical_match_pass": False,
                    "full_phase9_coverage_pass": False,
                    "qa_pass": False,
                    "note": f"Missing Phase 9 file: {phase9_path}",
                }
            ]
        )

    phase9 = pd.read_csv(phase9_path)
    _require_columns(phase9, ["case_id", "patient_id", "phase5", *METRICS], "Phase 9 QA input")
    phase9["phase5"] = phase9["phase5"].map(_normalise_stage)
    phase9_case = (
        phase9.groupby(["case_id", "patient_id", "phase5"], as_index=False)
        .agg(n_windows_phase9=("phase5", "size"), **{metric: (metric, "mean") for metric in METRICS})
    )
    merged = empirical.merge(
        phase9_case,
        on=["case_id", "patient_id", "phase5"],
        how="outer",
        indicator=True,
        suffixes=("_v2", "_phase9"),
    )
    n_matched = int((merged["_merge"] == "both").sum())
    window_counts_match = bool(
        (merged.loc[merged["_merge"] == "both", "n_windows"]
         == merged.loc[merged["_merge"] == "both", "n_windows_phase9"]).all()
    )
    empirical_rows_covered = bool(n_matched == len(empirical) and window_counts_match)
    full_coverage_pass = bool(
        empirical_rows_covered and n_matched == len(phase9_case)
    )
    rows: list[dict[str, object]] = []
    for metric in METRICS:
        both = merged[merged["_merge"] == "both"]
        diff = both[f"{metric}_v2"].to_numpy(dtype=float) - both[f"{metric}_phase9"].to_numpy(dtype=float)
        abs_diff = np.abs(diff)
        max_abs = float(abs_diff.max()) if len(abs_diff) else np.nan
        numerical_match = bool(
            empirical_rows_covered
            and np.isfinite(max_abs)
            and max_abs <= FLOAT_TOLERANCE
        )
        rows.append(
            {
                "metric": metric,
                "n_empirical_case_stage_rows": len(empirical),
                "n_phase9_case_stage_rows": len(phase9_case),
                "n_matched_rows": n_matched,
                "max_abs_diff": max_abs,
                "mean_abs_diff": float(abs_diff.mean()) if len(abs_diff) else np.nan,
                "rmse": float(np.sqrt(np.mean(diff**2))) if len(diff) else np.nan,
                "empirical_rows_covered_by_phase9": empirical_rows_covered,
                "numerical_match_pass": numerical_match,
                "full_phase9_coverage_pass": full_coverage_pass,
                "qa_pass": bool(numerical_match and full_coverage_pass),
                "note": "case-stage mean recomputation against Phase 9 window metrics",
            }
        )
    return pd.DataFrame(rows)


def audit_orientation_s(empirical: pd.DataFrame, null: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    orientation = null[null["null_model"] == "orientation_pair_swap"]
    merged = orientation.merge(
        empirical[["case_id", "patient_id", "phase5", "S"]],
        on=["case_id", "patient_id", "phase5"],
        how="left",
        suffixes=("_null", "_empirical"),
    )
    merged["abs_diff"] = (merged["S_null"] - merged["S_empirical"]).abs()
    per_iteration = (
        merged.groupby("iter", as_index=False)
        .agg(
            n_case_stage_rows=("abs_diff", "size"),
            n_missing_empirical=("S_empirical", lambda values: int(values.isna().sum())),
            max_abs_diff=("abs_diff", "max"),
            mean_abs_diff=("abs_diff", "mean"),
        )
    )
    per_iteration["invariance_pass"] = (
        per_iteration["n_missing_empirical"].eq(0)
        & per_iteration["max_abs_diff"].le(ORIENTATION_S_TOLERANCE)
    )
    overall = pd.DataFrame(
        [
            {
                "null_model": "orientation_pair_swap",
                "metric": "S",
                "n_iterations": int(per_iteration["iter"].nunique()),
                "n_case_stage_comparisons": int(len(merged)),
                "max_abs_diff": float(merged["abs_diff"].max()) if len(merged) else np.nan,
                "mean_abs_diff": float(merged["abs_diff"].mean()) if len(merged) else np.nan,
                "tolerance": ORIENTATION_S_TOLERANCE,
                "invariance_pass": bool(len(per_iteration) > 0 and per_iteration["invariance_pass"].all()),
                "interpretation": "invariant_by_construction",
            }
        ]
    )
    return per_iteration, overall


def audit_seed_manifest(
    manifests: pd.DataFrame,
    expected_iterations: int,
    empirical: pd.DataFrame,
    expected_master_seed: int,
) -> pd.DataFrame:
    if manifests.empty:
        return pd.DataFrame(
            [
                {
                    "manifest_rows": 0,
                    "manifest_files_present": False,
                    "required_columns_present": False,
                    "unique_stream_indices": False,
                    "observed_stream_indices_within_declared_range": False,
                    "stream_indices_exact_1_to_num_streams": False,
                    "num_streams_consistent": False,
                    "observed_model_case_window_rows_complete": False,
                    "substream_range_matches_B": False,
                    "sequence_position_is_iteration": False,
                    "case_coverage_complete": False,
                    "null_model_coverage_complete": False,
                    "master_seed_consistent": False,
                    "generator_consistent": False,
                    "master_seed_values": "",
                    "generator_values": "",
                    "seed_design_pass": False,
                    "full_stream_coverage_pass": False,
                    "seed_manifest_pass": False,
                    "note": "No seed_manifest_case*.csv files found",
                }
            ]
        )

    required = [
        "null_model",
        "case_id",
        "global_case_index",
        "window_index",
        "phase5",
        "global_window_index",
        "stream_index",
        "num_streams",
        "master_seed",
        "generator",
        "substream_first",
        "substream_last",
        "sequence_position",
    ]
    columns_present = all(column in manifests.columns for column in required)
    if not columns_present:
        missing = [column for column in required if column not in manifests.columns]
        return pd.DataFrame(
            [
                {
                    "manifest_rows": len(manifests),
                    "manifest_files_present": True,
                    "required_columns_present": False,
                    "unique_stream_indices": False,
                    "observed_stream_indices_within_declared_range": False,
                    "stream_indices_exact_1_to_num_streams": False,
                    "num_streams_consistent": False,
                    "observed_model_case_window_rows_complete": False,
                    "substream_range_matches_B": False,
                    "sequence_position_is_iteration": False,
                    "case_coverage_complete": False,
                    "null_model_coverage_complete": False,
                    "master_seed_consistent": False,
                    "generator_consistent": False,
                    "master_seed_values": "",
                    "generator_values": "",
                    "seed_design_pass": False,
                    "full_stream_coverage_pass": False,
                    "seed_manifest_pass": False,
                    "note": f"Missing manifest columns: {missing}",
                }
            ]
        )

    manifests = manifests.copy()
    manifests["case_id"] = manifests["case_id"].astype(str)
    manifests["null_model"] = manifests["null_model"].astype(str).str.strip()
    manifests["phase5"] = manifests["phase5"].map(_normalise_stage)
    numeric_manifest_columns = [
        "global_case_index",
        "window_index",
        "global_window_index",
        "stream_index",
        "num_streams",
        "master_seed",
        "substream_first",
        "substream_last",
    ]
    for column in numeric_manifest_columns:
        manifests[column] = pd.to_numeric(manifests[column], errors="coerce")

    empirical_windows = empirical[
        ["case_id", "phase5", "n_windows"]
    ].copy()
    empirical_windows["case_id"] = empirical_windows["case_id"].astype(str)
    case_window_totals = (
        empirical_windows.groupby("case_id", as_index=True)["n_windows"].sum().astype(int)
    )
    expected_case_ids = set(case_window_totals.index)
    total_global_windows = int(case_window_totals.sum())
    expected_observed_stream_rows = total_global_windows * len(NULL_MODELS)
    expected_num_streams = expected_observed_stream_rows

    unique_streams = not manifests.duplicated(["stream_index"]).any()
    unique_model_case_window_keys = not manifests.duplicated(
        ["null_model", "case_id", "window_index"]
    ).any()
    num_streams_consistent = bool(
        manifests["num_streams"].notna().all()
        and manifests["num_streams"].nunique() == 1
    )
    declared_num_streams = int(manifests["num_streams"].iloc[0]) if num_streams_consistent else -1
    declared_num_streams_exact = declared_num_streams == expected_num_streams
    stream_indices_within_range = bool(
        num_streams_consistent
        and manifests["stream_index"].notna().all()
        and manifests["stream_index"].between(1, declared_num_streams).all()
    )
    stream_indices_exact = bool(
        num_streams_consistent
        and len(manifests) == declared_num_streams
        and np.array_equal(
            np.sort(manifests["stream_index"].dropna().astype(int).unique()),
            np.arange(1, declared_num_streams + 1, dtype=int),
        )
    )
    observed_rows_complete = len(manifests) == int(expected_observed_stream_rows)
    substreams_ok = bool(
        manifests["substream_first"].eq(1).all()
        and manifests["substream_last"].eq(expected_iterations).all()
    )
    sequence_values = (
        manifests["sequence_position"].astype(str).str.strip().str.lower().str.replace(r"\s+", " ", regex=True)
    )
    sequence_ok = bool(sequence_values.eq("substream = iteration").all())
    cases_ok = set(manifests["case_id"].unique()) == expected_case_ids
    models_ok = set(manifests["null_model"].unique()) == set(NULL_MODELS)

    case_mapping = manifests[["case_id", "global_case_index"]].drop_duplicates()
    case_mapping_one_to_one = bool(
        case_mapping["global_case_index"].notna().all()
        and not case_mapping.duplicated("case_id").any()
        and not case_mapping.duplicated("global_case_index").any()
        and np.array_equal(
            np.sort(case_mapping["global_case_index"].astype(int).unique()),
            np.arange(1, len(expected_case_ids) + 1, dtype=int),
        )
    )

    window_indices_exact = True
    if cases_ok and models_ok:
        for (model, case_id), group in manifests.groupby(["null_model", "case_id"]):
            expected_windows = int(case_window_totals.loc[case_id])
            observed_windows = np.sort(group["window_index"].dropna().astype(int).unique())
            if not np.array_equal(observed_windows, np.arange(1, expected_windows + 1, dtype=int)):
                window_indices_exact = False
                break
    else:
        window_indices_exact = False

    phase_counts = (
        manifests.groupby(["null_model", "case_id", "phase5"], as_index=False)
        .size()
        .rename(columns={"size": "observed_windows"})
    )
    expected_phase_counts = pd.MultiIndex.from_product(
        [NULL_MODELS, sorted(expected_case_ids), STAGES],
        names=["null_model", "case_id", "phase5"],
    ).to_frame(index=False)
    expected_phase_counts = expected_phase_counts.merge(
        empirical_windows.rename(columns={"n_windows": "expected_windows"}),
        on=["case_id", "phase5"],
        how="left",
    ).merge(
        phase_counts,
        on=["null_model", "case_id", "phase5"],
        how="left",
    )
    expected_phase_counts[["expected_windows", "observed_windows"]] = expected_phase_counts[
        ["expected_windows", "observed_windows"]
    ].fillna(0)
    phase_counts_exact = bool(
        np.array_equal(
            expected_phase_counts["expected_windows"].to_numpy(dtype=int),
            expected_phase_counts["observed_windows"].to_numpy(dtype=int),
        )
    )

    mapping_formula_exact = False
    if case_mapping_one_to_one and cases_ok and models_ok:
        ordered_mapping = case_mapping.assign(
            global_case_index=case_mapping["global_case_index"].astype(int)
        ).sort_values("global_case_index")
        ordered_mapping["n_windows"] = ordered_mapping["case_id"].map(case_window_totals).astype(int)
        ordered_mapping["global_offset"] = ordered_mapping["n_windows"].cumsum().shift(fill_value=0)
        offset_map = ordered_mapping.set_index("case_id")["global_offset"]
        expected_global_window = (
            manifests["case_id"].map(offset_map).to_numpy(dtype=float)
            + manifests["window_index"].to_numpy(dtype=float)
        )
        model_offset = manifests["null_model"].map(
            {model: index * total_global_windows for index, model in enumerate(NULL_MODELS)}
        ).to_numpy(dtype=float)
        expected_stream = expected_global_window + model_offset
        mapping_formula_exact = bool(
            np.isfinite(expected_global_window).all()
            and np.array_equal(
                manifests["global_window_index"].to_numpy(dtype=float),
                expected_global_window,
            )
            and np.array_equal(
                manifests["stream_index"].to_numpy(dtype=float),
                expected_stream,
            )
        )

    master_seed_consistent = bool(
        manifests["master_seed"].notna().all()
        and manifests["master_seed"].nunique() == 1
        and int(manifests["master_seed"].iloc[0]) == int(expected_master_seed)
    )
    generator_values = manifests["generator"].astype(str).str.strip().str.lower()
    generator_consistent = bool(generator_values.eq("mrg32k3a").all())
    seed_design_pass = bool(
        columns_present
        and unique_streams
        and unique_model_case_window_keys
        and stream_indices_within_range
        and num_streams_consistent
        and declared_num_streams_exact
        and observed_rows_complete
        and substreams_ok
        and sequence_ok
        and cases_ok
        and models_ok
        and case_mapping_one_to_one
        and window_indices_exact
        and phase_counts_exact
        and mapping_formula_exact
        and master_seed_consistent
        and generator_consistent
    )
    full_stream_coverage_pass = bool(seed_design_pass and stream_indices_exact)
    passed = full_stream_coverage_pass
    return pd.DataFrame(
        [
            {
                "manifest_rows": len(manifests),
                "manifest_files_present": True,
                "required_columns_present": columns_present,
                "unique_stream_indices": unique_streams,
                "unique_model_case_window_keys": unique_model_case_window_keys,
                "observed_stream_indices_within_declared_range": stream_indices_within_range,
                "stream_indices_exact_1_to_num_streams": stream_indices_exact,
                "num_streams_consistent": num_streams_consistent,
                "declared_num_streams": declared_num_streams,
                "expected_num_streams": expected_num_streams,
                "declared_num_streams_exact": declared_num_streams_exact,
                "expected_observed_stream_rows": int(expected_observed_stream_rows),
                "observed_model_case_window_rows_complete": observed_rows_complete,
                "substream_range_matches_B": substreams_ok,
                "sequence_position_is_iteration": sequence_ok,
                "case_coverage_complete": cases_ok,
                "null_model_coverage_complete": models_ok,
                "case_global_index_mapping_one_to_one": case_mapping_one_to_one,
                "window_indices_exact_per_model_case": window_indices_exact,
                "phase_window_counts_match_empirical": phase_counts_exact,
                "global_window_and_stream_formula_exact": mapping_formula_exact,
                "expected_master_seed": int(expected_master_seed),
                "master_seed_consistent": master_seed_consistent,
                "generator_consistent": generator_consistent,
                "master_seed_values": ";".join(sorted(map(str, manifests["master_seed"].dropna().unique()))),
                "generator_values": ";".join(sorted(map(str, manifests["generator"].dropna().unique()))),
                "seed_design_pass": seed_design_pass,
                "full_stream_coverage_pass": full_stream_coverage_pass,
                "seed_manifest_pass": passed,
                "note": "one row per model-case-window stream; substream index equals absolute iteration",
            }
        ]
    )


def build_qa_checks(
    load_audit: dict[str, object],
    coverage: pd.DataFrame,
    uniqueness: pd.DataFrame,
    phase9_qa: pd.DataFrame,
    orientation_overall: pd.DataFrame,
    seed_audit: pd.DataFrame,
    convergence_acceptance: pd.DataFrame,
    expected_iterations: int,
) -> pd.DataFrame:
    checks = [
        {
            "check": "no_duplicate_input_keys",
            "pass": bool(
                int(load_audit.get("empirical_identical_duplicate_rows_dropped", -1)) == 0
                and int(load_audit.get("null_identical_duplicate_rows_dropped", -1)) == 0
            ),
            "detail": "Production inputs must contain no duplicate empirical or null key rows",
        },
        {
            "check": "null_coverage_and_iteration_ids",
            "pass": bool(len(coverage) == len(NULL_MODELS) and coverage["coverage_pass"].all()),
            "detail": f"Each model must contain exact iterations 1..{expected_iterations} and complete case-stage rows",
        },
        {
            "check": "all_iteration_signatures_unique",
            "pass": bool(
                len(uniqueness) == len(NULL_MODELS)
                and uniqueness["all_iteration_signatures_unique"].all()
                and uniqueness["n_unique_iteration_signatures"].eq(expected_iterations).all()
            ),
            "detail": f"Observed and unique realization signatures must both equal {expected_iterations}",
        },
        {
            "check": "empirical_recomputation_matches_phase9_numerically",
            "pass": bool(
                len(phase9_qa) == len(METRICS)
                and phase9_qa["numerical_match_pass"].all()
            ),
            "detail": f"All four case-stage metrics must match Phase 9 within {FLOAT_TOLERANCE:g}",
        },
        {
            "check": "empirical_full_phase9_coverage",
            "pass": bool(
                len(phase9_qa) == len(METRICS)
                and phase9_qa["full_phase9_coverage_pass"].all()
            ),
            "detail": "Production empirical output must cover every Phase 9 case-stage row and window count",
        },
        {
            "check": "orientation_pair_swap_preserves_S",
            "pass": bool(
                len(orientation_overall) == 1
                and bool(orientation_overall.iloc[0]["invariance_pass"])
            ),
            "detail": f"S must be invariant by construction within {ORIENTATION_S_TOLERANCE:g}",
        },
        {
            "check": "seed_manifest_design",
            "pass": bool(len(seed_audit) == 1 and bool(seed_audit.iloc[0]["seed_design_pass"])),
            "detail": "Observed model-case-window streams are unique/in-range and use absolute-iteration substreams",
        },
        {
            "check": "seed_manifest_full_stream_coverage",
            "pass": bool(len(seed_audit) == 1 and bool(seed_audit.iloc[0]["full_stream_coverage_pass"])),
            "detail": "Production union must contain every declared stream index exactly once",
        },
        {
            "check": "prefix_B_convergence_acceptance",
            "pass": bool(
                len(convergence_acceptance) == (
                    len(NULL_MODELS) * len(METRICS) * len(TARGET_STAGES) - len(TARGET_STAGES)
                )
                and convergence_acceptance["convergence_pass"].all()
            ),
            "detail": (
                "B/2-to-B stability of null means, central null intervals, raw/BH decisions, "
                "and Monte Carlo uncertainty"
            ),
        },
    ]
    return pd.DataFrame(checks)


def _format_number(value: object) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.6g}"


def write_report(
    paths: AnalysisPaths,
    empirical: pd.DataFrame,
    null: pd.DataFrame,
    load_audit: dict[str, object],
    summary: pd.DataFrame,
    coverage: pd.DataFrame,
    uniqueness: pd.DataFrame,
    phase9_qa: pd.DataFrame,
    orientation_overall: pd.DataFrame,
    seed_audit: pd.DataFrame,
    qa_checks: pd.DataFrame,
) -> None:
    overall_pass = bool(qa_checks["pass"].all())
    lines = [
        "# Phase 11 W-level null v2 analysis and QA",
        "",
        f"- Version tag: `{paths.tag}`",
        f"- Overall QA status: **{'PASS' if overall_pass else 'FAIL'}**",
        f"- Empirical files: {load_audit['empirical_files']}; null files: {load_audit['null_files']}; seed-manifest files: {load_audit['manifest_files']}.",
        f"- Empirical case-stage rows: {len(empirical)}; null case-stage rows: {len(null)}.",
        f"- Seizures: {empirical['case_id'].nunique()}; patients: {empirical['patient_id'].nunique()}.",
        "- Statistical hierarchy: case-stage means (already aggregated from windows) -> patient-stage means -> paired target minus pre-ictal differences -> group mean.",
        "- Patient bootstrap: 20,000 cluster resamples by default. Combined intervals add a Gaussian multiplier bootstrap for finite-B uncertainty of the null mean, using the full cross-patient iteration covariance.",
        "- Monte Carlo tail probabilities use the +1 correction. BH adjustment is across every valid two-sided W-level test in this analysis.",
        "- The reported two-sided Monte Carlo floor is `2 / (B + 1)` because the analysis doubles the smaller +1-corrected one-sided tail.",
        "- `orientation_pair_swap × S` is excluded from inference and labelled `invariant_by_construction`.",
        "",
        "## QA gates",
        "",
        "| Check | Pass | Detail |",
        "|---|---:|---|",
    ]
    for _, row in qa_checks.iterrows():
        lines.append(f"| {row['check']} | {'yes' if row['pass'] else 'no'} | {row['detail']} |")

    lines.extend(["", "## Iteration coverage and uniqueness", ""])
    for _, row in coverage.iterrows():
        unique_row = uniqueness[uniqueness["null_model"] == row["null_model"]].iloc[0]
        lines.append(
            f"- `{row['null_model']}`: observed B={row['observed_iterations']}, "
            f"unique signatures={unique_row['n_unique_iteration_signatures']}, "
            f"duplicate signatures={unique_row['n_duplicate_iteration_signatures']}, "
            f"coverage={'PASS' if row['coverage_pass'] else 'FAIL'}."
        )

    lines.extend(["", "## Independent QA", ""])
    max_phase9 = phase9_qa["max_abs_diff"].max() if "max_abs_diff" in phase9_qa else np.nan
    lines.append(f"- Maximum empirical absolute difference from Phase 9: {_format_number(max_phase9)}.")
    orientation_row = orientation_overall.iloc[0]
    lines.append(
        f"- Maximum `orientation_pair_swap` S difference: {_format_number(orientation_row['max_abs_diff'])}; "
        f"invariance={'PASS' if orientation_row['invariance_pass'] else 'FAIL'}."
    )
    lines.append(
        f"- Seed manifest: {'PASS' if seed_audit.iloc[0]['seed_manifest_pass'] else 'FAIL'} "
        f"({seed_audit.iloc[0]['note']})."
    )

    lines.extend(
        [
            "",
            "## Inferential summary",
            "",
            "| Null model | Metric | Stage | Empirical | Null mean | Empirical-null | Combined 95% CI | two-sided *p* | BH-adjusted *p* | Status |",
            "|---|---|---|---:|---:|---:|---|---:|---:|---|",
        ]
    )
    for _, row in summary.iterrows():
        interval = (
            f"[{_format_number(row['combined_bootstrap_ci_low'])}, "
            f"{_format_number(row['combined_bootstrap_ci_high'])}]"
        )
        lines.append(
            f"| {row['null_model']} | {row['metric']} | {row['target_stage']} | "
            f"{_format_number(row['empirical_effect'])} | {_format_number(row['null_mean'])} | "
            f"{_format_number(row['empirical_minus_null'])} | {interval} | "
            f"{_format_number(row['two_sided_p'])} | {_format_number(row['bh_adjusted_p'])} | "
            f"{row['inferential_status']} |"
        )

    lines.extend(
        [
            "",
            "## Output contract",
            "",
            f"- `{SUMMARY_FILENAME}`: complete inferential summary.",
            f"- `{FIGURE_SOURCE_FILENAME}`: stable Fig. 6 source-data contract.",
            "- `null_iteration_effects.csv`: group-level null effect for every iteration/test.",
            "- `patient_level_effects.csv`: empirical, mean-null, and residual effects by patient.",
            "- `convergence_checkpoints.csv`: prefix-B convergence diagnostics.",
            "- `convergence_acceptance.csv`: explicit B/2-to-B stability gate; a failed row requires extending B before figure rendering.",
            "- `iteration_hashes.csv` and `iteration_uniqueness_audit.csv`: exact realization-signature audit.",
            "- `empirical_phase9_qa.csv`, `orientation_S_invariance_audit.csv`, `seed_manifest_audit.csv`, and `qa_checks.csv`: independent QA evidence.",
            "",
            "No manuscript file was opened for writing or modified by this analysis.",
        ]
    )
    (paths.tag_dir / REPORT_FILENAME).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_analysis(
    project_root: Path,
    tag: str,
    expected_iterations: int = DEFAULT_EXPECTED_ITERATIONS,
    bootstrap_reps: int = DEFAULT_BOOTSTRAP_REPS,
    bootstrap_seed: int = DEFAULT_BOOTSTRAP_SEED,
    expected_master_seed: int = DEFAULT_EXPECTED_MASTER_SEED,
) -> tuple[Path, bool]:
    if not tag or tag in {".", ".."} or Path(tag).name != tag:
        raise ValueError("--tag must be one directory name, not a path")
    if expected_iterations < 2:
        raise ValueError("--expected-iterations must be at least 2")
    if bootstrap_reps < 100:
        raise ValueError("--bootstrap-reps must be at least 100")

    paths = AnalysisPaths(project_root=project_root.resolve(), tag=tag)
    if not paths.tag_dir.is_dir():
        raise FileNotFoundError(f"Versioned tag directory does not exist: {paths.tag_dir}")

    empirical, null, manifests, load_audit = load_inputs(paths)
    coverage = audit_coverage(empirical, null, expected_iterations)
    iteration_hashes, uniqueness = build_iteration_hash_audit(null)
    phase9_qa = audit_phase9_empirical(empirical, paths.phase9)
    orientation_per_iteration, orientation_overall = audit_orientation_s(empirical, null)
    seed_audit = audit_seed_manifest(
        manifests,
        expected_iterations,
        empirical,
        expected_master_seed,
    )

    empirical_stage, null_stage = aggregate_patient_stages(empirical, null)
    summary, null_effects, patient_effects, convergence = build_statistics(
        empirical_stage,
        null_stage,
        iteration_hashes,
        bootstrap_reps,
        bootstrap_seed,
    )
    convergence_acceptance = build_convergence_acceptance(
        convergence,
        expected_iterations,
    )
    qa_checks = build_qa_checks(
        load_audit,
        coverage,
        uniqueness,
        phase9_qa,
        orientation_overall,
        seed_audit,
        convergence_acceptance,
        expected_iterations,
    )

    # All writes are new, versioned derivatives inside <tag>; no legacy file is
    # overwritten and no manuscript-writing code is imported.
    empirical[EMPIRICAL_REQUIRED].to_csv(paths.tag_dir / "empirical_case_stage_consolidated.csv", index=False)
    empirical_stage.to_csv(paths.tag_dir / "empirical_patient_stage.csv", index=False)
    empirical_patient_differences(empirical_stage).to_csv(
        paths.tag_dir / "empirical_patient_differences.csv", index=False
    )
    summary.to_csv(paths.tag_dir / SUMMARY_FILENAME, index=False)

    figure_columns = [
        "null_model",
        "metric",
        "target_stage",
        "inferential_status",
        "n_patients",
        "n_null_iterations",
        "n_unique_null_iterations",
        "empirical_effect",
        "empirical_patient_bootstrap_ci_low",
        "empirical_patient_bootstrap_ci_high",
        "null_mean",
        "null_interval_low",
        "null_interval_high",
        "empirical_minus_null",
        "combined_bootstrap_ci_low",
        "combined_bootstrap_ci_high",
        "two_sided_p",
        "two_sided_mc_se",
        "bh_adjusted_p",
    ]
    figure_source = summary[figure_columns].copy()
    figure_source["stage_order"] = figure_source["target_stage"].map(
        {stage: index + 1 for index, stage in enumerate(TARGET_STAGES)}
    )
    figure_source.to_csv(paths.tag_dir / FIGURE_SOURCE_FILENAME, index=False)

    null_effects.to_csv(paths.tag_dir / "null_iteration_effects.csv", index=False)
    patient_effects.to_csv(paths.tag_dir / "patient_level_effects.csv", index=False)
    convergence.to_csv(paths.tag_dir / "convergence_checkpoints.csv", index=False)
    convergence_acceptance.to_csv(
        paths.tag_dir / "convergence_acceptance.csv", index=False
    )
    iteration_hashes.to_csv(paths.tag_dir / "iteration_hashes.csv", index=False)
    uniqueness.to_csv(paths.tag_dir / "iteration_uniqueness_audit.csv", index=False)
    coverage.to_csv(paths.tag_dir / "coverage_audit.csv", index=False)
    phase9_qa.to_csv(paths.tag_dir / "empirical_phase9_qa.csv", index=False)
    orientation_per_iteration.to_csv(paths.tag_dir / "orientation_S_invariance_by_iteration.csv", index=False)
    orientation_overall.to_csv(paths.tag_dir / "orientation_S_invariance_audit.csv", index=False)
    seed_audit.to_csv(paths.tag_dir / "seed_manifest_audit.csv", index=False)
    qa_checks.to_csv(paths.tag_dir / "qa_checks.csv", index=False)

    write_report(
        paths,
        empirical,
        null,
        load_audit,
        summary,
        coverage,
        uniqueness,
        phase9_qa,
        orientation_overall,
        seed_audit,
        qa_checks,
    )
    return paths.tag_dir, bool(qa_checks["pass"].all())


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True, help="Version directory name under stats/w_level_null_v2")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=DEFAULT_PROJECT_ROOT,
        help=f"Paper 1 root (default: {DEFAULT_PROJECT_ROOT})",
    )
    parser.add_argument(
        "--expected-iterations",
        type=int,
        default=DEFAULT_EXPECTED_ITERATIONS,
        help=f"Required B per null model (default: {DEFAULT_EXPECTED_ITERATIONS})",
    )
    parser.add_argument(
        "--bootstrap-reps",
        type=int,
        default=DEFAULT_BOOTSTRAP_REPS,
        help=f"Patient/joint bootstrap replicates (default: {DEFAULT_BOOTSTRAP_REPS})",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=DEFAULT_BOOTSTRAP_SEED,
        help=f"Reproducible analysis seed (default: {DEFAULT_BOOTSTRAP_SEED})",
    )
    parser.add_argument(
        "--expected-master-seed",
        type=int,
        default=DEFAULT_EXPECTED_MASTER_SEED,
        help=f"Required MATLAB master seed (default: {DEFAULT_EXPECTED_MASTER_SEED})",
    )
    parser.add_argument(
        "--allow-qa-fail",
        action="store_true",
        help="Write outputs but return exit code 0 when a formal QA gate fails (for diagnostics only)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir, passed = run_analysis(
        project_root=args.project_root,
        tag=args.tag,
        expected_iterations=args.expected_iterations,
        bootstrap_reps=args.bootstrap_reps,
        bootstrap_seed=args.bootstrap_seed,
        expected_master_seed=args.expected_master_seed,
    )
    print(f"W-level null v2 analysis outputs: {output_dir}")
    print(f"Formal QA status: {'PASS' if passed else 'FAIL'}")
    if not passed and not args.allow_qa_fail:
        print("One or more formal QA gates failed; inspect qa_checks.csv and the Markdown report.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
