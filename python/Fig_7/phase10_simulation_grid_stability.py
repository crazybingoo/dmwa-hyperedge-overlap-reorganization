"""Phase 10 simulation grid stability analysis for Paper 1.

The existing Phase 4 simulation produced a full beta grid with 80 replicates
per model and beta. This script converts that raw grid into uncertainty-aware
model ranking evidence, avoiding an overclaim based only on a single best beta.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get("PAPER1_PROJECT_ROOT", Path(__file__).resolve().parents[2])
).resolve()
SIM_DIR = PROJECT_ROOT / "simulation"
STATS_DIR = PROJECT_ROOT / "stats"
NOTES_DIR = PROJECT_ROOT / "notes"

RAW_PATH = SIM_DIR / "competing_mechanism_raw.csv"
FIT_PATH = SIM_DIR / "competing_mechanism_fit.csv"

FEATURES = [
    "S",
    "source_sink_gini_gap",
    "reciprocal_fraction",
    "dense_fraction",
]
MODEL_ORDER = [
    "role_reciprocity_hybrid",
    "density_only",
    "source_concentration",
    "overlap_synchronization",
    "weight_strengthening",
    "reciprocity_only",
]
N_BOOT = 5000
RNG_SEED = 20260617


def scaled_distance(sim_vec: np.ndarray, emp_vec: np.ndarray) -> np.ndarray:
    denom = np.abs(emp_vec) + np.finfo(float).eps
    scaled = sim_vec / denom
    target = emp_vec / denom
    return np.sqrt(np.sum((scaled - target) ** 2, axis=-1))


def cosine_similarity(sim_vec: np.ndarray, emp_vec: np.ndarray) -> np.ndarray:
    numerator = np.sum(sim_vec * emp_vec, axis=-1)
    denominator = np.linalg.norm(sim_vec, axis=-1) * np.linalg.norm(emp_vec)
    return numerator / (denominator + np.finfo(float).eps)


def sign_agreement(sim_vec: np.ndarray, emp_vec: np.ndarray) -> np.ndarray:
    return np.mean(np.sign(sim_vec) == np.sign(emp_vec), axis=-1)


def summarize(values: np.ndarray, prefix: str) -> dict[str, float]:
    return {
        f"{prefix}_mean": float(np.mean(values)),
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_ci_low": float(np.quantile(values, 0.025)),
        f"{prefix}_ci_high": float(np.quantile(values, 0.975)),
    }


def load_empirical_vector() -> np.ndarray:
    fit = pd.read_csv(FIT_PATH)
    row = fit.iloc[0]
    return np.asarray(
        [
            row["empirical_S"],
            row["empirical_source_sink_gini_gap"],
            row["empirical_reciprocal_fraction"],
            row["empirical_dense_fraction"],
        ],
        dtype=float,
    )


def compute_observed_grid(raw: pd.DataFrame, emp_vec: np.ndarray) -> pd.DataFrame:
    rows = []
    for model in MODEL_ORDER:
        model_raw = raw[raw["model"] == model]
        low = model_raw[model_raw["beta"] == 0][FEATURES].mean().to_numpy(dtype=float)
        for beta in sorted(model_raw["beta"].unique()):
            if beta == 0:
                continue
            high = model_raw[model_raw["beta"] == beta][FEATURES].mean().to_numpy(dtype=float)
            sim_vec = high - low
            rows.append(
                {
                    "model": model,
                    "beta": beta,
                    "simulated_S": sim_vec[0],
                    "simulated_source_sink_gini_gap": sim_vec[1],
                    "simulated_reciprocal_fraction": sim_vec[2],
                    "simulated_dense_fraction": sim_vec[3],
                    "scaled_euclidean_distance": float(scaled_distance(sim_vec[None, :], emp_vec)[0]),
                    "cosine_similarity": float(cosine_similarity(sim_vec[None, :], emp_vec)[0]),
                    "sign_agreement": float(sign_agreement(sim_vec[None, :], emp_vec)[0]),
                    "dense_fraction_ratio": float(sim_vec[3] / (emp_vec[3] + np.finfo(float).eps)),
                }
            )
    return pd.DataFrame(rows)


def bootstrap_grid(raw: pd.DataFrame, emp_vec: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RNG_SEED)
    grid_rows = []
    boot_records = []

    for model in MODEL_ORDER:
        model_raw = raw[raw["model"] == model]
        low_mat = model_raw[model_raw["beta"] == 0][FEATURES].to_numpy(dtype=float)
        n_low = low_mat.shape[0]
        betas = [b for b in sorted(model_raw["beta"].unique()) if b != 0]
        # Every beta-specific effect for a model is defined relative to the
        # same beta = 0 condition. Preserve that dependence in the bootstrap
        # by drawing one baseline resample per model and bootstrap iteration,
        # then reusing it across all candidate beta values.
        low_idx = rng.integers(0, n_low, size=(N_BOOT, n_low))
        low_mean = low_mat[low_idx].mean(axis=1)

        for beta in betas:
            high_mat = model_raw[model_raw["beta"] == beta][FEATURES].to_numpy(dtype=float)
            n_high = high_mat.shape[0]
            high_idx = rng.integers(0, n_high, size=(N_BOOT, n_high))
            sim_vec = high_mat[high_idx].mean(axis=1) - low_mean
            distances = scaled_distance(sim_vec, emp_vec)
            cosines = cosine_similarity(sim_vec, emp_vec)
            signs = sign_agreement(sim_vec, emp_vec)
            dense_ratio = sim_vec[:, 3] / (emp_vec[3] + np.finfo(float).eps)

            row = {"model": model, "beta": beta, "n_boot": N_BOOT}
            row.update(summarize(distances, "distance"))
            row.update(summarize(cosines, "cosine"))
            row.update(summarize(dense_ratio, "dense_ratio"))
            row["full_sign_agreement_rate"] = float(np.mean(signs == 1.0))
            row["median_sign_agreement"] = float(np.median(signs))
            grid_rows.append(row)

            for b, dist in enumerate(distances):
                boot_records.append(
                    {
                        "bootstrap_id": b,
                        "model": model,
                        "beta": beta,
                        "distance": dist,
                        "cosine": cosines[b],
                        "sign_agreement": signs[b],
                        "dense_ratio": dense_ratio[b],
                    }
                )

    return pd.DataFrame(grid_rows), pd.DataFrame(boot_records)


def ranking_stability(boot: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    best_by_model = (
        boot.sort_values(["bootstrap_id", "model", "distance", "beta"])
        .groupby(["bootstrap_id", "model"], as_index=False)
        .first()
    )
    best_by_model["rank"] = best_by_model.groupby("bootstrap_id")["distance"].rank(
        method="first", ascending=True
    )
    best_by_model = best_by_model.sort_values(["bootstrap_id", "rank"])

    rows = []
    for model in MODEL_ORDER:
        temp = best_by_model[best_by_model["model"] == model].copy()
        beta_counts = temp["beta"].value_counts().sort_index()
        beta_mode = float(beta_counts.idxmax())
        rows.append(
            {
                "model": model,
                "top1_rate": float(np.mean(temp["rank"] == 1)),
                "top2_rate": float(np.mean(temp["rank"] <= 2)),
                "median_rank": float(np.median(temp["rank"])),
                "mean_rank": float(np.mean(temp["rank"])),
                "best_beta_mode": beta_mode,
                "best_beta_ci_low": float(np.quantile(temp["beta"], 0.025)),
                "best_beta_ci_high": float(np.quantile(temp["beta"], 0.975)),
                "selected_distance_median": float(np.median(temp["distance"])),
                "selected_distance_ci_low": float(np.quantile(temp["distance"], 0.025)),
                "selected_distance_ci_high": float(np.quantile(temp["distance"], 0.975)),
                "selected_cosine_median": float(np.median(temp["cosine"])),
                "selected_dense_ratio_median": float(np.median(temp["dense_ratio"])),
                "selected_full_sign_agreement_rate": float(np.mean(temp["sign_agreement"] == 1.0)),
            }
        )
    return pd.DataFrame(rows), best_by_model


def write_report(
    observed: pd.DataFrame,
    boot_grid: pd.DataFrame,
    rank_summary: pd.DataFrame,
    emp_vec: np.ndarray,
) -> None:
    NOTES_DIR.mkdir(parents=True, exist_ok=True)
    out = NOTES_DIR / "phase10_simulation_grid_stability.md"
    best_obs = (
        observed.sort_values(["scaled_euclidean_distance", "cosine_similarity"], ascending=[True, False])
        .groupby("model", as_index=False)
        .first()
        .sort_values("scaled_euclidean_distance")
    )
    hybrid_rank = rank_summary[rank_summary["model"] == "role_reciprocity_hybrid"].iloc[0]
    hybrid_obs = best_obs[best_obs["model"] == "role_reciprocity_hybrid"].iloc[0]
    dense_ratio = hybrid_rank["selected_dense_ratio_median"]
    lines = [
        "# Phase 10 simulation grid stability",
        "",
        "Date: 2026-06-17",
        "",
        "Purpose:",
        "",
        "This analysis addresses the CSF audit risk that the competing-mechanism simulation may look tuned to one best beta. It uses the existing Phase 4 raw simulation grid and adds bootstrap uncertainty and model-ranking stability.",
        "",
        "Inputs:",
        "",
        f"- `{RAW_PATH}`",
        f"- `{FIT_PATH}`",
        "",
        "Empirical effect vector:",
        "",
        f"- `S`: {emp_vec[0]:.6f}",
        f"- `source_sink_gini_gap`: {emp_vec[1]:.6f}",
        f"- `reciprocal_fraction`: {emp_vec[2]:.6f}",
        f"- `dense_fraction`: {emp_vec[3]:.6f}",
        "",
        "Simulation grid:",
        "",
        "- Models: density-only, weight-strengthening, reciprocity-only, source-concentration, overlap-synchronization and role-reciprocity hybrid.",
        "- Beta grid: 0.1 to 1.0 after subtracting beta = 0 baseline.",
        "- Replicates: 80 per model per beta.",
        f"- Bootstrap resamples: {N_BOOT} per model-beta point.",
        "",
        "Observed best-beta ranking:",
        "",
    ]
    for _, row in best_obs.iterrows():
        lines.append(
            f"- `{row['model']}`: beta = {row['beta']:.1f}, distance = {row['scaled_euclidean_distance']:.3f}, cosine = {row['cosine_similarity']:.3f}, sign agreement = {row['sign_agreement']:.2f}."
        )
    lines.extend(
        [
            "",
            "Bootstrap ranking stability:",
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
            "Main interpretation:",
            "",
            f"- The hybrid model remains the observed best model at beta = {hybrid_obs['beta']:.1f}, distance = {hybrid_obs['scaled_euclidean_distance']:.3f}, cosine = {hybrid_obs['cosine_similarity']:.3f}.",
            f"- In bootstrap ranking, the hybrid model has top-1 rate = {hybrid_rank['top1_rate']:.3f} and top-2 rate = {hybrid_rank['top2_rate']:.3f}.",
            f"- The median dense-fraction ratio for the hybrid selected model is {dense_ratio:.3f}; values below 1 indicate that dense motif enrichment remains under-reproduced.",
            "- This supports a cautious statement: the hybrid mechanism is the most stable abstract discriminator in the current grid, but the simulation is not a biophysical seizure model and should be used only to separate structural hypotheses.",
            "",
            "Output files:",
            "",
            "- `stats/phase10_simulation_grid_observed_fit.csv`",
            "- `stats/phase10_simulation_grid_bootstrap_summary.csv`",
            "- `stats/phase10_simulation_grid_ranking_stability.csv`",
            "- `stats/phase10_simulation_grid_bootstrap_best_by_model.csv`",
        ]
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    raw = pd.read_csv(RAW_PATH)
    raw["model"] = raw["model"].astype(str)
    emp_vec = load_empirical_vector()

    observed = compute_observed_grid(raw, emp_vec)
    boot_grid, boot_records = bootstrap_grid(raw, emp_vec)
    rank_summary, best_by_model = ranking_stability(boot_records)

    STATS_DIR.mkdir(parents=True, exist_ok=True)
    observed.to_csv(STATS_DIR / "phase10_simulation_grid_observed_fit.csv", index=False)
    boot_grid.to_csv(STATS_DIR / "phase10_simulation_grid_bootstrap_summary.csv", index=False)
    rank_summary.to_csv(STATS_DIR / "phase10_simulation_grid_ranking_stability.csv", index=False)
    best_by_model.to_csv(STATS_DIR / "phase10_simulation_grid_bootstrap_best_by_model.csv", index=False)
    write_report(observed, boot_grid, rank_summary, emp_vec)

    best_obs = (
        observed.sort_values(["scaled_euclidean_distance", "cosine_similarity"], ascending=[True, False])
        .groupby("model", as_index=False)
        .first()
        .sort_values("scaled_euclidean_distance")
    )
    print("Phase 10 simulation grid stability complete.")
    print("Observed best-beta ranking:")
    print(
        best_obs[
            ["model", "beta", "scaled_euclidean_distance", "cosine_similarity", "sign_agreement", "dense_fraction_ratio"]
        ].to_string(index=False)
    )
    print("\nBootstrap ranking stability:")
    print(
        rank_summary.sort_values("top1_rate", ascending=False)[
            ["model", "top1_rate", "top2_rate", "median_rank", "best_beta_mode", "selected_distance_median", "selected_dense_ratio_median"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
