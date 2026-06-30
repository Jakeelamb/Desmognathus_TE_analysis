"""
Canonical diversity statistics for Desmognathus TE compositions.

Input breakdown tables store per-species percentages, so Simpson diversity is
computed from normalized proportions rather than from count combinatorics.
"""

from __future__ import annotations

import functools
import sys
from pathlib import Path

import numpy as np
import pandas as pd


sys.path.insert(0, str(Path(__file__).parent.parent))
from config import paths  # noqa: E402


INPUT_DATA_DIR = paths.results.data
OUTPUT_DATA_DIR = paths.results.data
FIGURE_DIR = paths.results.figures / "diversity_threshold_plots"

ORDER_CSV = INPUT_DATA_DIR / "dnaPipeTE_order_breakdown.csv"
SUPERFAMILY_CSV = INPUT_DATA_DIR / "dnaPipeTE_superfamily_breakdown.csv"
CANONICAL_ORDER_CSV = OUTPUT_DATA_DIR / "diversity_order_stats.csv"
CANONICAL_SUPERFAMILY_CSV = OUTPUT_DATA_DIR / "diversity_superfamily_stats.csv"
THRESHOLDS = [5.0, 1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0]


def positive_numeric_values(row: pd.Series) -> pd.Series:
    values = pd.to_numeric(row, errors="coerce")
    values = values[np.isfinite(values)]
    return values[values > 0]


def calculate_simpson_diversity(row: pd.Series) -> float:
    """Return Gini-Simpson diversity, 1 - sum(p_i^2)."""
    values = positive_numeric_values(row)
    total = values.sum()
    if values.empty or pd.isna(total) or total <= 0:
        return 0.0
    proportions = values / total
    return float(1.0 - np.sum(proportions * proportions))


def calculate_shannon_diversity(row: pd.Series) -> float:
    values = positive_numeric_values(row)
    total = values.sum()
    if values.empty or pd.isna(total) or total <= 0:
        return 0.0
    proportions = values / total
    return float(-np.sum(proportions * np.log(proportions)))


def calculate_pielou_evenness(row: pd.Series) -> float:
    values = positive_numeric_values(row)
    num_categories = len(values)
    if num_categories <= 1:
        return 0.0
    shannon_index = calculate_shannon_diversity(values)
    if pd.isna(shannon_index):
        return 0.0
    return float(shannon_index / np.log(num_categories))


def calculate_diversity_for_threshold(df: pd.DataFrame, threshold: float) -> tuple[pd.DataFrame, int]:
    print(f"--- Calculating diversity for threshold: {threshold}% ---")
    df_numeric = df.select_dtypes(include=np.number)
    with pd.option_context("mode.use_inf_as_na", True):
        cols_to_keep = df_numeric.columns[df_numeric.max(skipna=True) >= threshold]
    if len(cols_to_keep) == 0:
        print(f"Warning: no columns met threshold {threshold}%. Returning empty metrics.")
        metrics_df = pd.DataFrame(index=df.index)
        metrics_df["Simpson"] = 0.0
        metrics_df["Shannon"] = 0.0
        metrics_df["Pielou"] = 0.0
        return metrics_df, 0

    print(f"Keeping {len(cols_to_keep)} columns out of {len(df_numeric.columns)} numeric columns.")
    df_filtered = df_numeric[cols_to_keep]
    metrics_df = pd.DataFrame(index=df.index)
    metrics_df["Simpson"] = df_filtered.apply(calculate_simpson_diversity, axis=1)
    metrics_df["Shannon"] = df_filtered.apply(calculate_shannon_diversity, axis=1)
    metrics_df["Pielou"] = df_filtered.apply(calculate_pielou_evenness, axis=1)
    return metrics_df, len(cols_to_keep)


def plot_diversity_vs_threshold(long_df: pd.DataFrame, category_title: str, output_path: Path) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    print(f"Generating combined mean plot for {category_title} (0-5% threshold)...")
    mean_df = long_df.groupby("Threshold")[["Simpson", "Shannon", "Pielou"]].mean().reset_index()
    plot_df = pd.melt(
        mean_df,
        id_vars=["Threshold"],
        value_vars=["Simpson", "Shannon", "Pielou"],
        var_name="MetricType",
        value_name="MeanValue",
    )
    plt.figure(figsize=(8, 6))
    ax = sns.lineplot(data=plot_df, x="Threshold", y="MeanValue", hue="MetricType", marker="o")
    ax.set_xlabel("Minimum Abundance Threshold (%)")
    ax.set_ylabel("Mean Diversity Index Value")
    ax.set_title(f"Mean {category_title} Diversity vs. Abundance Threshold (0-5%)")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, title="Metric Type")
    ax.set_xlim(left=0, right=5.0)
    ax.set_xticks(sorted(plot_df["Threshold"].unique()))
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to: {output_path}")


def build_canonical_diversity_table(breakdown_df: pd.DataFrame, comparison_df: pd.DataFrame) -> pd.DataFrame:
    required_metric_cols = ["Simpson_0.0", "Shannon_0.0", "Pielou_0.0"]
    missing_cols = [col for col in required_metric_cols if col not in comparison_df.columns]
    if missing_cols:
        raise KeyError(
            "Canonical diversity summary requires threshold-0.0 metrics; "
            f"missing columns: {missing_cols}"
        )

    canonical_df = breakdown_df.reset_index().copy()
    aligned_metrics = comparison_df.loc[breakdown_df.index, required_metric_cols].rename(
        columns={
            "Simpson_0.0": "Simpson_Diversity",
            "Shannon_0.0": "Shannon_Diversity",
            "Pielou_0.0": "Pielou_Evenness",
        }
    )
    for col in aligned_metrics.columns:
        canonical_df[col] = aligned_metrics[col].to_numpy()
    return canonical_df


def analyze_rank_stability(
    long_df: pd.DataFrame,
    category_title: str,
    output_dir: Path,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    print(f"Analyzing rank stability for {category_title}...")
    metrics = ["Simpson", "Shannon", "Pielou"]
    rank_dfs: dict[str, pd.DataFrame] = {}
    stability_results: dict[str, pd.Series] = {}
    for metric in metrics:
        long_df[f"{metric}_Rank"] = long_df.groupby("Threshold")[metric].rank(ascending=False, method="min")
        rank_pivot = long_df.pivot_table(index="Species", columns="Threshold", values=f"{metric}_Rank")
        rank_dfs[metric] = rank_pivot
        stability_results[f"{metric}_Rank_StdDev"] = rank_pivot.std(axis=1)
        stability_results[f"{metric}_Rank_Constant"] = rank_pivot.nunique(axis=1) <= 1
        rank_pivot_path = output_dir / f"{category_title.lower()}_{metric}_ranks_by_threshold.csv"
        rank_pivot.to_csv(rank_pivot_path)
        print(f"- {metric} rank pivot table saved to: {rank_pivot_path}")

    stability_df = pd.DataFrame(stability_results, index=long_df["Species"].unique())
    stability_output_path = output_dir / f"{category_title.lower()}_rank_stability_summary.csv"
    stability_df.to_csv(stability_output_path)
    print(f"- Rank stability summary saved to: {stability_output_path}")
    return rank_dfs, stability_df


def plot_rank_changes(rank_dfs: dict[str, pd.DataFrame], category_title: str, figure_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    print(f"\nGenerating rank change plot for {category_title}...")
    metrics = ["Simpson", "Shannon", "Pielou"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(8 * len(metrics), 6), sharey=True)
    fig.suptitle(f"{category_title} Species Rank vs. Abundance Threshold", fontsize=16, y=1.03)

    for i, metric in enumerate(metrics):
        ax = axes[i]
        rank_pivot = rank_dfs[metric]
        rank_long = rank_pivot.reset_index().melt(id_vars="Species", var_name="Threshold", value_name="Rank")
        rank_long["Threshold"] = pd.to_numeric(rank_long["Threshold"])
        rank_long = rank_long.sort_values(by=["Species", "Threshold"])
        sns.lineplot(data=rank_long, x="Threshold", y="Rank", hue="Species", marker="o", ax=ax, legend=i == 0)
        ax.invert_yaxis()
        ax.set_xticks(sorted(rank_long["Threshold"].unique()))
        ax.tick_params(axis="x", rotation=45)
        ax.set_xlabel("Min Abundance Threshold (%)")
        ax.set_ylabel("Rank" if i == 0 else "")
        ax.set_title(f"{metric} Diversity Rank")

    legend = axes[0].get_legend()
    if legend:
        handles = getattr(legend, "legend_handles", None) or getattr(legend, "legendHandles", [])
        labels = [text.get_text() for text in legend.texts]
        legend.remove()
        valid_handles_labels = {label: handle for handle, label in zip(handles, labels) if not label.startswith("_")}
        if valid_handles_labels:
            fig.legend(
                valid_handles_labels.values(),
                valid_handles_labels.keys(),
                title="Species",
                bbox_to_anchor=(0.92, 0.5),
                loc="center left",
                ncol=2,
                fontsize="small",
            )

    plt.tight_layout(rect=[0, 0, 0.92, 1])
    plot_path = figure_dir / f"{category_title.lower()}_rank_changes_vs_threshold.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Rank change plot saved to: {plot_path}")


def calculate_rank_correlation(
    rank_pivot_df: pd.DataFrame,
    category_title: str,
    metric_name: str,
    output_dir: Path,
) -> pd.DataFrame:
    from scipy.stats import kendalltau, spearmanr

    print(f"\nCalculating rank correlations for {category_title} - {metric_name}...")
    thresholds = rank_pivot_df.columns
    results = []
    for i in range(len(thresholds)):
        for j in range(i + 1, len(thresholds)):
            thr1 = thresholds[i]
            thr2 = thresholds[j]
            valid_ranks = rank_pivot_df[[thr1, thr2]].dropna()
            if len(valid_ranks) < 2:
                spearman_corr, spearman_p = np.nan, np.nan
                kendall_corr, kendall_p = np.nan, np.nan
            else:
                spearman_corr, spearman_p = spearmanr(valid_ranks[thr1], valid_ranks[thr2])
                kendall_corr, kendall_p = kendalltau(valid_ranks[thr1], valid_ranks[thr2])
            results.append(
                {
                    "Threshold_1": thr1,
                    "Threshold_2": thr2,
                    "Spearman_Rho": spearman_corr,
                    "Spearman_P": spearman_p,
                    "Kendall_Tau": kendall_corr,
                    "Kendall_P": kendall_p,
                    "N_Compared": len(valid_ranks),
                }
            )

    correlation_df = pd.DataFrame(results)
    corr_output_path = output_dir / f"{category_title.lower()}_{metric_name}_rank_correlations.csv"
    correlation_df.to_csv(corr_output_path, index=False)
    print(f"Rank correlations saved to: {corr_output_path}")
    return correlation_df


def run_diversity_analysis(
    order_csv: Path = ORDER_CSV,
    superfamily_csv: Path = SUPERFAMILY_CSV,
    output_data_dir: Path = OUTPUT_DATA_DIR,
    figure_dir: Path = FIGURE_DIR,
) -> None:
    output_data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading Order data from: {order_csv}")
    order_df_orig = pd.read_csv(order_csv)
    order_df_orig.set_index(order_df_orig.columns[0], inplace=True)
    print(f"Loading Superfamily data from: {superfamily_csv}")
    superfamily_df_orig = pd.read_csv(superfamily_csv)
    superfamily_df_orig.set_index(superfamily_df_orig.columns[0], inplace=True)

    order_results_list = []
    superfamily_results_list = []
    wide_order_results = []
    wide_superfamily_results = []
    column_counts_log = []

    for threshold in THRESHOLDS:
        order_metrics, order_cols_kept = calculate_diversity_for_threshold(order_df_orig, threshold)
        column_counts_log.append({"Threshold": threshold, "Category": "Order", "NumColumnsKept": order_cols_kept})
        superfamily_metrics, superfamily_cols_kept = calculate_diversity_for_threshold(superfamily_df_orig, threshold)
        column_counts_log.append(
            {"Threshold": threshold, "Category": "Superfamily", "NumColumnsKept": superfamily_cols_kept}
        )

        order_metrics_thr = order_metrics.copy()
        order_metrics_thr["Threshold"] = threshold
        order_results_list.append(order_metrics_thr.reset_index())
        superfamily_metrics_thr = superfamily_metrics.copy()
        superfamily_metrics_thr["Threshold"] = threshold
        superfamily_results_list.append(superfamily_metrics_thr.reset_index())
        wide_order_results.append(order_metrics.rename(columns=lambda col: f"{col}_{threshold}"))
        wide_superfamily_results.append(superfamily_metrics.rename(columns=lambda col: f"{col}_{threshold}"))

    col_counts_df = pd.DataFrame(column_counts_log)
    order_long_df = pd.concat(order_results_list, ignore_index=True)
    superfamily_long_df = pd.concat(superfamily_results_list, ignore_index=True)
    species_col_name = order_df_orig.index.name
    order_long_df.rename(columns={species_col_name: "Species"}, inplace=True)
    superfamily_long_df.rename(columns={species_col_name: "Species"}, inplace=True)

    plot_diversity_vs_threshold(
        order_long_df,
        "Order",
        figure_dir / "order_mean_diversity_vs_threshold_combined_0_5pct.png",
    )
    plot_diversity_vs_threshold(
        superfamily_long_df,
        "Superfamily",
        figure_dir / "superfamily_mean_diversity_vs_threshold_combined_0_5pct.png",
    )

    order_compare_df = functools.reduce(
        lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how="outer"),
        [pd.DataFrame(index=order_df_orig.index)] + wide_order_results,
    )
    superfamily_compare_df = functools.reduce(
        lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how="outer"),
        [pd.DataFrame(index=superfamily_df_orig.index)] + wide_superfamily_results,
    )

    print("\n--- Wide Comparative Order Diversity Metrics (0-5% Thresholds) ---")
    print(order_compare_df)
    print("\n--- Wide Comparative Superfamily Diversity Metrics (0-5% Thresholds) ---")
    print(superfamily_compare_df)

    order_compare_df.to_csv(output_data_dir / "comparison_diversity_order_stats_granular_0_5pct.csv")
    superfamily_compare_df.to_csv(output_data_dir / "comparison_diversity_superfamily_stats_granular_0_5pct.csv")
    order_long_df.to_csv(output_data_dir / "long_format_diversity_order_stats.csv", index=False)
    superfamily_long_df.to_csv(output_data_dir / "long_format_diversity_superfamily_stats.csv", index=False)

    canonical_order_df = build_canonical_diversity_table(order_df_orig, order_compare_df)
    canonical_superfamily_df = build_canonical_diversity_table(superfamily_df_orig, superfamily_compare_df)
    canonical_order_df.to_csv(CANONICAL_ORDER_CSV, index=False)
    canonical_superfamily_df.to_csv(CANONICAL_SUPERFAMILY_CSV, index=False)
    print("\nCanonical diversity summary tables saved:")
    print(f"- Order: {CANONICAL_ORDER_CSV}")
    print(f"- Superfamily: {CANONICAL_SUPERFAMILY_CSV}")

    col_counts_df.to_csv(output_data_dir / "threshold_column_counts.csv", index=False)
    print(f"\nColumn counts saved to: {output_data_dir / 'threshold_column_counts.csv'}")

    order_rank_pivots, order_stability_df = analyze_rank_stability(order_long_df, "Order", output_data_dir)
    print("\nOrder Rank Stability Summary:")
    print(order_stability_df)
    superfamily_rank_pivots, superfamily_stability_df = analyze_rank_stability(
        superfamily_long_df,
        "Superfamily",
        output_data_dir,
    )
    print("\nSuperfamily Rank Stability Summary:")
    print(superfamily_stability_df)

    plot_rank_changes(order_rank_pivots, "Order", figure_dir)
    plot_rank_changes(superfamily_rank_pivots, "Superfamily", figure_dir)

    for metric, pivot_df in order_rank_pivots.items():
        calculate_rank_correlation(pivot_df, "Order", metric, output_data_dir)
    for metric, pivot_df in superfamily_rank_pivots.items():
        calculate_rank_correlation(pivot_df, "Superfamily", metric, output_data_dir)
    print("\nRank correlation analysis complete.")


def main() -> None:
    run_diversity_analysis()


if __name__ == "__main__":
    main()
