"""
Diversity Statistics Script

Calculates diversity indices (Shannon, Simpson, Pielou's evenness) for
transposable element compositions across Desmognathus species.
"""

import pandas as pd
import numpy as np
import os
import sys
import functools
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, kendalltau
from pathlib import Path

# Import centralized configuration
sys.path.insert(0, str(Path(__file__).parent.parent))  # Add scripts/ to path
from config import paths, PROJECT_ROOT

# --- Configuration ---
INPUT_DATA_DIR = paths.results.data
OUTPUT_DATA_DIR = paths.results.data
FIGURE_DIR = paths.results.figures / "diversity_threshold_plots"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

ORDER_CSV = INPUT_DATA_DIR / "dnaPipeTE_order_breakdown.csv"
SUPERFAMILY_CSV = INPUT_DATA_DIR / "dnaPipeTE_superfamily_breakdown.csv"

def calculate_simpson_diversity(row):
    row = row[row > 0]
    N = row.sum()
    if N <= 1 or pd.isna(N) or N * (N - 1) == 0:
        return 0.0
    simpson_dominance = np.sum(row * (row - 1)) / (N * (N - 1))
    return 1.0 - simpson_dominance

def calculate_shannon_diversity(row):
    row = row[row > 0]
    if row.empty or pd.isna(row.sum()) or row.sum() == 0:
        return 0.0
    proportions = row / row.sum()
    return -np.sum(proportions * np.log(proportions))

def calculate_pielou_evenness(row):
    row = row[row > 0]
    num_categories = len(row)
    if num_categories <= 1:
        return 0.0
    shannon_index = calculate_shannon_diversity(row)
    if pd.isna(shannon_index):
        return 0.0
    return shannon_index / np.log(num_categories)

def calculate_diversity_for_threshold(df, threshold):
    print(f"--- Calculating diversity for threshold: {threshold}% ---")
    df_numeric = df.select_dtypes(include=np.number)
    with pd.option_context('mode.use_inf_as_na', True):
        cols_to_keep = df_numeric.columns[df_numeric.max(skipna=True) >= threshold]
    if len(cols_to_keep) == 0:
        print(f"Warning: No columns met threshold {threshold}%. Returning empty metrics.")
        metrics_df = pd.DataFrame(index=df.index)
        metrics_df['Simpson'] = 0.0
        metrics_df['Shannon'] = 0.0
        metrics_df['Pielou'] = 0.0
        return metrics_df
    print(f"Keeping {len(cols_to_keep)} columns out of {len(df_numeric.columns)} numeric columns.")
    df_filtered = df_numeric[cols_to_keep]
    metrics_df = pd.DataFrame(index=df.index)
    metrics_df['Simpson'] = df_filtered.apply(calculate_simpson_diversity, axis=1)
    metrics_df['Shannon'] = df_filtered.apply(calculate_shannon_diversity, axis=1)
    metrics_df['Pielou'] = df_filtered.apply(calculate_pielou_evenness, axis=1)
    num_cols_kept = len(cols_to_keep)
    return metrics_df, num_cols_kept

def plot_diversity_vs_threshold(long_df, category_title, output_path):
    print(f'Generating combined mean plot for {category_title} (0-5% Threshold)...')
    mean_df = long_df.groupby('Threshold')[['Simpson', 'Shannon', 'Pielou']].mean().reset_index()
    plot_df = pd.melt(mean_df, id_vars=['Threshold'], value_vars=['Simpson', 'Shannon', 'Pielou'], var_name='MetricType', value_name='MeanValue')
    plt.figure(figsize=(8, 6))
    ax = sns.lineplot(data=plot_df, x='Threshold', y='MeanValue', hue='MetricType', marker='o')
    ax.set_xlabel('Minimum Abundance Threshold (%)')
    ax.set_ylabel('Mean Diversity Index Value')
    ax.set_title(f'Mean {category_title} Diversity vs. Abundance Threshold (0-5%)')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, title='Metric Type')
    ax.set_xlim(left=0, right=5.0)
    x_ticks = sorted(plot_df['Threshold'].unique())
    ax.set_xticks(x_ticks)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Plot saved to: {output_path}')

# Load data
print(f"Loading Order data from: {ORDER_CSV}")
order_df_orig = pd.read_csv(ORDER_CSV)
order_df_orig.set_index(order_df_orig.columns[0], inplace=True)
print(f"Loading Superfamily data from: {SUPERFAMILY_CSV}")
superfamily_df_orig = pd.read_csv(SUPERFAMILY_CSV)
superfamily_df_orig.set_index(superfamily_df_orig.columns[0], inplace=True)

thresholds = [5.0, 1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0]
order_results_list = []
superfamily_results_list = []
wide_order_results = []
wide_superfamily_results = []
column_counts_log = []

for thr in thresholds:
    order_metrics, order_cols_kept = calculate_diversity_for_threshold(order_df_orig, thr)
    column_counts_log.append({'Threshold': thr, 'Category': 'Order', 'NumColumnsKept': order_cols_kept})
    superfamily_metrics, superfamily_cols_kept = calculate_diversity_for_threshold(superfamily_df_orig, thr)
    column_counts_log.append({'Threshold': thr, 'Category': 'Superfamily', 'NumColumnsKept': superfamily_cols_kept})
    order_metrics_thr = order_metrics.copy()
    order_metrics_thr['Threshold'] = thr
    order_results_list.append(order_metrics_thr.reset_index())
    superfamily_metrics_thr = superfamily_metrics.copy()
    superfamily_metrics_thr['Threshold'] = thr
    superfamily_results_list.append(superfamily_metrics_thr.reset_index())
    wide_order_results.append(order_metrics.rename(columns=lambda c: f'{c}_{thr}'))
    wide_superfamily_results.append(superfamily_metrics.rename(columns=lambda c: f'{c}_{thr}'))

col_counts_df = pd.DataFrame(column_counts_log)
order_long_df = pd.concat(order_results_list, ignore_index=True)
superfamily_long_df = pd.concat(superfamily_results_list, ignore_index=True)
species_col_name = order_df_orig.index.name
order_long_df.rename(columns={species_col_name: 'Species'}, inplace=True)
superfamily_long_df.rename(columns={species_col_name: 'Species'}, inplace=True)

order_plot_path = os.path.join(FIGURE_DIR, 'order_mean_diversity_vs_threshold_combined_0_5pct.png')
superfamily_plot_path = os.path.join(FIGURE_DIR, 'superfamily_mean_diversity_vs_threshold_combined_0_5pct.png')
plot_diversity_vs_threshold(order_long_df, 'Order', order_plot_path)
plot_diversity_vs_threshold(superfamily_long_df, 'Superfamily', superfamily_plot_path)

order_compare_df = pd.DataFrame(index=order_df_orig.index)
superfamily_compare_df = pd.DataFrame(index=superfamily_df_orig.index)
order_compare_df = functools.reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how='outer'), [order_compare_df] + wide_order_results)
superfamily_compare_df = functools.reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how='outer'), [superfamily_compare_df] + wide_superfamily_results)

print('\n--- Wide Comparative Order Diversity Metrics (0-5% Thresholds) ---')
print(order_compare_df)
print('\n--- Wide Comparative Superfamily Diversity Metrics (0-5% Thresholds) ---')
print(superfamily_compare_df)

order_output_path = os.path.join(OUTPUT_DATA_DIR, 'comparison_diversity_order_stats_granular_0_5pct.csv')
superfamily_output_path = os.path.join(OUTPUT_DATA_DIR, 'comparison_diversity_superfamily_stats_granular_0_5pct.csv')
order_compare_df.to_csv(order_output_path)
superfamily_compare_df.to_csv(superfamily_output_path)
print(f'\nWide comparative diversity statistics saved:')
print(f'- Order: {order_output_path}')
print(f'- Superfamily: {superfamily_output_path}')

# Save long format data
order_long_output_path = os.path.join(OUTPUT_DATA_DIR, 'long_format_diversity_order_stats.csv')
superfamily_long_output_path = os.path.join(OUTPUT_DATA_DIR, 'long_format_diversity_superfamily_stats.csv')
order_long_df.to_csv(order_long_output_path, index=False)
superfamily_long_df.to_csv(superfamily_long_output_path, index=False)
print(f'\nLong format diversity statistics saved:')
print(f'- Order: {order_long_output_path}')
print(f'- Superfamily: {superfamily_long_output_path}')

print('\n--- Number of Columns Kept per Threshold ---')
print(col_counts_df.to_string())
col_counts_output_path = os.path.join(OUTPUT_DATA_DIR, 'threshold_column_counts.csv')
col_counts_df.to_csv(col_counts_output_path, index=False)
print(f'\nColumn counts saved to: {col_counts_output_path}')

print('\n--- Analyzing Rank Stability Across Thresholds ---')
def analyze_rank_stability(long_df, category_title, output_dir):
    print(f"Analyzing rank stability for {category_title}...")
    metrics = ['Simpson', 'Shannon', 'Pielou']
    rank_dfs = {}
    stability_results = {}
    for metric in metrics:
        long_df[f'{metric}_Rank'] = long_df.groupby('Threshold')[metric].rank(ascending=False, method='min')
        rank_pivot = long_df.pivot_table(index='Species', columns='Threshold', values=f'{metric}_Rank')
        rank_dfs[metric] = rank_pivot
        rank_std = rank_pivot.std(axis=1)
        stability_results[f'{metric}_Rank_StdDev'] = rank_std
        is_rank_constant = rank_pivot.nunique(axis=1) <= 1
        stability_results[f'{metric}_Rank_Constant'] = is_rank_constant
        rank_pivot_path = os.path.join(output_dir, f'{category_title.lower()}_{metric}_ranks_by_threshold.csv')
        rank_pivot.to_csv(rank_pivot_path)
        print(f"- {metric} rank pivot table saved to: {rank_pivot_path}")
    stability_df = pd.DataFrame(stability_results, index=long_df['Species'].unique())
    stability_output_path = os.path.join(output_dir, f'{category_title.lower()}_rank_stability_summary.csv')
    stability_df.to_csv(stability_output_path)
    print(f"- Rank stability summary saved to: {stability_output_path}")
    print(f"Rank stability analysis for {category_title} complete.")
    return rank_dfs, stability_df

order_rank_pivots, order_stability_df = analyze_rank_stability(order_long_df, 'Order', OUTPUT_DATA_DIR)
print("\nOrder Rank Stability Summary:")
print(order_stability_df)
superfamily_rank_pivots, superfamily_stability_df = analyze_rank_stability(superfamily_long_df, 'Superfamily', OUTPUT_DATA_DIR)
print("\nSuperfamily Rank Stability Summary:")
print(superfamily_stability_df)

def plot_rank_changes(rank_dfs, category_title, figure_dir):
    """Plots species rank changes across thresholds for each metric in a faceted plot."""
    print(f"\nGenerating rank change plot for {category_title}...")
    metrics = ['Simpson', 'Shannon', 'Pielou']
    num_metrics = len(metrics)
    fig, axes = plt.subplots(1, num_metrics, figsize=(8 * num_metrics, 6), sharey=True)
    fig.suptitle(f'{category_title} Species Rank vs. Abundance Threshold', fontsize=16, y=1.03)

    for i, metric in enumerate(metrics):
        ax = axes[i]
        rank_pivot = rank_dfs[metric]
        rank_long = rank_pivot.reset_index().melt(id_vars='Species', var_name='Threshold', value_name='Rank')
        rank_long['Threshold'] = pd.to_numeric(rank_long['Threshold'])
        rank_long = rank_long.sort_values(by=['Species', 'Threshold'])

        # Let seaborn create legend only on the first plot
        current_legend_setting = True if i == 0 else False
        sns.lineplot(data=rank_long, x='Threshold', y='Rank', hue='Species', marker='o', ax=ax, legend=current_legend_setting)

        ax.invert_yaxis()
        x_ticks = sorted(rank_long['Threshold'].unique())
        ax.set_xticks(x_ticks)
        ax.tick_params(axis='x', rotation=45)
        ax.set_xlabel('Min Abundance Threshold (%)')
        ax.set_ylabel('Rank' if i == 0 else '')
        ax.set_title(f'{metric} Diversity Rank')

    # Get legend from the first axis, remove it, and place it on the figure
    legend = axes[0].get_legend()
    if legend:
        handles = legend.legend_handles
        labels = [text.get_text() for text in legend.texts]
        legend.remove()
        if not handles or not labels:
            print("Warning: Legend found on axis 0, but no handles or labels extracted.")
        else:
             # Exclude internal labels (often start with _) if necessary, though usually handled by seaborn
            valid_handles_labels = {l: h for h, l in zip(handles, labels) if not l.startswith('_')}
            # Adjust bbox_to_anchor x-coordinate to move legend closer
            fig.legend(valid_handles_labels.values(), valid_handles_labels.keys(), title='Species', bbox_to_anchor=(.92, 0.5), loc='center left', ncol=2, fontsize='small')
    else:
        print("Warning: No legend found on the first axis to move to figure level.")

    # Adjust tight_layout rect to reduce space on the right
    plt.tight_layout(rect=[0, 0, 0.92, 1]) 
    plot_path = os.path.join(figure_dir, f'{category_title.lower()}_rank_changes_vs_threshold.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Rank change plot saved to: {plot_path}")

# Plot rank changes
plot_rank_changes(order_rank_pivots, 'Order', FIGURE_DIR)
plot_rank_changes(superfamily_rank_pivots, 'Superfamily', FIGURE_DIR)

def calculate_rank_correlation(rank_pivot_df, category_title, metric_name, output_dir):
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
            results.append({
                'Threshold_1': thr1,
                'Threshold_2': thr2,
                'Spearman_Rho': spearman_corr,
                'Spearman_P': spearman_p,
                'Kendall_Tau': kendall_corr,
                'Kendall_P': kendall_p,
                'N_Compared': len(valid_ranks)
            })
    correlation_df = pd.DataFrame(results)
    corr_output_path = os.path.join(output_dir, f'{category_title.lower()}_{metric_name}_rank_correlations.csv')
    correlation_df.to_csv(corr_output_path, index=False)
    print(f"Rank correlations saved to: {corr_output_path}")
    return correlation_df

all_correlations = {}
for metric, pivot_df in order_rank_pivots.items():
    all_correlations[f'Order_{metric}'] = calculate_rank_correlation(pivot_df, 'Order', metric, OUTPUT_DATA_DIR)
for metric, pivot_df in superfamily_rank_pivots.items():
    all_correlations[f'Superfamily_{metric}'] = calculate_rank_correlation(pivot_df, 'Superfamily', metric, OUTPUT_DATA_DIR)
print("\nRank correlation analysis complete.")