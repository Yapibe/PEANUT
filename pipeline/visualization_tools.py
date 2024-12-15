import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .settings import Settings

logger = logging.getLogger(__name__)

# Suppress debug messages from multipart and matplotlib.font_manager
logging.getLogger("multipart").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)


def plot_pathways_mean_scores(
    settings: Settings, all_pathways: dict, gsea_pathways: dict = None
) -> Path:
    """
    Plot mean scores of pathways across all conditions and save the plot as a PNG file.

    Parameters:
    - settings (Settings): General arguments and settings.
    - all_pathways (dict): Dictionary containing propagation pathway information across different conditions.
    - gsea_pathways (dict): Dictionary containing GSEA pathway information across different conditions (optional).

    Returns:
    - output_file_path (Path): Path to the output plot file.
    """
    if not all_pathways and not gsea_pathways:
        logger.info("No pathways to plot for either propagation or GSEA. Exiting function.")
        return None

    # Combine all_pathways and gsea_pathways into a single dictionary
    combined_pathways = {}
    if all_pathways:
        combined_pathways.update(all_pathways)
    if gsea_pathways:
        for pathway, gsea_conditions in gsea_pathways.items():
            if pathway not in combined_pathways:
                combined_pathways[pathway] = {}
            combined_pathways[pathway].update(gsea_conditions)

    # Initialize DataFrame with pathway data
    data_records = []
    for pathway, conditions in combined_pathways.items():
        for condition_name, condition_data in conditions.items():
            data_records.append({
                "Pathway": pathway,
                "Condition": condition_name,
                "Mean_Score": condition_data.get("Mean", 0),
                "P-value": condition_data.get("P-value", 1),
                "Significant": condition_data.get("P-value", 1) <= settings.FDR_THRESHOLD,
            })
    data_df = pd.DataFrame(data_records)

    if data_df.empty:
        logger.info("Data for plotting is empty after combining pathways. Exiting function.")
        return None

    # Filter pathways where at least one condition is significant
    significant_pathways = data_df[data_df["Significant"]]["Pathway"].unique()
    data_df = data_df[data_df["Pathway"].isin(significant_pathways)]

    if data_df.empty:
        logger.info("No significant pathways to plot. Exiting function.")
        return None

    # Sort pathways based on average rank across all conditions
    ranks_df = data_df.pivot(index="Pathway", columns="Condition", values="Mean_Score").rank(ascending=False, method="average", axis=0)
    mean_ranks = ranks_df.mean(axis=1)
    sorted_pathways = mean_ranks.sort_values().index.tolist()

    # Reorder DataFrame based on sorted pathways
    data_df["Pathway"] = pd.Categorical(data_df["Pathway"], categories=sorted_pathways, ordered=True)
    data_df.sort_values("Pathway", inplace=True)

    # Prepare data for plotting
    pathways = sorted_pathways
    conditions = data_df["Condition"].unique()
    num_conditions = len(conditions)
    num_pathways = len(pathways)
    positions = np.arange(num_pathways)
    bar_width = 0.8 / num_conditions

    # Create figure
    plt.figure(figsize=(16, max(6, num_pathways * 0.3)))
    ax = plt.subplot(111)

    # Generate a distinct color map for the conditions
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, num_conditions))

    # Plot bars for each condition
    for i, condition in enumerate(conditions):
        condition_data = data_df[data_df["Condition"] == condition]
        mean_scores = condition_data["Mean_Score"].values
        significance = condition_data["Significant"].values

        # Determine bar styles based on significance
        bar_styles = [
            {
                "color": colors[i],
                "edgecolor": "black",
                "hatch": "//" if not sig else None,
            }
            for sig in significance
        ]
        ax.barh(
            positions + bar_width * i,
            mean_scores,
            height=bar_width,
            align="center",
            color=[style["color"] for style in bar_styles],
            edgecolor=[style["edgecolor"] for style in bar_styles],
            hatch=[style["hatch"] for style in bar_styles],
            label=condition,
        )

    # Set y-axis labels
    ax.set_yticks(positions + bar_width * (num_conditions - 1) / 2)
    ax.set_yticklabels([pathway.replace("_", " ") for pathway in pathways], fontsize=10)

    # Label axes and set title
    ax.set_ylabel("Pathways", fontsize=12)
    ax.set_xlabel("Mean Scores", fontsize=12)
    ax.set_title("Pathway Mean Scores Across Different Conditions", fontsize=14)

    # Create legend
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles,
        labels,
        prop={"size": 10},
        title="Conditions",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )

    # Adjust layout and save plot
    plt.tight_layout()
    output_dir = settings.plot_output_path
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file_path = output_dir / f"{settings.experiment_name}_plot.png"
    plt.savefig(output_file_path, format="png", bbox_inches="tight", bbox_extra_artists=[legend])
    plt.close()

    return output_file_path





