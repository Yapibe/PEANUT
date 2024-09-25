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
    settings: Settings, all_pathways: dict
) -> Path:
    """
    Plot mean scores of pathways across all conditions and save the plot as a PNG file.

    Parameters:
    - settings (Settings): General arguments and settings.
    - all_pathways (dict): Dictionary containing pathway information across different conditions.

    Returns:
    - output_file_path (Path): Path to the output plot file.
    """
    if not all_pathways:
        logger.info("No pathways to plot. Exiting function.")
        return None

    # Initialize DataFrames with default values from settings
    data_records = []
    for pathway, conditions in all_pathways.items():
        for condition_name, condition_data in conditions.items():
            data_records.append({
                "Pathway": pathway,
                "Condition": condition_name,
                "Mean_Score": condition_data.get("Mean", 0),
                "P-value": condition_data.get("P-value", 1),
            })
    data_df = pd.DataFrame(data_records)

    if data_df.empty:
        logger.info("Data for plotting is empty. Exiting function.")
        return None

    # Filter pathways with p-value <= p_value_threshold for at least one condition
    significant_pathways = data_df[data_df["P-value"] <= settings.FDR_THRESHOLD]["Pathway"].unique()
    data_df = data_df[data_df["Pathway"].isin(significant_pathways)]

    if data_df.empty:
        logger.info(
            "No significant pathways to plot. Exiting function."
        )
        return None

    # Pivot the DataFrame to have pathways as index and conditions as columns
    mean_scores_df = data_df.pivot(
        index="Pathway", columns="Condition", values="Mean_Score"
    ).fillna(0)
    p_values_df = data_df.pivot(
        index="Pathway", columns="Condition", values="P-value"
    ).fillna(1)

    # Sort pathways by mean score
    mean_scores_df = mean_scores_df.loc[
        mean_scores_df.mean(axis=1).sort_values().index
    ]
    p_values_df = p_values_df.reindex(mean_scores_df.index)

    # Prepare data for plotting
    conditions = mean_scores_df.columns.tolist()
    pathways = mean_scores_df.index.tolist()
    num_conditions = len(conditions)
    num_pathways = len(pathways)
    positions = np.arange(num_pathways)
    bar_width = 0.8 / num_conditions

    # Create a figure
    plt.figure(figsize=(16, max(6, num_pathways * 0.3)))
    ax = plt.subplot(111)

    # Generate a color map for the conditions
    colors = plt.get_cmap("viridis")(
        np.linspace(0, 1, num_conditions)
    )

    # Plot each condition's mean scores for each pathway
    for i, condition in enumerate(conditions):
        mean_scores = mean_scores_df[condition].values
        p_values = p_values_df[condition].values
        # Determine bar styles based on p-value significance
        bar_styles = [
            {
                "color": colors[i],
                "edgecolor": "black",
                "hatch": (
                    "//" if p_val > settings.FDR_THRESHOLD else None
                ),
            }
            for p_val in p_values
        ]
        ax.barh(
            positions + bar_width * i,
            mean_scores,
            height=bar_width,
            align="center",
            color=[style["color"] for style in bar_styles],
            edgecolor=[style["edgecolor"] for style in bar_styles],
            hatch=[style["hatch"] for style in bar_styles],
        )

    # Set y-axis labels to be pathway names
    ax.set_yticks(positions + bar_width * (num_conditions - 1) / 2)
    ax.set_yticklabels(
        [pathway.replace("_", " ") for pathway in pathways],
        fontsize=10,
    )

    # Label axes and set title
    ax.set_ylabel("Pathways", fontsize=12)
    ax.set_xlabel("Mean Scores", fontsize=12)
    ax.set_title(
        "Pathway Mean Scores Across Different Conditions", fontsize=14
    )

    # Create a legend for the conditions
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=colors[i])
        for i in range(num_conditions)
    ]
    legend_handles.append(plt.Rectangle((0, 0), 1, 1, facecolor='none', edgecolor='black', hatch='//'))
    legend_labels = conditions + ["Non-significant"]
    legend = ax.legend(
        legend_handles,
        legend_labels,
        prop={"size": 10},
        title="Conditions",
        bbox_to_anchor=(1.05, 1), loc='upper left'
    )

    # Adjust subplot layout to avoid clipping of tick-labels
    plt.tight_layout()

    # Save the figure to a PNG file
    output_dir = settings.plot_output_path
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file_path = (
        output_dir / f"{settings.experiment_name}_plot.png"
    )
    plt.savefig(output_file_path, format="png", bbox_inches="tight", bbox_extra_artists=[legend])
    plt.close()
    return output_file_path
