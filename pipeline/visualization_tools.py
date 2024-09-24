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


def print_aggregated_pathway_information(
    settings: Settings, all_pathways: dict
) -> Path:
    """
    Print aggregated pathway information.

    Including P-values, trends, and significant genes
    for each pathway to a text file based on a given experiment.

    Parameters:
    - settings (Settings): General arguments and settings.
    - all_pathways (dict): Dictionary containing pathway information across different conditions.

    Returns:
    - file_path (Path): Path to the output file.
    """
    # Define the path for the output file
    output_dir = settings.root_folder / "Outputs" / "Text"
    output_dir.mkdir(parents=True, exist_ok=True)
    file_path = (
        output_dir
        / f"{settings.experiment_name}_{settings.pathway_file}_{settings.alpha}_aggregated.txt"
    )

    if not all_pathways:
        with file_path.open("w") as file:
            file.write("No significant pathways found.\n")
        logger.info(
            "No significant pathways found. File written with message."
        )
        return file_path

    # Create a DataFrame from all_pathways
    pathway_records = []
    for pathway, conditions in all_pathways.items():
        best_p_value = min(
            condition_data["P-value"]
            for condition_data in conditions.values()
        )
        trends = ", ".join(
            f"{condition_name}: {conditions[condition_name]['Trend']}"
            for condition_name in conditions
        )
        significant_genes = {}
        for _condition_name, condition_data in conditions.items():
            for gene_id, gene_info in condition_data.get(
                "significant_genes", {}
            ).items():
                if gene_id not in significant_genes:
                    significant_genes[gene_id] = {
                        "Symbol": gene_info["Symbol"],
                        "Scores": [],
                    }
                significant_genes[gene_id]["Scores"].append(
                    gene_info["Score"]
                )
        pathway_records.append(
            {
                "Pathway": pathway,
                "Best_P-value": best_p_value,
                "Trends": trends,
                "Significant_Genes": significant_genes,
            }
        )

    # Sort pathways by the best (lowest) P-value
    pathway_df = pd.DataFrame(pathway_records)
    pathway_df.sort_values(by="Best_P-value", inplace=True)

    # Write to the output file
    with file_path.open("w") as file:
        for _, row in pathway_df.iterrows():
            file.write(
                f"Pathway: {row['Pathway']} P-value: {row['Best_P-value']:.5f}\n"
            )
            file.write(f"Trends: {row['Trends']}\n")
            file.write("Significant Genes:\n")
            for _gene_id, gene_data in row[
                "Significant_Genes"
            ].items():
                scores_str = ", ".join(map(str, gene_data["Scores"]))
                file.write(
                    f"    {gene_data['Symbol']}: {scores_str}\n"
                )
            file.write("\n")
    logger.info(
        f"Aggregated pathway information written to {file_path}"
    )
    return file_path


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

    # Initialize DataFrames to store mean scores and p-values for each condition
    data_records = []
    for pathway, conditions in all_pathways.items():
        for condition_name, condition_data in conditions.items():
            data_records.append(
                {
                    "Pathway": pathway,
                    "Condition": condition_name,
                    "Mean_Score": condition_data.get("Mean", 0),
                    "P-value": condition_data.get("P-value", 1),
                }
            )
    data_df = pd.DataFrame(data_records)

    if data_df.empty:
        logger.info("Data for plotting is empty. Exiting function.")
        return None

    # Filter pathways with p-value <= 0.05 for at least one condition
    significant_pathways = data_df[data_df["P-value"] <= 0.05][
        "Pathway"
    ].unique()
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
    output_dir = settings.root_folder / "Outputs" / "Plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file_path = (
        output_dir / f"{settings.experiment_name}_plot.png"
    )
    plt.savefig(output_file_path, format="png", bbox_inches="tight", bbox_extra_artists=[legend])
    plt.close()
    return output_file_path
