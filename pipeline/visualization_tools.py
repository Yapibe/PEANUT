import logging
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

logger = logging.getLogger(__name__)

# Suppress debug messages from multipart and matplotlib.font_manager
logging.getLogger("multipart").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)


def passes_size_filter(pathway, pathway_sizes, settings):
    """
    Check if a pathway passes the size filter based on its gene count.

    Parameters:
    - pathway (str): Pathway name.
    - pathway_sizes (dict): Dictionary mapping pathways to their sizes.
    - settings: Settings object containing size thresholds.

    Returns:
    - bool: True if the pathway passes the filter, False otherwise.
    """
    pathway_size = pathway_sizes.get(pathway, 0)
    return settings.minimum_gene_per_pathway <= pathway_size <= settings.maximum_gene_per_pathway


def calculate_mean_score(pathway, genes_by_pathway, scores):
    """
    Calculate the mean score for a pathway.

    Parameters:
    - pathway (str): Pathway name.
    - genes_by_pathway (dict): Dictionary mapping pathways to their gene sets.
    - scores (dict): Dictionary mapping genes to their scores.

    Returns:
    - float: Mean score for the pathway.
    """
    pathway_genes = genes_by_pathway.get(pathway, set())
    pathway_gene_scores = [scores[gene_id] for gene_id in pathway_genes if gene_id in scores]
    return np.mean(pathway_gene_scores) if pathway_gene_scores else 0



# Main function: Plot pathways
def plot_pathways_mean_scores(settings, all_pathways, gsea_pathways=None, pathway_sizes=None, genes_by_pathway=None, scores=None) -> Path:
    """
    Plot mean scores of pathways across all conditions, ensuring that only relevant pathways are included.

    Parameters:
    - settings: Settings object containing plot settings and thresholds.
    - all_pathways (dict): Significant pathway results by condition.
    - gsea_pathways (dict): GSEA pathway results by condition (optional).
    - pathway_sizes (dict): Dictionary containing pathway sizes.
    - genes_by_pathway (dict): Dictionary mapping pathways to their gene sets.
    - scores (dict): Dictionary mapping genes to their scores.

    Returns:
    - Path: Path to the saved plot image.
    """
    if not all_pathways and not gsea_pathways:
        logger.info("No pathways to plot. Exiting.")
        return None

    # Combine pathway data
    combined_pathways = {**all_pathways, **(gsea_pathways or {})}

    # Identify relevant pathways
    relevant_pathways = {
        pathway
        for pathway, conditions in combined_pathways.items()
        if any(cond.get("P-value", 1) < settings.FDR_THRESHOLD for cond in conditions.values())
    }

    # Filter pathways based on relevance
    filtered_pathways = {
        pathway: conditions
        for pathway, conditions in combined_pathways.items()
        if pathway in relevant_pathways
    }

    # Identify all conditions
    all_conditions = sorted({cond for pathways in filtered_pathways.values() for cond in pathways})

    # Prepare data for plotting
    data_records = []
    for pathway, conditions in filtered_pathways.items():
        for condition in all_conditions:
            if condition in conditions:
                data_records.append({
                    "Pathway": pathway,
                    "Condition": condition,
                    "Mean_Score": conditions[condition].get("Mean", 0),
                    "P-value": conditions[condition].get("P-value", 1),
                    "Significant": True,
                })

            else:
                if passes_size_filter(pathway, pathway_sizes, settings):
                    mean_score = calculate_mean_score(pathway, genes_by_pathway, scores)
                    data_records.append({
                        "Pathway": pathway,
                        "Condition": condition,
                        "Mean_Score": mean_score,
                        "P-value": 1,
                        "Significant": False,
                    })

    # Create DataFrame and sort pathways
    data_df = pd.DataFrame(data_records)
    mean_ranks = data_df.pivot(index="Pathway", columns="Condition", values="Mean_Score").rank(
        ascending=False, method="average", axis=0
    ).mean(axis=1)
    sorted_pathways = mean_ranks.sort_values().index.tolist()
    data_df["Pathway"] = pd.Categorical(data_df["Pathway"], categories=sorted_pathways, ordered=True)
    data_df.sort_values("Pathway", inplace=True)

    # Plot data
    plt.figure(figsize=(16, max(6, len(sorted_pathways) * 0.3)))
    ax = plt.subplot(111)
    bar_width = 0.8 / len(all_conditions)
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(all_conditions)))

    for i, condition in enumerate(all_conditions):
        condition_data = data_df[data_df["Condition"] == condition]
        condition_data = condition_data.set_index("Pathway").reindex(sorted_pathways).reset_index()

        mean_scores = condition_data["Mean_Score"].values
        significance = condition_data["Significant"].values

        bar_styles = [
            {"color": colors[i], "edgecolor": "black", "hatch": "//" if sig else None}
            for sig in significance
        ]
        ax.barh(
            np.arange(len(sorted_pathways)) + bar_width * i,
            mean_scores,
            height=bar_width,
            align="center",
            color=[style["color"] for style in bar_styles],
            edgecolor=[style["edgecolor"] for style in bar_styles],
            hatch=[style["hatch"] for style in bar_styles],
            label=condition,
        )

    ax.set_yticks(np.arange(len(sorted_pathways)) + bar_width * (len(all_conditions) - 1) / 2)
    ax.set_yticklabels([pathway.replace("_", " ") for pathway in sorted_pathways], fontsize=10)
    ax.set_ylabel("Pathways", fontsize=12)
    ax.set_xlabel("Mean Scores", fontsize=12)
    ax.set_title(settings.figure_title, fontsize=14)

    # Add legend
    condition_patches = [Patch(color=colors[i], label=f"{condition}") for i, condition in enumerate(all_conditions)]
    significance_patch = Patch(facecolor="white", edgecolor="black", hatch="//", label="Significant Pathways")
    ax.legend(handles=condition_patches + [significance_patch], prop={"size": 10}, title="Legend",
              bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    output_dir = settings.plot_output_path
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file_path = output_dir / f"{settings.experiment_name}_significant_pathways_plot.png"
    plt.savefig(output_file_path, format="png", bbox_inches="tight")
    plt.close()

    return output_file_path








