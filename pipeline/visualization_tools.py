"""
Visualization tools for the PEANUT pipeline.

This module contains functions to visualize the results of pathway enrichment
analysis, including generating plots for pathway scores across conditions.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Set
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from .settings import Settings

logger = logging.getLogger(__name__)

# Suppress debug messages from multipart and matplotlib.font_manager
logging.getLogger("multipart").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)


def passes_size_filter(pathway: str, pathway_sizes: Dict[str, int], settings: Settings) -> bool:
    """
    Check if a pathway passes the size filter based on its gene count.

    Parameters:
    - pathway (str): Pathway name.
    - pathway_sizes (dict): Dictionary mapping pathways to their sizes.
    - settings (Settings): Settings object containing size thresholds.

    Returns:
    - bool: True if the pathway passes the filter, False otherwise.
    """
    pathway_size = pathway_sizes.get(pathway, 0)
    return settings.minimum_gene_per_pathway <= pathway_size <= settings.maximum_gene_per_pathway


def calculate_mean_score(pathway: str, genes_by_pathway: Dict[str, Set[str]], scores: Dict[str, float]) -> float:
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


def plot_pathways_mean_scores(
    settings: Settings, 
    all_pathways: Dict[str, Dict[str, Dict[str, float]]], 
    gsea_pathways: Optional[Dict[str, Dict[str, Dict[str, float]]]] = None, 
    pathway_sizes: Optional[Dict[str, int]] = None, 
    genes_by_pathway: Optional[Dict[str, Set[str]]] = None, 
    scores: Optional[Dict[str, float]] = None
) -> Optional[Path]:
    """
    Plot mean scores of pathways across all conditions, ensuring that only relevant pathways are included.

    Parameters:
    - settings (Settings): Settings object containing plot settings and thresholds.
    - all_pathways (Dict[str, Dict[str, Dict[str, float]]]): Significant pathway results by condition.
    - gsea_pathways (Dict[str, Dict[str, Dict[str, float]]]): GSEA pathway results by condition (optional).
    - pathway_sizes (Dict[str, int]): Dictionary containing pathway sizes.
    - genes_by_pathway (Dict[str, Set[str]]): Dictionary mapping pathways to their gene sets.
    - scores (Dict[str, float]): Dictionary mapping genes to their scores.

    Returns:
    - Optional[Path]: Path to the saved plot image, or None if no plot was created.
    """
    try:
        if not all_pathways and not gsea_pathways:
            logger.info("No pathways to plot. Exiting.")
            return None
            
        # Initialize default values for optional parameters
        if pathway_sizes is None:
            pathway_sizes = {}
        if genes_by_pathway is None:
            genes_by_pathway = {}
        if scores is None:
            scores = {}
            
        # Validate inputs
        if not isinstance(all_pathways, dict):
            logger.error(f"all_pathways must be a dict, got {type(all_pathways)}")
            return None
            
        logger.info("Starting pathway plotting...")

        # Combine pathway data
        combined_pathways = {**all_pathways, **(gsea_pathways or {})}

        # Identify relevant pathways
        relevant_pathways = {
            pathway
            for pathway, conditions in combined_pathways.items()
            if any(cond.get("P-value", 1) < settings.FDR_THRESHOLD for cond in conditions.values())
        }
        
        logger.info(f"Found {len(relevant_pathways)} relevant pathways")

        # Filter pathways based on relevance
        filtered_pathways = {
            pathway: conditions
            for pathway, conditions in combined_pathways.items()
            if pathway in relevant_pathways
        }
        
        if not filtered_pathways:
            logger.info("No significant pathways to plot after filtering.")
            return None

        # Identify all conditions
        all_conditions = sorted({cond for pathways in filtered_pathways.values() for cond in pathways})
        logger.info(f"Plotting data for {len(all_conditions)} conditions")

        # Prepare data for plotting
        data_records = []
        
        for pathway, conditions in filtered_pathways.items():
            for condition in all_conditions:
                if condition in conditions:
                    # Pathway is significant in this condition
                    data_records.append({
                        "Pathway": pathway,
                        "Condition": condition,
                        "Mean_Score": conditions[condition].get("Mean", 0),
                        "P-value": conditions[condition].get("P-value", 1),
                        "Significant": True,
                    })
                else:
                    # Pathway is not significant in this condition, but check if it passes size filter
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
        if data_df.empty:
            logger.warning("No data to plot after preparing records.")
            return None
            
        # Calculate mean ranks for sorting
        mean_ranks = data_df.pivot(index="Pathway", columns="Condition", values="Mean_Score").rank(
            ascending=False, method="average", axis=0
        ).mean(axis=1)
        sorted_pathways = mean_ranks.sort_values().index.tolist()
        
        # Convert Pathway column to categorical for consistent ordering
        data_df["Pathway"] = pd.Categorical(data_df["Pathway"], categories=sorted_pathways, ordered=True)
        data_df.sort_values("Pathway", inplace=True)

        # Create the plot
        plt.figure(figsize=(16, max(6, len(sorted_pathways) * 0.3)))
        ax = plt.subplot(111)
        bar_width = 0.8 / len(all_conditions)
        colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(all_conditions)))

        # Plot each condition's data
        for i, condition in enumerate(all_conditions):
            condition_data = data_df[data_df["Condition"] == condition]
            condition_data = condition_data.set_index("Pathway").reindex(sorted_pathways).reset_index()

            mean_scores = condition_data["Mean_Score"].values
            significance = condition_data["Significant"].values

            # Apply different styles based on significance
            bar_styles = [
                {"color": colors[i], "edgecolor": "black", "hatch": "//" if sig else None}
                for sig in significance
            ]
            
            # Draw the bars
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

        # Set up plot labels and titles
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

        # Adjust layout and save
        plt.tight_layout()
        output_dir = settings.plot_output_path
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file_path = output_dir / f"{settings.experiment_name}_significant_pathways_plot.png"
        
        plt.savefig(output_file_path, format="png", bbox_inches="tight", dpi=300)
        plt.close()
        
        logger.info(f"Plot saved to {output_file_path}")
        return output_file_path

    except Exception as e:
        logger.error(f"Error in plot_pathways_mean_scores: {e}", exc_info=True)
        plt.close()  # Make sure to close any open figure
        return None








