import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

def plot_tau_summary(history):
    """
    Plots phosphorylation count (P > 0.5), average phosphorylation probability, and aggregation state over time.
    """
    times = [entry["minute"] for entry in history]
    phospho_counts = [entry["phospho_count"] for entry in history]
    avg_probs = [entry.get("avg_prob", 0) for entry in history]
    aggregation_states = [entry["aggregation_state"] for entry in history]
    agg_state_numeric = {"monomer": 0, "oligomer": 1, "fibril": 2}
    aggregation_numeric = [agg_state_numeric[state] for state in aggregation_states]

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set_xlabel("Time Step")
    ax1.set_ylabel("Phosphorylation Count (P > 0.5)", color="tab:blue")
    ax1.plot(times, phospho_counts, color="tab:blue", label="Phosphorylation Count")
    ax1.plot(
        times,
        avg_probs,
        color="tab:orange",
        label="Average Phosphorylation Probability",
    )
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    ax2 = ax1.twinx()
    ax2.set_ylabel("Aggregation State", color="tab:green")
    ax2.step(
        times,
        aggregation_numeric,
        color="tab:green",
        where="post",
        label="Aggregation State",
    )
    ax2.set_yticks([0, 1, 2])
    ax2.set_yticklabels(["Monomer", "Oligomer", "Fibril"])
    ax2.tick_params(axis="y", labelcolor="tab:green")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    plt.title("Tau Phosphorylation and Aggregation Over Time")
    fig.tight_layout()
    plt.show()

    print(f"Final average phosphorylation probability: {avg_probs[-1]:.3f}")
    print(f"Final aggregation state: {aggregation_states[-1]}")
    print(
        f"Maximum phosphorylation count: {max(phospho_counts)} out of {max(phospho_counts) if phospho_counts else 0} sites"
    )

def plot_site_probabilities(site_probabilities, timepoints):
    """
    Plots the probability trajectory for each phosphorylation site as a line plot.
    """
    plt.figure(figsize=(12, 6))
    for site, probs in site_probabilities.items():
        plt.plot(timepoints, probs, label=f"Site {site}")
    plt.xlabel("Time Step")
    plt.ylabel("Phosphorylation Probability")
    plt.title("Site-specific Phosphorylation Probability Over Time")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small", ncol=2)
    plt.tight_layout()
    plt.show()

def plot_phosphorylation_heatmap(site_probabilities, timepoints):
    """
    Plots a heatmap of phosphorylation probabilities for all sites over time.
    """
    df = pd.DataFrame(site_probabilities, index=timepoints)
    plt.figure(figsize=(14, 6))
    sns.heatmap(df.T, cmap="viridis", cbar_kws={"label": "Phosphorylation Probability"})
    plt.xlabel("Time Step")
    plt.ylabel("Site")
    plt.title("Phosphorylation Probability Heatmap (Sites x Time)")
    plt.tight_layout()
    plt.show() 