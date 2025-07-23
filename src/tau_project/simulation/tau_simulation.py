"""
tau_simulation.py
Runs a simulation of tau protein phosphorylation and aggregation over time under specified environmental conditions, and visualizes the results.
"""

import numpy as np
from ..models.tau_protein import TauProtein
from ..environment import Environment
from ..plot_utils import plot_tau_summary, plot_site_probabilities, plot_phosphorylation_heatmap

# --- Main Simulation and Visualization Functions ---


def run_and_plot_simulation(env=None, plot_sites=False, plot_heatmap=False):
    """
    Runs the tau protein simulation and generates summary and (optionally) site-specific visualizations.
    Args:
        env: Environment object (optional). If None, uses default parameters.
        plot_sites: If True, plot all site probability trajectories.
        plot_heatmap: If True, plot a heatmap of all site probabilities.
    """
    # Set up environment if not provided
    if env is None:
        env = Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)
    # Initialize tau protein and run simulation
    tau = TauProtein(isoform="4R")
    timepoints = np.arange(100)
    site_probabilities = tau.update_state(env, timepoints)
    # Main summary plot
    plot_tau_summary(tau.history)
    # Optional: plot all site probability trajectories
    if plot_sites:
        plot_site_probabilities(site_probabilities, timepoints)
    # Optional: plot heatmap of all site probabilities
    if plot_heatmap:
        plot_phosphorylation_heatmap(site_probabilities, timepoints)
    # Print tau protein object summary
    print(tau)


# If run as a script, show the main summary plot and advanced visualizations
if __name__ == "__main__":
    # By default, show both site-specific and heatmap visualizations
    run_and_plot_simulation(plot_sites=True, plot_heatmap=True)
