"""
tau_simulation.py
Runs a simulation of tau protein phosphorylation and aggregation over time under specified environmental conditions, and visualizes the results.
"""

import numpy as np
from TauProtein import TauProtein
from Environment import Environment
import matplotlib.pyplot as plt
import re

# --- Environment and TauProtein Setup ---
# Create an environment with specified parameters (healthy state example)
env = Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)
# Initialize a TauProtein object (4R isoform by default)
tau = TauProtein(isoform='4R')

# --- Simulation Loop ---
# Simulate tau protein state changes over 100 time steps
timepoints = np.arange(100)
tau.update_state(env, timepoints)

# --- Data Collection for Visualization ---
# Extract history data for plotting
print(f"Final aggregation state: {tau.aggregation_state}")
print(f"Phosphorylation sites phosphorylated: {sum(tau.phosphorylation_sites.values())}")
print("History:", tau.history)

ages = [entry['age'] for entry in tau.history]
phospho_counts = [entry['phospho_count'] for entry in tau.history]
aggregation_states = [entry['aggregation_state'] for entry in tau.history]

# Map aggregation state strings to numeric values for plotting
agg_state_numeric = {'monomer': 0, 'oligomer': 1, 'fibril': 2}
aggregation_numeric = [agg_state_numeric[state] for state in aggregation_states]

# --- Visualization ---
# Plot phosphorylation count and aggregation state over time
fig, ax1 = plt.subplots(figsize=(10,6))

ax1.set_xlabel('Time Step')
ax1.set_ylabel('Phosphorylation Count', color='tab:blue')
ax1.plot(ages, phospho_counts, color='tab:blue', label='Phosphorylation Count')
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax2 = ax1.twinx()
ax2.set_ylabel('Aggregation State', color='tab:red')
ax2.plot(ages, aggregation_numeric, color='tab:red', linestyle='--', label='Aggregation State')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.set_yticks([0,1,2])
ax2.set_yticklabels(['Monomer', 'Oligomer', 'Fibril'])

plt.title('Tau Phosphorylation & Aggregation Over Time')
fig.tight_layout()
plt.show()

# Example: Print TauProtein summary
print(tau)

def plot_phosphorylation_count(history):
    # Placeholder: plot number of phosphorylated sites over time
    pass

def plot_aggregation_state(history):
    # Placeholder: plot aggregation state over time
    pass

def plot_site_probabilities(probabilities, timepoints):
    import matplotlib.pyplot as plt
    for site in probabilities[0].keys():
        plt.plot(timepoints, [p[site] for p in probabilities], label=site)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylation Probability')
    plt.title('Site-specific Phosphorylation Probability Over Time')
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_phosphorylation_heatmap(probabilities, timepoints):
    # Placeholder: heatmap of site phosphorylation over time
    pass 