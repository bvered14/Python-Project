import numpy as np
from TauProtein import TauProtein
# from Environment import Environment  # Uncomment if you have this class
import matplotlib.pyplot as plt
import re

# Placeholder Environment class if not defined elsewhere
class Environment:
    def __init__(self, temperature=37, kinase_level=1.0, oxidative_stress=0.0):
        self.temperature = temperature
        self.kinase_level = kinase_level
        self.oxidative_stress = oxidative_stress

# --- Simulation Setup ---
env = Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)
tau = TauProtein(isoform='4R')

for t in range(100):
    tau.update_state(env)

# --- View final state and history ---
print(f"Final aggregation state: {tau.aggregation_state}")
print(f"Phosphorylation sites phosphorylated: {sum(tau.phosphorylation_sites.values())}")
print("History:", tau.history)

# --- Visualization ---
ages = [entry['age'] for entry in tau.history]
phospho_counts = [entry['phospho_count'] for entry in tau.history]
aggregation_states = [entry['aggregation_state'] for entry in tau.history]

agg_state_numeric = {'monomer': 0, 'oligomer': 1, 'fibril': 2}
aggregation_numeric = [agg_state_numeric[state] for state in aggregation_states]

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