# Tau Protein Simulator

This Python-based simulator models the behavior of tau protein isoforms and their key molecular transitions under biological conditions. It simulates phosphorylation, aggregation, microtubule binding, and truncation, providing visual outputs for research and education.

---

## What the Simulator Does

- Models tau protein isoforms (`3R`, `4R`) and their behavior over time.
- Simulates molecular events:
  - **Phosphorylation** at multiple sites
  - **Aggregation** into oligomers and fibrils
  - **Binding** and unbinding to microtubules
  - **Truncation** at defined amino acid sites
- Allows customization of biological environment (e.g., kinase activity, oxidative stress).
- Tracks changes in state and visualizes dynamics with plots.

---

## How to Run It

1. **Set up the environment and tau protein:**

from TauProtein import TauProtein
from Environment import Environment
import numpy as np

# Define biological conditions
env = Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)

# Create a tau protein instance
tau = TauProtein(isoform='4R')

# Define time steps
timepoints = np.arange(100)

# Simulate
tau.update_state(env, timepoints)

2. **Visualize the results:**

Run the visualization script (`tau_simulation.py`) to plot phosphorylation counts and aggregation state changes over time.

---

## Code Structure and Key Functions

- **`TauProtein.py`**  
  Defines the `TauProtein` class, which models tau protein behavior.

  **Key Methods:**
  - `phosphorylate(site)` – Phosphorylates the tau protein at a specific site.
  - `aggregate()` – Transitions the protein into aggregated states based on phosphorylation.
  - `bind_microtubule()` – Toggles the microtubule-bound state.
  - `truncate(site)` – Truncates the tau protein at a given residue position.
  - `update_state(environment, timepoints)` – Runs the simulation over time.

- **`Environment.py`**  
  Defines the `Environment` class, modeling biological conditions.

  **Attributes:**
  - `temperature`, `kinase_level`, `phosphatase_level`, `protease_level`, `oxidative_stress`

- **`truncation.py`**  
  Contains the `ProteinTruncator` class.

  **Key Function:**
  - `truncate(sequence, site)` – Truncates the protein sequence at a given site (e.g., D421).

- **`tau_simulation.py`**  
  Script that initializes the simulation, runs it, and generates visual output.

- **`AA.py`**  
  Defines the `AminoAcid` class with biochemical properties:
  - `name`, `three_letter`, `one_letter`, `polarity`, `charge`, `r_group`, `codon_list`

---

## References

1. Chang, C. W., Shao, E., & Mucke, L. (2021). Tau: enabler of diverse brain disorders and target of rapidly evolving therapeutic strategies. Science 371: eabb8255.

2. Cisek, K., L Cooper, G., J Huseby, C., & Kuret, J. (2014). Structure and mechanism of action of tau aggregation inhibitors. Current Alzheimer Research, 11(10), 918-927.

3. Hernández, F., Ferrer, I., Pérez, M., Zabala, J. C., Del Rio, J. A., & Avila, J. (2023). Tau aggregation. Neuroscience, 518, 64-69.

4. Lane‐Donovan, C., Smith, A. W., Saloner, R., Miller, B. L., Casaletto, K. B., & Kao, A. W. (2025). Tau phosphorylation at Alzheimer's disease biomarker sites impairs its cleavage by lysosomal proteases. Alzheimer's & Dementia, 21(6), e70320.

5. Ren, Q. G., Liao, X. M., Chen, X. Q., Liu, G. P., & Wang, J. Z. (2007). Effects of tau phosphorylation on proteasome activity. FEBS letters, 581(7), 1521-1528.

6. Alavi Naini, S. M., & Soussi-Yanicostas, N. (2015). Tau hyperphosphorylation and oxidative stress, a critical vicious circle in neurodegenerative tauopathies?. Oxidative medicine and cellular longevity, 2015(1), 151979.

7. Salazar, C., & Höfer, T. (2006). Kinetic models of phosphorylation cycles: a systematic approach using the rapid-equilibrium approximation for protein–protein interactions. Biosystems, 83(2-3), 195-206.