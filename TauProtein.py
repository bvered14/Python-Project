"""
TauProtein.py
Defines the TauProtein class, which models tau protein isoforms and their molecular transitions (phosphorylation, aggregation, truncation, etc.).
"""
import random
import math
import numpy as np
from typing import Optional
from Protein import Protein
from truncation import ProteinTruncator
from AA import AminoAcid
import re
from Environment import Environment as env

# TauProtein class models tau-specific logic, including isoforms, phosphorylation, aggregation, truncation, and simulation state
class TauProtein(Protein):
    def __init__(self, name="Tau", isoform="4R", sequence: Optional[np.array]=None,
                 weight=None, length=None, organism=None, location=None, expression_level=None):
        # Initialize tau protein properties, including isoform, phosphorylation, aggregation, and history
        super().__init__(name, sequence, weight, length, organism, location, expression_level)
        # Site-specific phosphorylation booleans
        self.S202_T205_phosphorylation = False
        self.T231_phosphorylation = False
        self.S396_S404_phosphorylation = False
        self.phosphorylation_rate = 0.1
        self.phosphatase_activity = 0.05
        self.temperature = 37
        self.aggregation_level = 0.0
        self.microtubule_binding = 1.0
        # Advanced/General attributes
        self.isoform = isoform
        self.phosphorylation_sites = {}  # e.g., {site: True/False}
        self.aggregation_state = "monomer"  # or "oligomer", "fibril"
        self.truncated_site = None
        self.soluble = True
        self.history = []
        self.age = 0
        self.is_truncated = False
        self.pathological = False
        self.sequence = []
        # If sequence is not provided, initialize as empty list
        if self.sequence is None:
            self.sequence = []

    def define_isoform(self, exon):
        # Define the tau isoform based on exon information
        if exon == "R2 at 10":
            self.isoform = "4R"
        if exon == "R1 at 10":
            self.isoform = "3R"
        return self.isoform

    # --- Site-specific phosphorylation logic ---
    def phosphorylate(self, site):
        # Site-specific phosphorylation logic for key tau residues
        if site == "S202_T205":
            self.S202_T205_phosphorylation = random.random() < self.phosphorylation_rate
            self.phosphorylation_sites[site] = self.S202_T205_phosphorylation
        elif site == "T231":
            if self.S202_T205_phosphorylation:
                self.T231_phosphorylation = random.random() < (self.phosphorylation_rate * 1.5)
            else:
                self.T231_phosphorylation = random.random() < self.phosphorylation_rate
            self.phosphorylation_sites[site] = self.T231_phosphorylation
        elif site == "S396_S404":
            if self.T231_phosphorylation:
                self.S396_S404_phosphorylation = random.random() < (self.phosphorylation_rate * 1.5)
            else:
                self.S396_S404_phosphorylation = random.random() < self.phosphorylation_rate
            self.phosphorylation_sites[site] = self.S396_S404_phosphorylation
        else:
            # General site logic
            self.phosphorylation_sites[site] = random.random() < self.phosphorylation_rate

    def dephosphorylate(self, site):
        # Site-specific dephosphorylation logic for key tau residues
        if site == "S202_T205":
            if random.random() < self.phosphatase_activity:
                self.S202_T205_phosphorylation = False
                self.phosphorylation_sites[site] = False
        elif site == "T231":
            if random.random() < self.phosphatase_activity:
                self.T231_phosphorylation = False
                self.phosphorylation_sites[site] = False
        elif site == "S396_S404":
            if random.random() < self.phosphatase_activity:
                self.S396_S404_phosphorylation = False
                self.phosphorylation_sites[site] = False
        else:
            if random.random() < self.phosphatase_activity:
                self.phosphorylation_sites[site] = False

    def simulate_phosphorylation(self, condition, temperature=None):
        # Simulate phosphorylation under different biological conditions
        if temperature is not None:
            self.temperature = temperature
        temperature_factor = 1 + (self.temperature - 37) * 0.01
        self.phosphorylation_rate *= temperature_factor
        if condition == "healthy":
            self.phosphorylation_rate = 0.1
            self.phosphatase_activity = 0.05
        elif condition == "pathological":
            self.phosphorylation_rate = 0.5
            self.phosphatase_activity = 0.01
        elif condition == "recovery":
            self.phosphorylation_rate = 0.05
            self.phosphatase_activity = 0.1
        self.phosphorylate("S202_T205")
        self.phosphorylate("T231")
        self.phosphorylate("S396_S404")
        self.update_aggregation_and_binding()

    def update_aggregation_and_binding(self):
        # Update aggregation and microtubule binding state based on phosphorylation
        total_phosphorylation = sum([self.S202_T205_phosphorylation, self.T231_phosphorylation, self.S396_S404_phosphorylation])
        if total_phosphorylation >= 0.7:
            self.aggregation_level = 1.0
            self.microtubule_binding = 0.0
        elif total_phosphorylation >= 0.5:
            self.aggregation_level = 0.7
            self.microtubule_binding = 0.3
        elif total_phosphorylation >= 0.3:
            self.aggregation_level = 0.3
            self.microtubule_binding = 0.6
        else:
            self.aggregation_level = 0.0
            self.microtubule_binding = 1.0

    def display_phosphorylation(self):
        # Print the current phosphorylation and aggregation state
        print(f"S202_T205 Phosphorylation: {self.S202_T205_phosphorylation}")
        print(f"T231 Phosphorylation: {self.T231_phosphorylation}")
        print(f"S396_S404 Phosphorylation: {self.S396_S404_phosphorylation}")
        print(f"Aggregation Level: {self.aggregation_level}")
        print(f"Microtubule Binding: {self.microtubule_binding}")

    # --- Advanced/General aggregation logic ---
    def count_phosphorylated_residues(self):
        # Count the number of phosphorylated residues in the sequence
        return sum(1 for aa in self.sequence if hasattr(aa, 'PTM') and aa.PTM == "Phospho")
    
    def detect_aggregation_motifs(self):
        # Detect aggregation-prone motifs in the tau sequence
        sequence_str = ''.join(aa.one_letter for aa in self.sequence if hasattr(aa, 'one_letter'))
        motifs = ['VQIINK', 'VQIVYK']
        count = 0
        for motif in motifs:
            count += len(re.findall(motif, sequence_str))
        return count
    
    def compute_aggregation_score(self):
        # Compute an aggregation score based on phosphorylation and motifs
        phospho_count = self.count_phosphorylated_residues()
        motif_count = self.detect_aggregation_motifs()
        score = 0
        score += phospho_count * 1.5
        score += motif_count * 2
        if self.is_truncated:
            score += 2
        if self.isoform == '4R':
            score += 1
        return score
    
    def update_aggregation_state(self):
        # Update the aggregation state (monomer, oligomer, fibril) based on score
        score = self.compute_aggregation_score()
        if score >= 7:
            self.aggregation_state = 'fibril'
        elif score >= 4:
            self.aggregation_state = 'oligomer'
        else:
            self.aggregation_state = 'monomer'

    def update_state(self, environment, timepoints: np.array):
        """
        Main orchestrator: updates tau protein state over a series of timepoints by checking temperature, kinase, phosphatase, protease, and oxidative stress effects.
        """
        probabilities = []
        for t in timepoints:
            site_probs = {}
            # Check all effects
            temp_effect = self.check_temp(environment)
            kinase_effect = self.check_kinase(environment)
            phosphatase_effect = self.check_phosphatase(environment)
            protease_effect = self.check_protease(environment)
            oxidative_effect = self.check_oxidative_stress(environment)
            # Combine effects for each site
            for site in self.phosphorylation_sites:
                # Combine k_p and k_d effects multiplicatively
                k_p = temp_effect * kinase_effect * oxidative_effect
                k_d = phosphatase_effect * protease_effect
                P = 1 if self.phosphorylation_sites[site] else 0
                delta_P = (k_p * (1 - P) - k_d * P) * (t if t > 0 else 1)
                prob = P + delta_P
                site_probs[site] = min(max(prob, 0), 1)
            probabilities.append(site_probs)
            # Update state for each site
            for site, prob in site_probs.items():
                self.phosphorylation_sites[site] = random.random() < prob
            # Log state at this timepoint
            self.history.append({
                'age': t,
                'phospho_count': self.count_phosphorylated_residues(),
                'aggregation_state': self.aggregation_state,
                'is_truncated': self.is_truncated,
                'pathological': self.pathological
            })
        return probabilities
      
    def check_temp(self, environment):
        healthy_temp_range = (36, 38)
        if healthy_temp_range[0] <= environment.temperature <= healthy_temp_range[1]:
            return 1.0
        elif environment.temperature < healthy_temp_range[0]:
            return 0.7
        else:
            return 1.3

    def check_kinase(self, environment):
        return 0.05 * environment.kinase_level

    def check_phosphatase(self, environment):
        return 0.02 * environment.phosphatase_level

    def check_protease(self, environment):
        if 0.4 < environment.protease_level < 0.8:
            return 0.7
        elif environment.protease_level <= 0.4:
            return 0.2
        elif environment.protease_level <= 1:
            return 1.5
        return 1.0

    def check_oxidative_stress(self, environment):
        if 0.4 < environment.oxidative_stress < 0.8:
            return 0.7
        elif environment.oxidative_stress <= 0.4:
            return 0.2
        elif environment.oxidative_stress <= 1:
            return 1.5
        return 1.0

    def truncate(self, site):
        # Truncate the tau protein sequence at a given site
        if self.sequence is not None:
            truncated_seq, trunc_aa = ProteinTruncator.truncate(self.sequence, site)
            self.sequence = truncated_seq
            self.truncated_site = site
            self.is_truncated = True
            self.truncation_aa = trunc_aa
        else:
            raise ValueError("No sequence to truncate.") 
        
tau = TauProtein()
envnrmt = env(temperature = 38)
timestamps = np.array([0, 30, 60])
probs = tau.update_state(envnrmt, timestamps)    
print(probs)

    # def check_temp(self, environment):
    # effects = {}
    # for site in self.phosphorylation_sites:
    #     if environment.temperature < 36:
    #         effects[site] = {'k_d': 0.5}  # 50% dephosphorylation
    #     elif environment.temperature > 38:
    #         effects[site] = {'k_d': 1.2}  # 20% faster cleanup
    #     else:
    #         effects[site] = {'k_d': 1.0}
    # return effects

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import Environment as env

    # Example environment and tau protein
    environment = env.Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)
    tau = TauProtein(isoform='4R')

    # Simulate over 100 timepoints
    timepoints = np.arange(1, 101)
    tau.update_state(environment, timepoints)

    # Example: View final state and history
    print(f"Final aggregation state: {tau.aggregation_state}")
    print(f"Phosphorylation sites phosphorylated: {sum(tau.phosphorylation_sites.values())}")
    print("History:", tau.history)

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

    tau = TauProtein()
    print(tau)
