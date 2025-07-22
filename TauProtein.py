import random
import math
import numpy as np
from typing import Optional
from Protein import Protein
from truncation import ProteinTruncator
from AA import AminoAcid
import re

class TauProtein(Protein):
    def __init__(self, name="Tau", isoform="4R", sequence: Optional[np.array]=None,
                 weight=None, length=None, organism=None, location=None, expression_level=None):
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
        # If sequence is not provided, initialize as empty list
        if self.sequence is None:
            self.sequence = []

    def define_isoform(self, exon):
        if exon == "R2 at 10":
            self.isoform = "4R"
        if exon == "R1 at 10":
            self.isoform = "3R"
        return self.isoform

    # --- Site-specific phosphorylation logic ---
    def phosphorylate(self, site):
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
        print(f"S202_T205 Phosphorylation: {self.S202_T205_phosphorylation}")
        print(f"T231 Phosphorylation: {self.T231_phosphorylation}")
        print(f"S396_S404 Phosphorylation: {self.S396_S404_phosphorylation}")
        print(f"Aggregation Level: {self.aggregation_level}")
        print(f"Microtubule Binding: {self.microtubule_binding}")

    # --- Advanced/General aggregation logic ---
    def count_phosphorylated_residues(self):
        return sum(1 for aa in self.sequence if hasattr(aa, 'PTM') and aa.PTM == "Phospho")
    
    def detect_aggregation_motifs(self):
        sequence_str = ''.join(aa.one_letter for aa in self.sequence if hasattr(aa, 'one_letter'))
        motifs = ['VQIINK', 'VQIVYK']
        count = 0
        for motif in motifs:
            count += len(re.findall(motif, sequence_str))
        return count
    
    def compute_aggregation_score(self):
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
        score = self.compute_aggregation_score()
        if score >= 7:
            self.aggregation_state = 'fibril'
        elif score >= 4:
            self.aggregation_state = 'oligomer'
        else:
            self.aggregation_state = 'monomer'

    def update_state(self, environment, time_step=1):
        # Placeholder: implement these methods as needed
        if hasattr(self, 'apply_phosphorylation'):
            self.apply_phosphorylation(environment)
        if hasattr(self, 'apply_dephosphorylation'):
            self.apply_dephosphorylation(environment)
        if hasattr(self, 'apply_truncation'):
            self.apply_truncation(environment)
        self.update_aggregation_state()
        if hasattr(self, 'apply_clearance'):
            self.apply_clearance(environment)
        self.age += time_step
        # Log history
        self.history.append({
            'age': self.age,
            'phospho_count': self.count_phosphorylated_residues(),
            'aggregation_state': self.aggregation_state,
            'is_truncated': self.is_truncated,
            'pathological': self.pathological
        })

    def truncate(self, site):
        if self.sequence is not None:
            truncated_seq, trunc_aa = ProteinTruncator.truncate(self.sequence, site)
            self.sequence = truncated_seq
            self.truncated_site = site
            self.is_truncated = True
            self.truncation_aa = trunc_aa
        else:
            raise ValueError("No sequence to truncate.") 