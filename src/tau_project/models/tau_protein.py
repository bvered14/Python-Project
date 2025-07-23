"""
tau_protein.py
Defines the TauProtein class, which models tau protein isoforms and their molecular transitions (phosphorylation, aggregation, truncation, etc.).
"""
import random
import math
import numpy as np
from typing import Optional
from .protein import Protein
from .truncation import ProteinTruncator
from .aa import AminoAcid
import re
from ..environment import Environment as env
from collections import defaultdict

class TauProtein(Protein):
    """
    Represents a tau protein isoform and its molecular transitions.
    Extends Protein and implements additional logic for phosphorylation, aggregation, and truncation.
    Demonstrates the Strategy pattern for simulation logic.
    """
    def __init__(self, name="Tau", isoform="4R", sequence: Optional[np.array]=None,
                 weight=None, length=None, organism=None, location=None, expression_level=None):
        """
        Initialize a TauProtein instance.
        Args:
            name (str): Name of the protein.
            isoform (str): Isoform type (e.g., '4R').
            sequence (np.array, optional): Amino acid sequence.
            weight (float, optional): Molecular weight.
            length (int, optional): Sequence length.
            organism (str, optional): Source organism.
            location (str, optional): Cellular location.
            expression_level (float, optional): Expression level.
        """
        super().__init__(name, sequence, weight, length, organism, location, expression_level)
        self.S202_T205_phosphorylation = False
        self.T231_phosphorylation = False
        self.S396_S404_phosphorylation = False
        self.phosphorylation_rate = 0.1
        self.phosphatase_activity = 0.05
        self.temperature = 37
        self.aggregation_level = 0.0
        self.microtubule_binding = 1.0
        self.isoform = isoform
        self.phosphorylation_sites = {i: np.array([np.random.rand()]) for i in range(1, 80)}
        self.aggregation_state = "monomer"
        self.truncated_site = None
        self.soluble = True
        self.history = []
        self.age = 0
        self.is_truncated = False
        self.pathological = False
        self.sequence = []
        if self.sequence is None:
            self.sequence = []

    def define_isoform(self, exon):
        """
        Define the tau isoform based on exon information.
        Args:
            exon (str): Exon descriptor.
        Returns:
            str: Isoform type.
        """
        if exon == "R2 at 10":
            self.isoform = "4R"
        if exon == "R1 at 10":
            self.isoform = "3R"
        return self.isoform

    def phosphorylate(self, site):
        """
        Site-specific phosphorylation logic for key tau residues.
        Args:
            site (str): Site identifier.
        """
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
            self.phosphorylation_sites[site] = random.random() < self.phosphorylation_rate

    def dephosphorylate(self, site):
        """
        Site-specific dephosphorylation logic for key tau residues.
        Args:
            site (str): Site identifier.
        """
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
        """
        Simulate phosphorylation under different biological conditions.
        Args:
            condition (str): Biological condition ('healthy', 'pathological', 'recovery').
            temperature (float, optional): Temperature in Celsius.
        """
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
        """
        Update aggregation and microtubule binding state based on phosphorylation.
        """
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
        """
        Print the current phosphorylation and aggregation state.
        """
        print(f"S202_T205 Phosphorylation: {self.S202_T205_phosphorylation}")
        print(f"T231 Phosphorylation: {self.T231_phosphorylation}")
        print(f"S396_S404 Phosphorylation: {self.S396_S404_phosphorylation}")
        print(f"Aggregation Level: {self.aggregation_level}")
        print(f"Microtubule Binding: {self.microtubule_binding}")

    def count_phosphorylated_residues(self):
        """
        Count the number of phosphorylated sites (probability > 0.5).
        Returns:
            int: Number of phosphorylated sites.
        """
        return sum(1 for p in self.phosphorylation_sites.values() if p > 0.5)

    def detect_aggregation_motifs(self):
        """
        Detect aggregation-prone motifs in the tau sequence.
        Returns:
            int: Number of detected motifs.
        """
        sequence_str = ''.join(aa.one_letter for aa in self.sequence if hasattr(aa, 'one_letter'))
        motifs = ['VQIINK', 'VQIVYK']
        count = 0
        for motif in motifs:
            count += len(re.findall(motif, sequence_str))
        return count

    def compute_aggregation_score(self):
        """
        Compute an aggregation score based on phosphorylation and motifs.
        Returns:
            float: Aggregation score.
        """
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
        """
        Update the aggregation state (monomer, oligomer, fibril) based on score.
        """
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
        Populates self.history with a dict for each timepoint.
        Args:
            environment (Environment): Simulation environment.
            timepoints (np.array): Array of timepoints.
        Returns:
            dict: Site probabilities over time.
        """
        from collections import defaultdict
        self.history = []  # Reset history at the start
        site_probabilities = {site: [float(self.phosphorylation_sites[site][0])] for site in self.phosphorylation_sites}
        temp_effect = self.check_temp(environment)
        kinase_effect = self.check_kinase(environment)
        phosphatase_effect = self.check_phosphatase(environment)
        protease_effect = self.check_protease(environment)
        oxidative_effect = self.check_oxidative_stress(environment)
        k_p = temp_effect * kinase_effect * oxidative_effect
        k_d = phosphatase_effect * protease_effect
        for i, time in enumerate(timepoints):
            for site in site_probabilities:
                prev = site_probabilities[site][-1]
                if i == 0:
                    continue  # already initialized
                P_new = prev + (k_p * (1 - prev) - k_d * prev)
                if P_new < 0:
                    P_new = 0.0
                elif P_new > 1:
                    P_new = 1.0
                site_probabilities[site].append(P_new)
            probs_now = [site_probabilities[site][i] for site in site_probabilities]
            phospho_count = sum(p > 0.5 for p in probs_now)
            avg_prob = float(np.mean(probs_now))
            self.update_aggregation_state()
            entry = {
                'age': int(time),
                'minute': int(time),
                'phospho_count': phospho_count,
                'aggregation_state': self.aggregation_state,
                'avg_prob': avg_prob
            }
            self.history.append(entry)
        return site_probabilities

    def check_temp(self, environment):
        """
        Check the effect of temperature on the simulation.
        Args:
            environment (Environment): Simulation environment.
        Returns:
            float: Temperature effect multiplier.
        """
        healthy_temp_range = (36, 38)
        if healthy_temp_range[0] <= environment.temperature <= healthy_temp_range[1]:
            return 1.0
        elif environment.temperature < healthy_temp_range[0]:
            return 0.7
        else:
            return 1.3

    def check_kinase(self, environment):
        """
        Check the effect of kinase level on the simulation.
        Args:
            environment (Environment): Simulation environment.
        Returns:
            float: Kinase effect multiplier.
        """
        healthy_kinase_range = (0.8, 1)
        if healthy_kinase_range[0] <= environment.temperature <= healthy_kinase_range[1]:
            return 1.0
        else:
            return 0.05 * environment.kinase_level

    def check_phosphatase(self, environment):
        """
        Check the effect of phosphatase level on the simulation.
        Args:
            environment (Environment): Simulation environment.
        Returns:
            float: Phosphatase effect multiplier.
        """
        phosphatase_range = (0.8, 1)
        if phosphatase_range[0] <= environment.kinase_level <= phosphatase_range[1]:
            return environment.kinase_level
        else:
            return 0.02 * environment.phosphatase_level

    def check_protease(self, environment):
        """
        Check the effect of protease level on the simulation.
        Args:
            environment (Environment): Simulation environment.
        Returns:
            float: Protease effect multiplier.
        """
        if 0.4 < environment.protease_level < 0.8:
            return 0.7
        elif environment.protease_level <= 0.4:
            return 0.2
        elif environment.protease_level <= 1:
            return 1.5
        return 1.0

    def check_oxidative_stress(self, environment):
        """
        Check the effect of oxidative stress on the simulation.
        Args:
            environment (Environment): Simulation environment.
        Returns:
            float: Oxidative stress effect multiplier.
        """
        if 0.4 < environment.oxidative_stress < 0.8:
            return 0.7
        elif environment.oxidative_stress <= 0.4:
            return 0.2
        elif environment.oxidative_stress <= 1:
            return 1.5
        return 1.0

    def truncate(self, site):
        """
        Truncate the tau protein sequence at a given site.
        Args:
            site (str): Truncation site (e.g., 'D421').
        Raises:
            ValueError: If no sequence is present.
        """
        if self.sequence is not None:
            truncated_seq, trunc_aa = ProteinTruncator.truncate(self.sequence, site)
            self.sequence = truncated_seq
            self.truncated_site = site
            self.is_truncated = True
            self.truncation_aa = trunc_aa
        else:
            raise ValueError("No sequence to truncate.")
