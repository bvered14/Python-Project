import numpy as np
from typing import Optional
from Protein import Protein
import AA

class TauProtein(Protein):
    def __init__(self, name="Tau", isoform="4R", sequence: Optional[np.array]=None,
                 weight=None, length=None, organism=None, location=None, expression_level=None):
        
        super().__init__(name, sequence, weight, length, organism, location, expression_level)
        
        self.isoform = isoform # add isoform defining method!
        self.phosphorylation_sites = {}  # tbd - maybe not necessary here
        self.aggregation_state = "monomer"  # or "oligomer", "fibril"
        self.truncated_site = None
        self.soluble = True

    def define_isoform(self, exon):
        if exon == "R2 at 10":
            self.isoform = "4R"
        if exon == "R1 at 10":
            self.isoform = "3R"
        return self.isoform

    def aggregation(self, ):
        if self.is_level_of_phosphorylation_healthy():
            return KeyError("This tau protein will not aggregate, it is healthy!")
        else:
            agg_score = 0

        phospho_count = len(self.phosphorylated_sites)
        agg_score += phospho_count * 1.5  # Increase weight of phospho

        if self.isoform == "4R":
            agg_score += 2

        if self.is_truncated:
            agg_score += 2

        motif_score = self.count_motifs(['VQIINK', 'VQIVYK'])
        agg_score += motif_score * 2

        if self.temperature > 37:
            agg_score += 0.5

        if self.molecular_weight > self.native_mass + 1000:
            agg_score += 1

        if agg_score >= 7:
            self.aggregation_state = 'fibril'
        elif agg_score >= 4:
            self.aggregation_state = 'oligomer'
        else:
            self.aggregation_state = 'monomer'

    
    def is_level_of_phosphorylation_healthy():
        # if phosphorylate.phosphorylate()
        return True

"""
    def phosphorylate(self, site, kinase="generic_kinase"):
        if self.sequence is None or site >= len(self.sequence):
            print(f"[!] Invalid site: {site}")
            return

        aa = self.sequence[site]
        if not getattr(aa, 'phosphorylatable', False):
            print(f"[!] {aa.name} at position {site} is not phosphorylatable.")
            return

        self.phosphorylation_sites[site] = {"state": True, "kinase": kinase}
        print(f"[+] Phosphorylated {aa.name} at site {site} by {kinase}.")
"""

tau = TauProtein()
print(tau)