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

tau = TauProtein()
print(tau)