import numpy as np
from typing import Optional
import AA

class Protein:
    def __init__(self, name, sequence: Optional[np.array]=None, weight=None, length=None, organism=None, location=None, expression_level=None):
        self.name = name
        self.sequence = sequence
        self.weight = weight
        self.organism = organism
        self.location = location
        self.expression_level = expression_level

        # Compute or assign length
        if length is not None:
            self.length = length
        elif sequence is not None:
            self.length = len(sequence)
    
    def bond_AA(self, aa):
        np.append(self.sequence,aa)

    def translate(self,DNA_seq,AA_lib):
        AA_lib=AA_lib()
        if DNA_seq is str:
            split=[DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
        else:
            split=DNA_seq
        for codon in split:
            for aa in AA_lib:
                if codon in self.codon_list:
                    self.bond_AA(aa)
                    break