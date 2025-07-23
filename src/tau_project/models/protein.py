"""
protein.py
Defines the base Protein class, which provides basic sequence and translation logic for proteins.
"""
import numpy as np
from typing import Optional
from .aa import AminoAcid, AminoAcid_data, AminoAcid_library

class Protein:
    """
    Represents a generic protein and its sequence.
    Implements the Factory pattern for amino acid creation.
    """

    def __init__(self, name, sequence: Optional[str] = None, weight=None, length=None, organism=None, location=None, expression_level=None):
        """
        Initialize a Protein instance.
        Args:
            name (str): Name of the protein.
            sequence (str, optional): Amino acid sequence.
            weight (float, optional): Molecular weight.
            length (int, optional): Sequence length.
            organism (str, optional): Source organism.
            location (str, optional): Cellular location.
            expression_level (float, optional): Expression level.
        """
        self.name = name
        self.sequence = sequence
        self.weight = weight
        self.organism = organism
        self.location = location
        self.expression_level = expression_level
        if length is not None:
            self.length = length
        elif sequence is not None:
            self.length = len(sequence)
            self.ribosome(self.sequence)

    def bond_AA(self, aa):
        """
        Add an amino acid to the protein sequence.
        Args:
            aa (AminoAcid): Amino acid to add.
        """
        np.append(self.sequence, aa)

    def translate(self, DNA_seq, AA_lib):
        """
        Translate a DNA sequence into a protein sequence using the amino acid library.
        Args:
            DNA_seq (str or list): DNA sequence.
            AA_lib (list): Amino acid library.
        """
        AA_lib = AA_lib()
        if DNA_seq is str:
            split = [DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
        else:
            split = DNA_seq
        for codon in split:
            for aa in AA_lib:
                if codon in self.codon_list:
                    self.bond_AA(aa)
                    break

    def ribosome(self, AA_seq, AA_lib=AminoAcid_library):
        """
        Convert a sequence of identifiers into AminoAcid objects using the library.
        Args:
            AA_seq (str or list): Sequence of identifiers.
            AA_lib (list): Amino acid library.
        """
        created_seq = np.empty(self.length, dtype=AminoAcid)
        for i in range(self.length):
            iden = AA_seq[i]
            for aa in AA_lib:
                if aa.is_identifier(iden):
                    created_seq[i] = aa
                    break
        self.sequence = created_seq
