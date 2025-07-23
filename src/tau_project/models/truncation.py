"""
truncation.py
Provides the ProteinTruncator class for truncating protein sequences at specific sites.
Implements the Static Method pattern for utility functionality.
"""
from .aa import AminoAcid

class ProteinTruncator:
    """
    Provides static methods for truncating protein sequences at specific sites.
    Demonstrates the Static Method pattern.
    """
    @staticmethod
    def truncate(sequence, site):
        """
        Truncate a protein sequence at the given site (e.g., 'D421').
        Args:
            sequence (str or list): Protein sequence as a string or list of AminoAcid objects.
            site (str): Truncation site, e.g., 'D421' (residue and position).
        Returns:
            tuple: (truncated sequence, AminoAcid at truncation site)
        Raises:
            ValueError: If the residue at the site does not match.
            TypeError: If the sequence is not a string or list.
        """
        residue = site[0]
        pos = int(site[1:])
        if isinstance(sequence, str):
            if sequence[pos-1] != residue:
                raise ValueError(f"Residue at position {pos} is not {residue}")
            trunc_aa = AminoAcid(
                name="Aspartic acid" if residue == 'D' else f"Residue {residue}",
                three_letter="Asp" if residue == 'D' else "",
                one_letter=residue,
                polarity="polar" if residue == 'D' else "",
                charge="-1" if residue == 'D' else "",
                r_group="CH2COOH" if residue == 'D' else "",
                codon_list=["GAU", "GAC"] if residue == 'D' else []
            )
            return sequence[:pos], trunc_aa
        elif isinstance(sequence, list):
            if sequence[pos-1].one_letter != residue:
                raise ValueError(f"Residue at position {pos} is not {residue}")
            return sequence[:pos], sequence[pos-1]
        else:
            raise TypeError("Sequence must be a string or list of AminoAcid objects.")
