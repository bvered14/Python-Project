from AA import AminoAcid

class ProteinTruncator:
    @staticmethod
    def truncate(sequence, site):
        """
        Truncate a protein sequence at the given site (e.g., 'D421').
        sequence: a string or list of AminoAcid objects
        site: string, e.g., 'D421' (residue and position)
        Returns the truncated sequence and the AminoAcid at the truncation site.
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