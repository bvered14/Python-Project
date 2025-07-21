import numpy as np
from typing import Optional
import AA

class Protein:
    AminoAcid_data= [
    {
        "name": "Glycine",
        "three_letter": "Gly",
        "one_letter": "G",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "H",
        "codon_list": ["GGT", "GGC", "GGA", "GGG"],
        "pI": 6.06,
        "volume": 60.1,
        "weight": 75.07,
    },
    {
        "name": "Alanine",
        "three_letter": "Ala",
        "one_letter": "A",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CH3",
        "codon_list": ["GCT", "GCC", "GCA", "GCG"],
        "pI": 6.11,
        "volume": 88.6,
        "weight": 89.09,
    },
    {
        "name": "Valine",
        "three_letter": "Val",
        "one_letter": "V",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CHCH3CH3",
        "codon_list": ["GTT", "GTC", "GTA", "GTG"],
        "pI": 6.00,
        "volume": 140.0,
        "weight": 117.15,
    },
    {
        "name": "Leucine",
        "three_letter": "Leu",
        "one_letter": "L",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CH2CHCH3CH3",
        "codon_list": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "pI": 6.01,
        "volume": 166.7,
        "weight": 131.17,
    },
    {
        "name": "Isoleucine",
        "three_letter": "Ile",
        "one_letter": "I",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CHCH3CH2CH3",
        "codon_list": ["ATT", "ATC", "ATA"],
        "pI": 6.05,
        "volume": 168.8,
        "weight": 131.17,
    },
    {
        "name": "Methionine",
        "three_letter": "Met",
        "one_letter": "M",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CH2CH2SCH3",
        "codon_list": ["ATG"],
        "pI": 5.74,
        "volume": 162.9,
        "weight": 149.21,
    },
    {
        "name": "Phenylalanine",
        "three_letter": "Phe",
        "one_letter": "F",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CH2C6H5",
        "codon_list": ["TTT", "TTC"],
        "pI": 5.48,
        "volume": 189.9,
        "weight": 165.19,
    },
    {
        "name": "Tryptophan",
        "three_letter": "Trp",
        "one_letter": "W",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CH2C8H6N",
        "codon_list": ["TGG"],
        "pI": 5.89,
        "volume": 227.8,
        "weight": 204.23,
    },
    {
        "name": "Proline",
        "three_letter": "Pro",
        "one_letter": "P",
        "polarity": "nonpolar",
        "charge": 0,
        "r_group": "CH2CH2CH2 (cyclic to backbone N)",
        "codon_list": ["CCT", "CCC", "CCA", "CCG"],
        "pI": 6.30,
        "volume": 112.7,
        "weight": 115.13,
    },
    {
        "name": "Serine",
        "three_letter": "Ser",
        "one_letter": "S",
        "polarity": "polar",
        "charge": 0,
        "r_group": "CH2OH",
        "codon_list": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "pI": 5.68,
        "volume": 89.0,
        "weight": 105.09,
    },
    {
        "name": "Threonine",
        "three_letter": "Thr",
        "one_letter": "T",
        "polarity": "polar",
        "charge": 0,
        "r_group": "CHOHCH3",
        "codon_list": ["ACT", "ACC", "ACA", "ACG"],
        "pI": 5.60,
        "volume": 116.1,
        "weight": 119.12,
    },
    {
        "name": "Cysteine",
        "three_letter": "Cys",
        "one_letter": "C",
        "polarity": "polar",
        "charge": 0,
        "r_group": "CH2SH",
        "codon_list": ["TGT", "TGC"],
        "pI": 5.07,
        "volume": 108.5,
        "weight": 121.15,
    },
    {
        "name": "Tyrosine",
        "three_letter": "Tyr",
        "one_letter": "Y",
        "polarity": "polar",
        "charge": 0,
        "r_group": "CH2C6H4OH",
        "codon_list": ["TAT", "TAC"],
        "pI": 5.66,
        "volume": 193.6,
        "weight": 181.19,
    },
    {
        "name": "Asparagine",
        "three_letter": "Asn",
        "one_letter": "N",
        "polarity": "polar",
        "charge": 0,
        "r_group": "CH2CONH2",
        "codon_list": ["AAT", "AAC"],
        "pI": 5.41,
        "volume": 114.1,
        "weight": 132.12,
    },
    {
        "name": "Glutamine",
        "three_letter": "Gln",
        "one_letter": "Q",
        "polarity": "polar",
        "charge": 0,
        "r_group": "CH2CH2CONH2",
        "codon_list": ["CAA", "CAG"],
        "pI": 5.65,
        "volume": 143.8,
        "weight": 146.15,
    },
    {
        "name": "Aspartic Acid",
        "three_letter": "Asp",
        "one_letter": "D",
        "polarity": "acidic",
        "charge": -1,
        "r_group": "CH2COOH",
        "codon_list": ["GAT", "GAC"],
        "pI": 2.77,
        "volume": 111.1,
        "weight": 133.10,
    },
    {
        "name": "Glutamic Acid",
        "three_letter": "Glu",
        "one_letter": "E",
        "polarity": "acidic",
        "charge": -1,
        "r_group": "CH2CH2COOH",
        "codon_list": ["GAA", "GAG"],
        "pI": 3.22,
        "volume": 138.4,
        "weight": 147.13,
    },
    {
        "name": "Lysine",
        "three_letter": "Lys",
        "one_letter": "K",
        "polarity": "basic",
        "charge": +1,
        "r_group": "CH2CH2CH2CH2NH2",
        "codon_list": ["AAA", "AAG"],
        "pI": 9.74,
        "volume": 171.3,
        "weight": 146.19,
    },
    {
        "name": "Arginine",
        "three_letter": "Arg",
        "one_letter": "R",
        "polarity": "basic",
        "charge": +1,
        "r_group": "CH2CH2CH2NHCNHNH2",
        "codon_list": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "pI": 10.76,
        "volume": 202.0,
        "weight": 174.20,
    },
    {
        "name": "Histidine",
        "three_letter": "His",
        "one_letter": "H",
        "polarity": "basic",
        "charge": +1,
        "r_group": "CH2C3H3N2",
        "codon_list": ["CAT", "CAC"],
        "pI": 7.59,
        "volume": 153.2,
        "weight": 155.16,
    },
]

    AminoAcid_library=[]
    for aa in AminoAcid_data:
        amino_acid = AA(
            name=aa["name"],
            three_letter=aa["three_letter"],
            one_letter=aa["one_letter"],
            polarity=aa["polarity"],
            charge=aa["charge"],
            r_group=aa["r_group"],
            codon_list=aa["codon_list"],
            pKa=aa.get("pKa"),
            volume=aa.get("volume"),
            weight=aa.get("weight"),
            PTM=aa.get("PTM", False)
        )
        AminoAcid_library.append(amino_acid)

    def __init__(self, name, sequence: Optional[str]=None, weight=None, length=None, organism=None, location=None, expression_level=None):
        self.name = name
        self.weight = weight
        self.organism = organism
        self.location = location
        self.expression_level = expression_level

        # Compute or assign length
        if length is not None:
            self.length = length
        elif sequence is not None:
            self.length = len(sequence)
            self.ribosome(self.sequence)

    def bond_AA(self, aa):
        np.append(self.sequence,aa)
   
    def translate(self,DNA_seq,AA_lib=AminoAcid_library):
        if DNA_seq is str:
            split=[DNA_seq[i:i+3] for i in range(0, len(DNA_seq), 3)]
        else:
            split=DNA_seq
        for codon in split:
            for aa in AA_lib:
                if codon in self.codon_list:
                    self.bond_AA(aa)
                    break

    def ribosome(self,AA_seq, AA_lib=AminoAcid_library):
        created_seq=np.empty(self.length,dtype=AA)
        for i in range(self.length):
            iden=AA_seq[i]
            for aa in AA_lib:
                if aa.is_identifier(iden):
                    created_seq[i]=aa
                    break
        self.sequence=created_seq
    
   