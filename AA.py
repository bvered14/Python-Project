import math
import re

class AminoAcid:
    def __init__(self, name, three_letter, one_letter, polarity, charge, r_group, codon_list, pKa=None, volume=None, weight=None, PTM=False):
        self.name = name
        self.three_letter = three_letter
        self.one_letter = one_letter
        self.polarity = polarity
        self.charge = charge
        self.r_group = r_group
        self.codon_list = codon_list
        self.pKa = pKa
        self.volume = volume
        self.PTM = PTM

        if weight is not None:
            self.weight = weight
        else:
            self.weight = self.calculate_weight()

    def __get_charge__(self, pH=None):
        if pH is not None and self.pKa is not None:
            change = int(math.copysign(1, self.pKa - pH)) if self.pKa != pH else 1
            self.charge=self.charge-(change<0)
        return self.charge  
    
    def __has_codon__(self, codon):
        return codon.upper() in self.codon_list
    
    def __is_identifier__(self, identifier):

        identifier = identifier.upper()
        return identifier in {self.name.upper(), self.three_letter.upper(), self.one_letter.upper()}
    
    def calculate_weight(self):

        atomic_weights = {
            'C': 12.01,
            'H': 1.008,
            'N': 14.01,
            'O': 16.00,
            'S': 32.06,
            'P': 30.97,
        }
        weight = 74.06
        match self.PTM:
            case "Ubi":
                self.r_group=self.r_group[:-4]
                weight+=8565
            case "GlcNAc":
                self.r_group=self.r_group[:-7]
                weight+=203
        cleaned = re.sub(r'[^A-Za-z0-9]', '', self.r_group)
        pattern = r'([A-Z][a-z]*)(\d*)'
        counts = {}
        
        for (elem, count) in re.findall(pattern, cleaned):
            count = int(count) if count else 1
            counts[elem] = counts.get(elem, 0) + count
        
        
        for atom, cnt in counts.items():
            w = atomic_weights.get(atom)
            if w is None:
                raise ValueError(f"Unknown element '{atom}' in formula")
            weight += w * cnt
        self.weight=weight
        return weight
    
    def add_PTM(self, modification):
        if self.one_letter in {"A","V","L","I","F","W"}:
            raise ("Amino acid does not undergo PTM")
        match modification:
            case "Phosphorylation"|"p"|"Phospo"|"P":
                if self.three_letter not in {"Ser", "Thr", "Tyr", "His", "Asp", "Glu", "Arg", "Lys", "Cys"}:
                    raise ({f"{self.name}does not undergo Phosphorylation"})
                self.r_group=self.r_group[:-1]+"-PO3"
                self.PTM="Phospho"
                self.calculate_weight()
            case "Acetylation"|"a"|"A"|"Acetyl":
                if self.three_letter not in {"Lys", "Met"}:
                    raise ({f"{self.name}does not undergo Acetylation"})
                self.r_group=self.r_group[:-1]+"-COCH3"
                self.PTM="Acetyl"
                self.calculate_weight()
            case "Methylation"|"m"|"M"|"Methyl":
                if self.three_letter not in {"Lys", "Arg", "His"}:
                    raise ({f"{self.name}does not undergo Methylation"})
                self.r_group=self.r_group[:-1]+"-CH3"
                self.PTM="Methyl"
                self.calculate_weight()
            case "Ubiquitination"|"u"|"U"|"Ubi":
                if self.three_letter not in {"Lys", "Met"}:
                    raise ({f"{self.name}does not undergo Ubiquitination"})
                self.r_group=self.r_group[:-1]+"-UBI"
                self.PTM="Ubi"
                self.calculate_weight()
            case "O-GlcNAcylation"|"O-Glc"|"GlcNAc":
                if self.three_letter not in {"Ser", "Thr"}:
                    raise ({f"{self.name}does not undergo O-GlcNAcylation"})
                self.r_group=self.r_group[:-1]+"-GlcNAc"
                self.PTM="GlcNAc"
                self.calculate_weight()

    def __str__(self):
        return (
            f"Amino Acid: {self.name} "
            f"({self.one_letter}, {self.three_letter})\n"
            f"  Polarity: {self.polarity}\n"
            f"  Charge: {self.charge}\n"
            f"  R-group: {self.r_group}\n"
            f"  Codons: {', '.join(self.codon_list)}\n"
            f"  pKa: {self.pKa if self.pKa is not None else 'N/A'}\n"
            f"  Volume: {self.volume if self.volume is not None else 'N/A'}\n"
        )
        
