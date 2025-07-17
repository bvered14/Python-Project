import math

class AminoAcid:
    def __init__(self, name, three_letter, one_letter, polarity, charge, r_group, codon_list, pKa=None, volume=None):
        self.name = name
        self.three_letter = three_letter
        self.one_letter = one_letter
        self.polarity = polarity
        self.charge = charge
        self.r_group = r_group
        self.codon_list = codon_list
        self.pKa = pKa
        self.volume = volume

    def __get_charge__(self, pH=None):
        if pH is not None and self.pKa is not None:
            change = int(math.copysign(1, self.pKa - pH)) if self.pKa != pH else 1
        return self.charge-(change<0)  
    
    def __has_codon__(self, codon):
        return codon.upper() in self.codon_list
    
    def __is_identifier__(self, identifier):
        """
        Check if input matches any form: full name, three-letter code, or one-letter code.
        Case-insensitive.
        """
        identifier = identifier.upper()
        return identifier in {self.name.upper(), self.three_letter.upper(), self.one_letter.upper()}

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
        
