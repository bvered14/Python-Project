"""
phospho_utils.py
Utility functions for phosphorylation and dephosphorylation logic, shared across the project.
"""
import random
import numpy as np

def phosphorylation_constants(protein, list_of_PTMs=None):
    initial_phospo_constant = 0.005
    initial_dephospo_constant = 0.02
    acetyl_constant_ranges = [2, 5, 0.8, 1]
    methyl_constant_ranges = [0.5, 0.9, 1, 1]
    ubi_constant_ranges = [0.5, 1, 1.2, 2]
    glyco_constant_ranges = [0.1, 0.3, 1.5, 2.5]
    phospho_residues = 1
    possible_phopho = 0
    acetyl_counts = 0
    methyl_counts = 0
    ubi_counts = 0
    glyco_counts = 0

    if list_of_PTMs is not None:
        for residue in list_of_PTMs:
            protein.sequence[residue[0] - 1].add_PTM(residue[1])
    for aa in protein.sequence:
        match aa.PTM:
            case "Phospho":
                phospho_residues += 1
                possible_phopho += 1
            case "Acetyl":
                acetyl_counts += 1
            case "Methyl":
                methyl_counts += 1
            case "Ubi":
                ubi_counts += 1
            case "GlcNAc":
                glyco_counts += 1
            case __ if aa.three_letter in {"Ser", "Thr", "Tyr"}:
                possible_phopho += 1

    aK = random.uniform(acetyl_constant_ranges[0], acetyl_constant_ranges[1])
    adK = random.uniform(acetyl_constant_ranges[2], acetyl_constant_ranges[3])
    mK = random.uniform(methyl_constant_ranges[0], methyl_constant_ranges[1])
    mdK = random.uniform(methyl_constant_ranges[2], methyl_constant_ranges[3])
    uK = random.uniform(ubi_constant_ranges[0], ubi_constant_ranges[1])
    udK = random.uniform(ubi_constant_ranges[2], ubi_constant_ranges[3])
    gK = random.uniform(glyco_constant_ranges[0], glyco_constant_ranges[1])
    gdK = random.uniform(glyco_constant_ranges[2], glyco_constant_ranges[3])

    pK = (
        initial_phospo_constant
        * (aK**acetyl_counts)
        * (mK**methyl_counts)
        * (uK**ubi_counts)
        * (gK**glyco_counts)
    )
    dpK = (
        initial_dephospo_constant
        * (adK**acetyl_counts)
        * (mdK**methyl_counts)
        * (udK**ubi_counts)
        * (gdK**glyco_counts)
    )

    return pK, dpK

def phosphorylation(protein, phospho_k, dephospho_k):
    phospho_residues = 0
    possible_phopho = 0
    for aa in protein.sequence:
        match aa.PTM:
            case "Phospho":
                phospho_residues += 1
                possible_phopho += 1
            case __ if aa.three_letter in {"Ser", "Thr", "Tyr"}:
                possible_phopho += 1
    dP = (
        phospho_k * (1 - (phospho_residues / possible_phopho))
        - dephospho_k * phospho_residues
    )

    phospho_residues = 0
    possible_phopho = 0
    for aa in protein.sequence:
        match aa.PTM:
            case "Phospho":
                if random.random() < dephospho_k:
                    aa.remove_PTM()
                    possible_phopho += 1
                else:
                    phospho_residues += 1
                    possible_phopho += 1
            case __ if aa.three_letter in {"Ser", "Thr", "Tyr"}:
                if random.random() < dP:
                    aa.add_PTM("Phosphorylation")
                    phospho_residues += 1
                    possible_phopho += 1
                else:
                    possible_phopho += 1

    return phospho_residues / possible_phopho * 100

def phospo_over_time(protein, time, list_of_PTMs=None):
    p_percentage = np.zeros(time)
    for i in range(time):
        pK, dpK = phosphorylation_constants(protein, list_of_PTMs)
        p_percentage[i] = phosphorylation(protein, pK, dpK)
    return p_percentage 