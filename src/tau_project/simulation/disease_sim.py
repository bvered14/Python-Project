from ..models.aa import AminoAcid
from ..models.protein import Protein
from ..phospho_utils import phosphorylation_constants, phosphorylation, phospo_over_time
from ..plot_utils import plot_tau_summary, plot_site_probabilities, plot_phosphorylation_heatmap  # Shared utilities available
import numpy as np
import matplotlib.pyplot as plt


def run_and_plot_disease_simulation():
    tau_seq = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
    tau_prot = Protein("Tau", tau_seq)
    acetyl_sites = [
        (163, "Acetyl"),
        (174, "Acetyl"),
        (180, "Acetyl"),
    ]

    # phospho_data=phospo_over_time(tau_prot,180, acetyl_sites)
    phospho_data = phospo_over_time(tau_prot, 180)

    plt.figure(figsize=(10, 6))
    plt.plot(phospho_data, label="Phosphorylation %", color="mediumblue", linewidth=2)
    plt.title("Tau Phosphorylation Simulation Over Time")
    plt.xlabel("Time Steps")
    plt.ylabel("Phosphorylated Residues (%)")
    plt.ylim(0, 5)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    print("Close the plot window to return to the menu.")
    plt.show()

# Shared plot utilities (plot_tau_summary, plot_site_probabilities, plot_phosphorylation_heatmap) are available for future use.
