import numpy as np
from typing import Optional
from Protein import Protein
import AA
import Environment
import matplotlib.pyplot as plt

class TauProtein(Protein):
    def __init__(self, name="Tau", isoform="4R", sequence: Optional[np.array]=None,
                 weight=None, length=None, organism=None, location=None, expression_level=None):
        
        super().__init__(name, sequence, weight, length, organism, location, expression_level)
        
        self.isoform = isoform # add isoform defining method!
        self.phosphorylation_sites = {}  # tbd - maybe not necessary here
        self.aggregation_state = "monomer"  # or "oligomer", "fibril"
        self.truncated_site = None
        self.soluble = True

    def define_isoform(self, exon):
        if exon == "R2 at 10":
            self.isoform = "4R"
        if exon == "R1 at 10":
            self.isoform = "3R"
        return self.isoform

    """
    def aggregation(self, ):
        if self.is_level_of_phosphorylation_healthy():
            return KeyError("This tau protein will not aggregate, it is healthy!")
        else:
            agg_score = 0

        phospho_count = len(self.phosphorylated_sites)
        agg_score += phospho_count * 1.5  # Increase weight of phospho

        if self.isoform == "4R":
            agg_score += 2

        if self.is_truncated:
            agg_score += 2

        motif_score = self.count_motifs(['VQIINK', 'VQIVYK'])
        agg_score += motif_score * 2

        if self.temperature > 37:
            agg_score += 0.5

        if self.molecular_weight > self.native_mass + 1000:
            agg_score += 1

        if agg_score >= 7:
            self.aggregation_state = 'fibril'
        elif agg_score >= 4:
            self.aggregation_state = 'oligomer'
        else:
            self.aggregation_state = 'monomer'
    """

    def count_phosphorylated_residues(self):
        return sum(1 for aa in self.sequence if aa.PTM == "Phospho")
    
    def detect_aggregation_motifs(self):
        sequence_str = ''.join(aa.one_letter for aa in self.sequence)
        motifs = ['VQIINK', 'VQIVYK']
        count = 0
        for motif in motifs:
            count += len(re.findall(motif, sequence_str))
        return count
    
    def compute_aggregation_score(self):
        phospho_count = self.count_phosphorylated_residues()
        motif_count = self.detect_aggregation_motifs()
        score = 0

        score += phospho_count * 1.5
        score += motif_count * 2

        if self.is_truncated:
            score += 2

        if self.isoform == '4R':
            score += 1

        return score
    
    def update_aggregation_state(self):
        score = self.compute_aggregation_score()

        if score >= 7:
            self.aggregation_state = 'fibril'
        elif score >= 4:
            self.aggregation_state = 'oligomer'
        else:
            self.aggregation_state = 'monomer'

    
    def is_level_of_phosphorylation_healthy():
        # if phosphorylate.phosphorylate()
        return True
    
    def update_state(self, environment, time_step=1):
        self.apply_phosphorylation(environment)
        self.apply_dephosphorylation(environment)
        self.apply_truncation(environment)
        self.update_aggregation_state()
        self.apply_clearance(environment)
        self.age += time_step

        # Log history
        self.history.append({
            'age': self.age,
            'phospho_count': sum(self.phosphorylation_sites.values()),
            'aggregation_state': self.aggregation_state,
            'is_truncated': self.is_truncated,
            'pathological': self.pathological
        })

env = Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)
tau = TauProtein(isoform='4R')

for t in range(100):
    tau.update_state(env)

# Example: View final state and history
print(f"Final aggregation state: {tau.aggregation_state}")
print(f"Phosphorylation sites phosphorylated: {sum(tau.phosphorylation_sites.values())}")
print("History:", tau.history)

ages = [entry['age'] for entry in tau.history]
phospho_counts = [entry['phospho_count'] for entry in tau.history]
aggregation_states = [entry['aggregation_state'] for entry in tau.history]

agg_state_numeric = {'monomer': 0, 'oligomer': 1, 'fibril': 2}
aggregation_numeric = [agg_state_numeric[state] for state in aggregation_states]

fig, ax1 = plt.subplots(figsize=(10,6))

ax1.set_xlabel('Time Step')
ax1.set_ylabel('Phosphorylation Count', color='tab:blue')
ax1.plot(ages, phospho_counts, color='tab:blue', label='Phosphorylation Count')
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax2 = ax1.twinx()
ax2.set_ylabel('Aggregation State', color='tab:red')
ax2.plot(ages, aggregation_numeric, color='tab:red', linestyle='--', label='Aggregation State')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.set_yticks([0,1,2])
ax2.set_yticklabels(['Monomer', 'Oligomer', 'Fibril'])

plt.title('Tau Phosphorylation & Aggregation Over Time')
fig.tight_layout()
plt.show()

tau = TauProtein()
print(tau)