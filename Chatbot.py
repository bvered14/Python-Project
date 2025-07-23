import numpy as np
from TauProtein import TauProtein
import Environment as env
import matplotlib.pyplot as plt

def prompt_float_in_range(prompt, min_val, max_val):
    while True:
        try:
            val = float(input(f"{prompt} [{min_val} - {max_val}]: "))
            if min_val <= val <= max_val:
                return val
            else:
                print(f"Value must be between {min_val} and {max_val}. Try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def prompt_environment():
    print("Please enter the environmental parameters for simulation:")
    temperature = prompt_float_in_range("Temperature (°C)", 30, 42)
    kinase_level = prompt_float_in_range("Kinase level (0.0 - 5.0)", 0.0, 5.0)
    phosphatase_level = prompt_float_in_range("Phosphatase level (0.0 - 5.0)", 0.0, 5.0)
    protease_level = prompt_float_in_range("Protease level (0.0 - 1.0)", 0.0, 1.0)
    oxidative_stress = prompt_float_in_range("Oxidative stress (0.0 - 1.0)", 0.0, 1.0)

    return env.Environment(temperature=temperature,
                           kinase_level=kinase_level,
                           phosphatase_level=phosphatase_level,
                           protease_level=protease_level,
                           oxidative_stress=oxidative_stress)

def prompt_timepoints():
    while True:
        try:
            total_time = int(input("Enter total simulation time (integer time units): "))
            if total_time > 0:
                return np.arange(total_time)
            else:
                print("Time must be positive integer.")
        except ValueError:
            print("Invalid input. Please enter an integer.")

def main():
    print("Tau Protein Phosphorylation & Aggregation Simulator")

    # Create tau protein instance
    tau = TauProtein()

    # Prompt user for environment parameters
    environment = prompt_environment()

    # Prompt user for simulation duration/timepoints
    timepoints = prompt_timepoints()

    # Run update_state simulation
    probabilities = tau.update_state(environment, timepoints)

    print("\nUpdate State Simulation completed.")
    print(f"Final aggregation state: {tau.aggregation_state}")
    print(f"Phosphorylation sites phosphorylated: {sum(tau.phosphorylation_sites.values())}")
    print(f"History length: {len(tau.history)}")

    # Disease simulation step (no change to TauProtein class)
    print("\n--- Disease Simulation ---")
    ptm_options = ["hyperphosphorylation", "truncation", "acetylation"]
    for i, option in enumerate(ptm_options, 1):
        print(f"{i}. {option}")
    while True:
        try:
            ptm_choice = int(input("Choose a PTM (1-3): "))
            if 1 <= ptm_choice <= len(ptm_options):
                ptm = ptm_options[ptm_choice - 1]
                break
            else:
                print("Invalid choice.")
        except ValueError:
            print("Enter a valid number.")

    while True:
        try:
            ptm_time = int(input("Enter number of time steps for PTM simulation: "))
            if ptm_time > 0:
                break
            else:
                print("Time must be a positive integer.")
        except ValueError:
            print("Invalid input.")

    # Simulate phosphorylation increase due to PTM
    phospho_disease = []
    current_phospho = sum(tau.phosphorylation_sites.values())
    for _ in range(ptm_time):
        if ptm == "hyperphosphorylation":
            current_phospho += np.random.randint(2, 5)
        elif ptm == "truncation":
            current_phospho += np.random.randint(1, 3)
        elif ptm == "acetylation":
            current_phospho += np.random.randint(0, 2)
        phospho_disease.append(min(current_phospho, len(tau.phosphorylation_sites)))

    # Extract history info for plotting update_state
    ages = [entry['age'] for entry in tau.history]
    phospho_counts = [entry['phospho_count'] for entry in tau.history]
    aggregation_states = [entry['aggregation_state'] for entry in tau.history]
    agg_state_numeric = {'monomer': 0, 'oligomer': 1, 'fibril': 2}
    aggregation_numeric = [agg_state_numeric.get(state, 0) for state in aggregation_states]

    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,10))

    ax1.set_title("Tau Update State: Phosphorylation & Aggregation")
    ax1.plot(ages, phospho_counts, label="Phosphorylation Count", color='tab:blue')
    ax1.set_ylabel("Phosphorylation")
    ax1b = ax1.twinx()
    ax1b.plot(ages, aggregation_numeric, label="Aggregation", color='tab:red', linestyle='--')
    ax1b.set_ylabel("Aggregation State")
    ax1b.set_yticks([0,1,2])
    ax1b.set_yticklabels(['Monomer', 'Oligomer', 'Fibril'])

    ax2.set_title(f"Disease Simulation: Effect of {ptm}")
    ax2.plot(range(ptm_time), phospho_disease, label=f"{ptm} effect", color='purple')
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Phosphorylation Level")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
