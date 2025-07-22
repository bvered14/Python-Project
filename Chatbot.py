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
    
    print("\nSimulation completed.")
    print(f"Final aggregation state: {tau.aggregation_state}")
    print(f"Phosphorylation sites phosphorylated: {sum(tau.phosphorylation_sites.values())}")
    print(f"History length: {len(tau.history)}")

    # Extract history info for plotting
    ages = [entry['age'] for entry in tau.history]
    phospho_counts = [entry['phospho_count'] for entry in tau.history]
    aggregation_states = [entry['aggregation_state'] for entry in tau.history]

    agg_state_numeric = {'monomer': 0, 'oligomer': 1, 'fibril': 2}
    aggregation_numeric = [agg_state_numeric.get(state, 0) for state in aggregation_states]

    # Plot phosphorylation count and aggregation state over time
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

if __name__ == "__main__":
    main()
