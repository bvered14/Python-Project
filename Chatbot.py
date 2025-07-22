def get_user_input():
    print("=== Tau Protein Setup ===")
    
    # Isoform input
    isoform = input("Enter tau isoform (3R or 4R): ").strip().upper()
    while isoform not in ("3R", "4R"):
        isoform = input("Invalid input. Please enter '3R' or '4R': ").strip().upper()

    # Phosphorylation sites input
    raw_sites = input("Enter phosphorylation sites (comma-separated, e.g., S202,T205): ")
    phosphorylation_sites = [site.strip().upper() for site in raw_sites.split(",") if site.strip()]
    
    # Binding state input
    bound_input = input("Is the protein bound to microtubules? (yes/no): ").strip().lower()
    while bound_input not in ("yes", "no"):
        bound_input = input("Please enter 'yes' or 'no': ").strip().lower()
    bound_state = (bound_input == "yes")
    
    # Aggregation state (optional or default)
    aggregation_state = "monomer"

    return {
        "isoform": isoform,
        "phosphorylation_sites": phosphorylation_sites,
        "bound_state": bound_state,
        "aggregation_state": aggregation_state
    }

# Example usage
if __name__ == "__main__":
    user_data = get_user_input()
    print("\nYou entered:")
    for k, v in user_data.items():
        print(f"{k}: {v}")
