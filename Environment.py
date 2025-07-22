"""
Environment.py
Defines the Environment class, which models the biological environment for tau protein simulation (temperature, kinase, phosphatase, protease, oxidative stress).
"""
# Environment class models the simulation environment for tau protein
class Environment:
    def __init__(self, temperature=37, kinase_level=1.0, phosphatase_level=1.0, protease_level=1.0, oxidative_stress=0.0):
        # Initialize environment properties
        self.temperature = temperature
        self.kinase_level = kinase_level
        self.phosphatase_level = phosphatase_level
        self.protease_level = protease_level
        self.oxidative_stress = oxidative_stress 