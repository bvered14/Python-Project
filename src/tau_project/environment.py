"""
environment.py
Defines the Environment class, which models the biological environment for tau protein simulation (temperature, kinase, phosphatase, protease, oxidative stress).
"""

class Environment:
    """
    Models the simulation environment for tau protein.
    Encapsulates biological parameters such as temperature, kinase, phosphatase, protease, and oxidative stress.
    """
    def __init__(self, temperature=37, kinase_level=1.0, phosphatase_level=1.0, protease_level=1.0, oxidative_stress=0.0):
        """
        Initialize an Environment instance.
        Args:
            temperature (float): Temperature in Celsius.
            kinase_level (float): Kinase activity level.
            phosphatase_level (float): Phosphatase activity level.
            protease_level (float): Protease activity level.
            oxidative_stress (float): Oxidative stress level.
        """
        self.temperature = temperature
        self.kinase_level = kinase_level
        self.phosphatase_level = phosphatase_level
        self.protease_level = protease_level
        self.oxidative_stress = oxidative_stress
