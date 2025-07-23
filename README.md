# Tau Protein Simulation Project

## Overview
This project simulates tau protein phosphorylation, aggregation, and disease-related modifications using modern Python best practices. It is designed for research, education, and demonstration of scientific Python development.

---

## Features
- Object-oriented design for proteins, amino acids, and environment
- Modular, PEP8-compliant codebase (autoformatted with `black`)
- Integration of the scientific Python stack: numpy, matplotlib, seaborn, pandas
- Interactive CLI chatbot for simulation and visualization
- Comprehensive pytest-based test suite
- Ready for packaging and PyPI distribution

---

## Project Structure
```
Tau-Project/
├── src/
│   └── tau_project/
│       ├── models/
│       │   ├── aa.py
│       │   ├── protein.py
│       │   ├── tau_protein.py
│       │   └── truncation.py
│       ├── simulation/
│       │   ├── tau_simulation.py
│       │   └── disease_sim.py
│       ├── environment.py
│       ├── phosphorylation.py
│       └── chatbot.py
├── tests/
│   ├── test_tau_simulation.py
│   └── test_disease_sim.py
├── requirements.txt
├── setup.py
├── .gitignore
└── README.md
```

---

## Installation
1. Clone the repository:
   ```sh
   git clone <your-repo-url>
   cd Tau-Project
   ```
2. Install dependencies:
   ```sh
   pip install -r requirements.txt
   ```
3. (Optional) Install in development mode:
   ```sh
   pip install -e .
   ```

---

## Usage
- **Run the interactive chatbot:**
  ```sh
  python -m src.tau_project.chatbot
  ```
- **Run a simulation directly:**
  ```sh
  python -m src.tau_project.simulation.tau_simulation
  ```
- **Run tests:**
  ```sh
  pytest
  ```

---

## Scientific Stack
- [numpy](https://numpy.org/): numerical computation
- [matplotlib](https://matplotlib.org/): plotting
- [seaborn](https://seaborn.pydata.org/): advanced visualization
- [pandas](https://pandas.pydata.org/): data handling
- [pytest](https://docs.pytest.org/): testing

---

## Design Patterns & Coding Paradigms
- **Object-Oriented Programming (OOP):**
  - `AminoAcid`, `Protein`, `TauProtein`, and `Environment` are all classes.
- **Factory Pattern:**
  - Used in `Protein` to create amino acid objects from data.
- **Strategy Pattern:**
  - Simulation logic is modular and can be extended for new behaviors.
- **Static Methods:**
  - Used for sequence truncation and utility functions.

---

## Testing
- All code is tested with `pytest`.
- To run all tests:
  ```sh
  pytest
  ```
- Tests are located in the `tests/` directory.

---

## Contributing
- Use git for version control.
- Follow PEP8 and autoformat with `black`.
- Add/expand tests for new features.
- Document all public classes and functions with docstrings.

---

## Packaging & PyPI
- Install in development mode with `pip install -e .`
- To build for PyPI:
  ```sh
  python setup.py sdist bdist_wheel
  twine upload dist/*
  ```
- Update `setup.py` with your name, email, and project details before publishing.

---

## License
MIT License

---

## Acknowledgments
- Developed as a demonstration of best practices in scientific Python development.
- Integrates modern design, testing, and documentation standards.
