import pytest
import numpy as np
from src.tau_project.simulation.disease_sim import (
    phosphorylation_constants,
    phosphorylation,
    phospo_over_time,
)


class DummyAA:
    def __init__(self, three_letter, PTM=None):
        self.three_letter = three_letter
        self.PTM = PTM or ""

    def add_PTM(self, ptm):
        self.PTM = ptm

    def remove_PTM(self):
        self.PTM = ""


class DummyProteinClass:
    def __init__(self, sequence):
        self.sequence = [DummyAA(aa) for aa in sequence.split("-")]


@pytest.fixture(autouse=True)
def patch_random(monkeypatch):
    monkeypatch.setattr("random.uniform", lambda a, b: (a + b) / 2)
    monkeypatch.setattr("random.random", lambda: 1)
    yield


def test_phosphorylation_constants_no_ptms():
    prot = DummyProteinClass("Ser-Thr-Gly-Tyr-Ala")
    pK, dpK = phosphorylation_constants(prot, list_of_PTMs=None)
    assert pytest.approx(pK, rel=1e-6) == 0.005
    assert pytest.approx(dpK, rel=1e-6) == 0.02


def test_phosphorylation_constants_with_ptms():
    prot = DummyProteinClass("Ser-Ala-Thr-Gly")
    ptms = [(2, "Acetyl"), (4, "Methyl")]
    pK, dpK = phosphorylation_constants(prot, list_of_PTMs=ptms)
    assert prot.sequence[1].PTM == "Acetyl"
    assert prot.sequence[3].PTM == "Methyl"
    expected_pK = 0.005 * (3.5**1) * (0.7**1)
    expected_dpK = 0.02 * ((0.8 + 1) / 2) ** 1 * ((1 + 1) / 2) ** 1
    assert pytest.approx(pK, rel=1e-6) == expected_pK
    assert pytest.approx(dpK, rel=1e-6) == expected_dpK


def test_phosphorylation_zero_rates():
    prot = DummyProteinClass("Ser-Ser-Tyr")
    result = phosphorylation(prot, phospho_k=0, dephospho_k=0)
    assert result == 0.0


def test_phospo_over_time_length_and_values(monkeypatch):
    monkeypatch.setattr(
        "src.tau_project.phospho_utils.phosphorylation_constants",
        lambda protein, list_of_PTMs=None: (0.1, 0.05),
    )
    monkeypatch.setattr(
        "src.tau_project.phospho_utils.phosphorylation",
        lambda protein, pk, dpk: 42.0,
    )
    prot = DummyProteinClass("Ser-Ser")
    time_steps = 10
    data = phospo_over_time(prot, time_steps)
    assert isinstance(data, np.ndarray)
    assert len(data) == time_steps
    assert all(val == 42.0 for val in data)
