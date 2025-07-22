import pytest
import numpy as np
import matplotlib.pyplot as plt
import tau_simulation
from tau_simulation import env, tau, agg_state_numeric, plot_site_probabilities

@pytest.fixture(autouse=True)
def no_show(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: None)

def test_update_state_populates_history():
    tau.history = []
    timepoints = np.arange(10)
    tau.update_state(env, timepoints)
    assert len(tau.history) == 10
    for entry in tau.history:
        assert 'age' in entry and isinstance(entry['age'], (int, float))
        assert 'phospho_count' in entry and isinstance(entry['phospho_count'], int)
        assert 'aggregation_state' in entry and entry['aggregation_state'] in agg_state_numeric

def test_final_states_consistent_with_history():
    last = tau.history[-1]
    assert tau.aggregation_state == last['aggregation_state']
    total_sites = sum(tau.phosphorylation_sites.values())
    assert total_sites == last['phospho_count']

def test_aggregation_numeric_mapping():
    assert agg_state_numeric == {'monomer': 0, 'oligomer': 1, 'fibril': 2}
    vals = set(agg_state_numeric.values())
    assert vals == {0, 1, 2}

def test_plot_site_probabilities_structure(tmp_path):
    probs = [
        {'S1': 0.1, 'S2': 0.9},
        {'S1': 0.2, 'S2': 0.8},
        {'S1': 0.3, 'S2': 0.7},
    ]
    timepoints = [0,1,2]
    plot_site_probabilities(probs, timepoints)
    fig = plt.gcf()
    ax = plt.gca()
    lines = ax.get_lines()
    assert len(lines) == 2
    for line in lines:
        xdata, ydata = line.get_data()
        assert list(xdata) == timepoints
        assert len(ydata) == 3
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert set(labels) == {'S1', 'S2'}
    plt.close(fig)
