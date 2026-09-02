import json, os, numpy as np, pytest
HERE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'audit')
def test_a2_invariants():
    p = os.path.join(HERE, 'a2_data.npz')
    if not os.path.exists(p): pytest.skip('run audit/A1 and A2 first')
    d = np.load(p); assert sorted(int(x) for x in d['degrees']) == [8, 14, 16, 22] and d['inv'].shape == (4, 78)
def test_a6_stratification():
    p = os.path.join(HERE, 'a6_data.npz')
    if not os.path.exists(p): pytest.skip('run audit/A6 first')
    d = np.load(p); counts = sorted(int(c) for c in d['counts'])
    assert len(counts) == 15 and counts == [2] * 6 + [6] * 9
