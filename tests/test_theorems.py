"""Frozen values: each test fails if a theorem script's output drifts.  Run the scripts first (./run_all.sh)."""
import json, os, pytest
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'theorems', 'out')
def load(n):
    p = os.path.join(OUT, n + '.json')
    if not os.path.exists(p): pytest.skip(f'{n}.json missing: run theorems/{n}_*.py first')
    return json.load(open(p))
def test_T01_rule():
    d = load('T01')
    assert d['orbit_size'] == 4 and d['reversal_conjugator'] == 'A' and d['det'] == -1
    assert d['M2'] == [[2, 1], [1, 1]] and d['reversal_closed'] is True
    assert d['complexity'] == [2, 3, 4, 5, 6, 7, 8, 9] and d['bb_has_preimage'] is False and d['aaa_preimage'] == ['bbb']
def test_T02_double_tick():
    d = load('T02')
    assert d['detF'] == -1 and d['commute'] is True and d['m000_orientable'] is False and d['cover_is_m004'] is True
    assert all(v == [False, True] for v in d['one_tick_bundles'].values())
    assert [-1, -1, -1, 0] in d['det_minus1_centralisers']
def test_T03_symmetry():
    d = load('T03')
    assert d['sym_order'] == 8 and d['amphichiral'] is True and d['n_isometries'] == 8 and d['n_orientation_reversing'] == 4
    assert d['n_conjugate_realisations'] == 132 and d['first'][0][:2] == ['aabAB', 'abaBAba']
    (x, y, z) = d['character']; assert abs(x[0] + 1.5) < 1e-9 and abs(x[1] - 3**.5/2) < 1e-9 and abs(z[0] + 2) < 1e-9
def test_T04_lattice():
    d = load('T04')
    assert d['cusp_group_order'] == 4 and d['abelian'] is True and d['linear_characters'] == 4
    assert abs(d['CS']) < 1e-9 and abs(d['Vol'] - 2.029883213) < 1e-8 and d['T4_mirror_odd'] and d['T7_time_odd']
def test_T05_lookups():
    d = load('T05')
    assert d['labels'] == {'1': 'A_0', '3': 'E6', '5': 'E8'} and d['imaginary_labelled'] == [3]
    assert d['n_sharing'] == 14 and sorted(d['sharing']) == sorted(['m003','m004','m202','m203','m206','m207','m208','m410','m412','s118','s119','s594','s595','s596'])
    assert d['bianchi_indices']['m004'] == 12.0 and d['bianchi_indices']['m410'] == 30.0 and d['bianchi_indices']['m202'] == 24.0
def test_T06_object():
    d = load('T06')
    assert sorted(d['separators']) == ['H1', 'covers', 'shape'] and d['n_chiral'] == 8
    assert d['rows']['m004']['covers'] == [1, 1, 2, 4, 11] and d['rows']['m003']['covers'] == [1, 1, 2, 8, 7]
    assert d['fillings']['m004']['5'][2] == '0' and d['fillings']['m003']['5'][2] == 'Z/35'
    assert d['fillings']['m004']['5'][1] < 0 < d['fillings']['m003']['5'][1]
    assert all(r < 1e-9 for r in d['curve_residuals'].values()) and d['curves']['m004'] != d['curves']['m003']
def test_T07_dictionary():
    d = load('T07')
    assert d['classes'] == 4 and sorted(d['admissible']) == ['c=C,gamma5=T', 'c=P,gamma5=C', 'c=P,gamma5=T']
    t = d['admissible']['c=P,gamma5=T']; assert t['T4'] == 'W' and t['T7'] == 'K' and t['absent(odd,odd)'] == 'E'
