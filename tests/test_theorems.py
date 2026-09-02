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
    assert d['same_type_ratio_is_even'] is True and d['even_in_image'] is False and d['predictions'] == []
def test_T08_spectrum():
    d = load('T08')
    assert d['sites'] == 10947 and d['max_label_residual'] < 3 / d['sites']
    assert all(v[0] for v in d['k4_same'].values()) and d['four_letter_max_residual'] > 1e-3
    assert d['trace_map_matches_fricke'] is True and d['fricke_spread'] < 1e-12
    assert sorted(abs(r[3]) for r in d['labels']) == [1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 6, 6]
def test_T08_spectrum():
    d = load('T08')
    assert d['sites'] == 10946 and d['max_label_residual'] < 3 / d['sites'] and d['four_letter_max_residual'] > 1e-3
    assert all(v[0] for v in d['k4_same'].values()) and d['trace_map_matches_fricke'] is True
    assert abs(d['fricke_invariant'] - (4 + 4 * d['V'] ** 2)) < 1e-9 and d['fricke_spread'] == 0.0
def test_T09_object_point():
    d = load('T09')
    assert d['fixed_curve'] == {'y': 'x/(x - 1)', 'z': 'x/(x - 1)'} and d['minpoly'] == 'x**2 - 3*x + 3' and d['discriminant'] == -3
    assert abs(d['tr_comm'][0] + 2) < 1e-9 and abs(d['tr_comm'][1]) < 1e-9
    assert d['conj_residual'] < 1e-9 and abs(abs(d['trT'][0]) - 2) < 1e-6 and abs(d['trT'][1]) < 1e-6 and d['cusp_abelian_residual'] < 1e-9
    m1 = d['metallic_fixed_points']['1']; assert '3 - 3*x^1 + 1*x^2' in m1 and '3 + 3*x^1 + 1*x^2' in m1
    assert all(len(d['metallic_fixed_points'][str(m)]) >= 1 for m in range(2, 5))
def test_T10_level_sets():
    d = load('T10')
    assert d['chain_kappa'] == '4*(V**2 + 1)' and d['n_fillings'] > 50 and d['real_longitude_fillings'] == []
    assert abs(d['longitude_trace_complete'][0] + 2) < 1e-8

def test_T13_metallic_fields():
    d = load('T13')
    assert all(d['geometric_found'][str(m)] for m in range(1, 5))
    assert d['quadratic_discriminants']['1'] == [-3, -3, -3, -3] or -3 in d['quadratic_discriminants']['1']
    assert all(-3 not in d['quadratic_discriminants'][str(m)] for m in (2, 3, 4))
def test_T13b_bundles():
    d = load('T13b'); c = load('T13b_closing')
    assert d['1']['is_m004'] is True and [3, -3, 1] in d['1']['trace_minpolys'].values()
    assert [8, 0, -4, 0, 1] in d['2']['trace_minpolys'].values() and [8, -4, 1] in d['2']['square_minpolys'].values()
    assert d['3']['trace_degrees'] == [1, 8] and d['4']['square_degrees'] == [1, 4]
    assert all(c[str(m)] is True for m in range(1, 5))
def test_T14_gap_law():
    d = load('T14')
    slopes = [v[0] for v in d['fits'].values()]
    assert all(0.8 < s < 1.15 for s in slopes) and d['sites'] == 10946
    assert d['widths']['1'][-1] > d['widths']['1'][0] * 10
def test_T15_metallic_chains():
    d = load('T15')
    for m in ('1', '2', '3'):
        assert d[m]['max_residual'] < 3 / d[m]['sites'] and abs(d[m]['omega'] - d[m]['freq_b']) < 1e-6
    assert d['2']['max_residual_wrong_omega'] > 1e-2 and d['3']['max_residual_wrong_omega'] > 1e-2

def test_T17_beta_odd_bit():
    d = load('T17')
    assert d['sym_order'] == 8 and abs(d['cs_m004']) < 1e-9
    assert d['types'] == {'T4 (chirality row)': 'W', 'T7 = T3 (time row)': 'K', 'CS (absent axis)': 'E'}
    assert d['chiral_census'][0] == 594 and d['chiral_census'][1] == 1
    assert 'theta-bar' in d['prediction'] and 'm004 -> 0' in d['prediction']
    assert d['field_class']['m004']['cs'] == 0.0 and d['field_class']['m003']['cs'] == 0.25
