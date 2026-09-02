"""Frozen checks for the audit scripts A7 and A8 (the record's downstream results recomputed on this tree's e6)."""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.join(os.path.dirname(HERE), 'audit')

def run(script):
    p = subprocess.run([sys.executable, os.path.join(AUDIT, script)], capture_output=True, text=True, timeout=600)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout

def test_a7_hypercharge_weight_level():
    out = run('A7_hypercharge_at_trinification_landing.py')
    assert 'A7 PASS' in out
    assert 'carrying the SM hypercharge multiset: 18' in out
    assert 'pure on one factor among them: 0' in out
    assert 'pure on one factor: 9 (target has 2)' in out

def test_a8_hatch_centralizers():
    if not os.path.exists(os.path.join(AUDIT, 'e6_data.npz')):
        run('A1_build_e6.py')
    out = run('A8_trinification_hatch.py')
    assert 'A8 PASS' in out
    assert 'centralizer dim 16 (expected 16), derived dim 16, rank 4' in out
    assert 'centralizer dim 35 (expected 35), derived dim 35, rank 5' in out
