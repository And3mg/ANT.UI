"""Get U: gas-phase N, N-1, N+1 inputs and the Get_U.py analysis (Referee 1, point 5)."""
import os
import subprocess
import sys

import numpy as np
import pytest

import _DFT
import _OUT
from _DICS import parse_u_line
from _SUMMARY import describe_export
from antui_testlib import Gjf, basis_ecp_problems, electron_count, python_syntax_error, spin_consistent

HARTREE = 27.211386
# PBE/lanl2dz glycerol, the gas-phase runs behind Molecules/Glycerol.xyz (cluster, Glycerol-0002/DFTU)
GLYCEROL = {'N': (-344.344679117, [-0.23463, -0.21791], [0.00960, 0.05388], None),
            'cation': (-344.009018908, [-0.49707, -0.47955], [-0.19510, -0.15124], ([-0.48563, -0.46911], [-0.45252, -0.19195])),
            'anion': (-344.236564229, [0.20526], [0.24378, 0.25443], ([-0.04003, -0.01720], [0.22723, 0.25926]))}


def file_atoms(mol):
    lines = open(f'Molecules/{mol}.xyz').read().split('\n')
    return [l.split() for l in lines[2:2 + int(lines[0])]]


def fake_single_point_log(path, energy, occ, virt, beta = None, normal = True):
    """The parts of a Gaussian single-point log that Get_U.py reads."""
    out = [f' SCF Done:  E(RPBE-PBE) =  {energy:.9f}     A.U. after   20 cycles',
           ' Population analysis using the SCF Density.',
           ' Alpha  occ. eigenvalues -- ' + ''.join(f'{x:10.5f}' for x in occ),
           ' Alpha virt. eigenvalues -- ' + ''.join(f'{x:10.5f}' for x in virt)]
    if beta:
        out += ['  Beta  occ. eigenvalues -- ' + ''.join(f'{x:10.5f}' for x in beta[0]),
                '  Beta virt. eigenvalues -- ' + ''.join(f'{x:10.5f}' for x in beta[1])]
    if normal:
        out.append(' Normal termination of Gaussian 09 at Sun Aug 18 18:09:02 2026.')
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')


def run_get_u(folder, mol, data = GLYCEROL, broken = ()):
    for tag, (energy, occ, virt, beta) in data.items():
        fake_single_point_log(f'{folder}/{mol}_{tag}.log', energy, occ, virt, beta, normal = tag not in broken)
    return subprocess.run([sys.executable, 'Get_U.py'], cwd = folder, capture_output = True, text = True)


@pytest.mark.parametrize('mol', ['Benzene', 'Glycerol', 'Toluene'])
def test_three_charge_states_are_written(junction, outdir, dft_opts, mol):
    objs = junction(mol = mol)
    d, name = outdir()
    _OUT.Create_getU_output(objs, d, name, dft_opts(Functional = 'PBE'))
    charges = {}
    for tag in ('N', 'cation', 'anion'):
        g = Gjf(f'{d}/{mol}_{tag}.gjf')
        charges[tag] = g.charge
        assert basis_ecp_problems(g) == []
        assert spin_consistent(g), f'{tag}: {electron_count(g)} electrons, charge {g.charge}, multiplicity {g.mult}'
        assert ('rPBEPBE' if g.mult == 1 else 'uPBEPBE') in g.route and 'nosymm' in g.route
        assert '%Subst' not in g.text                      # plain Gaussian, no ANT
        assert [[e] + list(map(float, r)) for e, r in zip(g.elements, g.coords)] == \
               [[a[0]] + list(map(float, a[1:4])) for a in file_atoms(mol)]
    assert charges == {'N': 0, 'cation': 1, 'anion': -1}
    launch = open(f'{d}/launchGaussian.sh').read()
    assert all(f'{mol}_{tag}.gjf' in launch for tag in charges) and launch.rstrip().endswith('Get_U.py')
    assert python_syntax_error(f'{d}/Get_U.py') is None
    assert _OUT.verify_output(d) == 0


def test_odd_molecule_gets_a_doublet_and_singlet_ions(dft_opts):
    states = {tag: (charge, mult, shell) for tag, charge, mult, shell in _DFT.gas_states(41, dft_opts())}
    assert states == {'N': (0, 2, 'u'), 'cation': (1, 1, 'r'), 'anion': (-1, 1, 'r')}


def test_u_from_the_glycerol_runs(junction, outdir, dft_opts):
    """U = (IP - EA) - KS gap, the correction ANT splits as -U/2 / +U/2. HSE06 is added; the PBE value stays."""
    d, name = outdir()
    _OUT.Create_getU_output(junction(mol = 'Glycerol'), d, name, dft_opts(Functional = 'HSE06'))
    r = run_get_u(d, 'Glycerol')
    assert r.returncode == 0, r.stderr
    E = {tag: v[0] for tag, v in GLYCEROL.items()}
    IP, EA = (E['cation'] - E['N'])*HARTREE, (E['N'] - E['anion'])*HARTREE
    gap = (0.00960 + 0.21791)*HARTREE
    U = IP - EA - gap
    assert U == pytest.approx(5.885, abs = 1e-3)
    report = open(f'{d}/U.Glycerol.dat').read()
    assert f'U = Sigma_occ + Sigma_unocc = {U:.3f} eV' in report
    assert f'IP - EA = {IP-EA:.3f} eV' in report
    lines = open(f'{d}/Glycerol.xyz').read().split('\n')
    assert parse_u_line(lines[1]) == {'PBE': '11.86', 'HSE06': f'{U:.2f}'}
    assert [l.split() for l in lines[2:-1]] == file_atoms('Glycerol')


def test_new_u_replaces_the_old_one_for_the_same_functional(junction, outdir, dft_opts):
    d, name = outdir()
    _OUT.Create_getU_output(junction(mol = 'Glycerol'), d, name, dft_opts(Functional = 'PBE'))
    assert run_get_u(d, 'Glycerol').returncode == 0
    assert parse_u_line(open(f'{d}/Glycerol.xyz').read().split('\n')[1]) == {'PBE': '5.88'}


def test_failed_gaussian_run_gives_no_u(junction, outdir, dft_opts):
    d, name = outdir()
    _OUT.Create_getU_output(junction(mol = 'Glycerol'), d, name, dft_opts(Functional = 'PBE'))
    r = run_get_u(d, 'Glycerol', broken = ('anion',))
    assert r.returncode != 0 and 'Glycerol_anion.log' in r.stderr
    assert not os.path.exists(f'{d}/Glycerol.xyz')


@pytest.mark.parametrize('comment,expected', [('\tEnergy:      10.5736998', {}),
                                              ('PBE 11.86', {'PBE': '11.86'}),
                                              ('PBE 11.86|HSE06 4.1\n', {'PBE': '11.86', 'HSE06': '4.1'}),
                                              ('PBE abc|XC 3.0', {}), ('', {})])
def test_comment_line_keeps_only_functional_u_pairs(comment, expected):
    assert parse_u_line(comment) == expected


def test_summary_lists_the_three_calculations(junction, dft_opts):
    lines, warnings = describe_export(junction(mol = 'Glycerol'), 'getU', DFT_opts = dft_opts(Functional = 'PBE'))
    text = '\n'.join(lines)
    assert 'Glycerol_N.gjf: neutral, N electrons; charge 0, multiplicity 1 (rPBEPBE)' in text
    assert 'Glycerol_cation.gjf: cation, N-1; charge 1, multiplicity 2 (uPBEPBE)' in text
    assert 'Glycerol_anion.gjf: anion, N+1; charge -1, multiplicity 2 (uPBEPBE)' in text
    assert 'already has U = 11.86 eV for PBE' in text
    assert 'U_ELECTRON_SCREENING' in text
