"""Single-junction exports: Create ANT, SOC ANT, Create xyz."""
import os

import numpy as np
import pytest

import _OUT
from antui_testlib import Gjf, basis_ecp_problems, electron_count, python_syntax_error, spin_consistent

METALS = ['Au', 'Ag', 'Cu', 'Al', 'Ni', 'Pt']


def _ant_case(metal, drop):
    # X2: Pt and Ni have an even electron count per atom with the shipped ECPs
    return pytest.param(metal, drop, id=f'{metal}-{40 - drop}atoms')


@pytest.mark.parametrize('metal,drop', [_ant_case(m, d) for m in METALS for d in (0, 1)])
def test_ant_input_is_complete(junction, outdir, opt_data, dft_opts, metal, drop):
    objs = junction(top=metal, bottom=metal, drop_top=drop)
    n_atoms = sum(o.N for o in objs)
    d, name = outdir()
    _OUT.Create_output(objs, d, name, opt_data(), dft_opts())
    g = Gjf(f'{d}/input{name}.gjf')
    assert len(g.elements) == n_atoms
    assert basis_ecp_problems(g) == []
    assert spin_consistent(g), f'{electron_count(g)} electrons, charge {g.charge}, multiplicity {g.mult}'
    for f in ('.ant', '.xyz'):
        assert os.path.exists(f'{d}/input{name}{f}' if f == '.ant' else f'{d}/{name}{f}')
    assert os.path.exists(f'{d}/launchANT.sh')


def _soc_case(bottom, top):
    # R03: mixed electrodes lost the active-metal SOC ECP when the last metal had none (Cu, Al)
    return pytest.param(bottom, top, id=f'{bottom}-{top}')


@pytest.mark.parametrize('bottom,top', [_soc_case(*p) for p in
                                        [('Au', 'Au'), ('Pt', 'Pt'), ('Ag', 'Ag'), ('Cu', 'Cu'), ('Al', 'Al'),
                                         ('Au', 'Cu'), ('Cu', 'Au'), ('Ag', 'Cu'), ('Au', 'Pt')]])
def test_soc_input_is_complete(junction, outdir, opt_data, dft_opts, bottom, top):
    objs = junction(top=top, bottom=bottom)
    n_atoms = sum(o.N for o in objs)
    d, name = outdir()
    _OUT.Create_output_soc(objs, d, name, opt_data(), dft_opts())
    g = Gjf(f'{d}/input{name}.gjf')
    assert len(g.elements) == n_atoms
    assert basis_ecp_problems(g) == []
    assert spin_consistent(g)


def test_soc_factors_reach_every_active_electrode_atom(junction, outdir, opt_data, dft_opts):
    """The 2026 recalibration: each active atom gets its own metal's factors, the Bethe region none."""
    objs = junction(bottom='Au', top='Pt')
    d, name = outdir()
    _OUT.Create_output_soc(objs, d, name, opt_data(), dft_opts())
    lines = open(f'{d}/input{name}.ant').read().split('SOCFACATOM')[1].split('\n')
    n = int(lines[1])
    factors = [l.split(',', 1)[1].split() for l in lines[2:2 + n]]
    assert n == len(factors) and n > 0
    assert {' '.join(f) for f in factors} == {'201.3, 43.2, 12.6, 2.0', '229.4, 41.6, 9.3, 2.0'}


def test_soc_with_optimisation_writes_valid_connector(junction, outdir, opt_data, dft_opts):
    d, name = outdir()
    _OUT.Create_output_soc(junction(), d, name, opt_data(Molecule='all'), dft_opts())
    assert python_syntax_error(f'{d}/link_opt.py') is None


def test_bethe_layers_move_rigidly(junction, outdir, opt_data, dft_opts):
    d, name = outdir()
    _OUT.Create_output(junction(), d, name, opt_data(Top_Bethe='z', Bottom_Bethe='z'), dft_opts())
    g = Gjf(f'{d}/input{name}.gjf')
    # Point_FCC_[111] has two Bethe layers per electrode -> one shared z per layer.
    assert len(g.variables) == 4, sorted(g.variables)


@pytest.fixture
def bromobenzene():
    path = 'Molecules/zzpytest_bromobenzene.xyz'
    with open(path, 'w') as f:
        f.write('12\n\nC 0.000 1.395 0\nC 1.208 0.698 0\nC 1.208 -0.698 0\nC 0.000 -1.395 0\n'
                'C -1.208 -0.698 0\nC -1.208 0.698 0\nH 0 2.48 0\nH 2.148 1.24 0\nH 2.148 -1.24 0\n'
                'H 0 -2.48 0\nH -2.148 -1.24 0\nBr -2.85 1.65 0\n')
    yield 'zzpytest_bromobenzene'
    os.remove(path)


def test_every_molecule_element_is_written(junction, outdir, opt_data, dft_opts, bromobenzene):
    objs = junction(mol=bromobenzene)
    d, name = outdir()
    _OUT.Create_output(objs, d, name, opt_data(), dft_opts())
    g = Gjf(f'{d}/input{name}.gjf')
    assert len(g.elements) == sum(o.N for o in objs)
    assert basis_ecp_problems(g) == []


def test_absolute_energy_scripts_compile(junction, outdir, opt_data, dft_opts):
    d, name = outdir()
    _OUT.Create_output(junction(), d, name, opt_data(Molecule='all'), dft_opts(EABS=True))
    for script in ('link_opt.py', 'Get_EABS.py'):
        assert python_syntax_error(f'{d}/{script}') is None


def test_settings_do_not_leak_between_exports(junction, outdir, opt_data, dft_opts):
    d1, n1 = outdir()
    _OUT.Create_output(junction(), d1, n1, opt_data(Molecule='all'), dft_opts(EABS=True))
    d2, n2 = outdir()
    _OUT.Create_output_soc(junction(), d2, n2, opt_data(), dft_opts(EABS=False))
    assert 'EABS' not in open(f'{d2}/launchANT.sh').read()


def test_gaussian_header_has_no_config_comments(junction, outdir, opt_data, dft_opts):
    d, name = outdir()
    _OUT.Create_output(junction(), d, name, opt_data(Molecule='all'), dft_opts())
    lines = open(f'{d}/input{name}.gjf').read().split('\n')
    route = next(i for i, l in enumerate(lines) if l.lstrip().startswith('#'))
    assert all('#' not in l for l in lines[:route]), lines[:route]


def test_xyz_export_has_every_atom(junction, outdir):
    objs = junction()
    d, name = outdir()
    _OUT.Create_xyz_output(objs, d, name)
    lines = open(f'{d}/{name}.xyz').read().split('\n')
    assert int(lines[0]) == sum(o.N for o in objs)
    assert len([l for l in lines[2:] if l.strip()]) == int(lines[0])
