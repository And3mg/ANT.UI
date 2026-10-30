"""Pull / Grid / Rotation assistants.

The propagation tests stand in for the cluster: they write the Gaussian .log
an optimisation would leave if no atom moved, run the generated link_opt.py,
and require the propagated geometries to equal the ones ANT.UI writes for the
same scan without optimisation.
"""
import os

import numpy as np
import pytest

import _OUT
from antui_testlib import Gjf, fake_gaussian_log, python_syntax_error, run_link_opt

NUMPY2 = int(np.__version__.split('.')[0]) >= 2
CHAINS = ['none', 'first step', 'every step']


def generate(kind, objs, d, name, chains, opt_data, dft_opts, **kw):
    opt = opt_data(Chains=chains, Molecule='fix' if chains == 'none' else 'all')
    if kind == 'pull':
        _OUT.Create_pull_output(objs, d, name, kw['dz'], kw['nz'], kw['both'], opt, dft_opts())
        return kw['nz']
    if kind == 'grid':
        _OUT.Create_grid_output(objs, d, name, kw['grid'], opt, dft_opts())
        return kw['grid'][0][2] * kw['grid'][1][2]
    if kind == 'rot':
        angles = np.arange(kw['frm'], kw['to'] + kw['by'], kw['by'])
        _OUT.Create_rot_output(objs, d, name, angles, kw.get('which', (False, True, False)), False, opt, dft_opts())
        return len(angles)
    raise ValueError(kind)


def step_folders(d, name):
    return sorted(f for f in os.listdir(d) if f.startswith(name + '_') and os.path.isdir(f'{d}/{f}'))


PULL = [dict(dz=0.25, nz=4, both=False), dict(dz=0.25, nz=4, both=True),
        dict(dz=0.1, nz=3, both=False), dict(dz=0.1, nz=3, both=True)]
GRID = [dict(grid=[[-1.0, 1.0, 3], [-0.5, 0.5, 2], 9.0])]
ROT = [dict(frm=10, to=70, by=20), dict(frm=0, to=45, by=15)]


def _case(kind, chains, kw):
    label = '-'.join(f'{k}{v}' for k, v in kw.items() if k != 'grid')
    return pytest.param(kind, chains, kw, id=f'{kind}-{chains.replace(" ", "")}-{label}')


# Regressions guarded here: R01 (pull, first step), R06 (rotation, first step), R07 (grid, first step,
# numpy>=2), R08 (symmetric pull dz=0.1 nz=3).
CASES = ([_case('pull', c, kw) for c in CHAINS for kw in PULL] +
         [_case('grid', c, kw) for c in CHAINS for kw in GRID] +
         [_case('rot', c, kw) for c in CHAINS for kw in ROT])


@pytest.mark.parametrize('kind,chains,kw', CASES)
def test_sequence_is_complete(junction, outdir, opt_data, dft_opts, kind, chains, kw):
    d, name = outdir()
    n = generate(kind, junction(), d, name, chains, opt_data, dft_opts, **kw)
    folders = step_folders(d, name)
    assert len(folders) == n
    assert os.path.exists(f'{d}/launchANT.sh')
    film = open(f'{d}/{name}.xyz').read()
    assert film.count('Timestep') == n
    assert os.path.exists(f'{d}/{folders[0]}/input{folders[0]}.gjf')
    for f in folders:
        assert os.path.exists(f'{d}/{f}/input{f}.ant')
        if chains == 'none':
            assert os.path.exists(f'{d}/{f}/input{f}.gjf')
        if os.path.exists(f'{d}/{f}/link_opt.py'):
            assert python_syntax_error(f'{d}/{f}/link_opt.py') is None


def _prop(kind, chains, kw, *bugs):
    """bugs: the regressions this case guards (see PaperRevision plan / PR history)."""
    label = '-'.join(f'{k}{v}' for k, v in kw.items() if k not in ('grid', 'mol_xy'))
    if 'mol_xy' in kw:
        label += '-offcentre'
    return pytest.param(kind, chains, kw, id=f'{kind}-{chains.replace(" ", "")}-{label}')


# X1: first-step chains moved one extra increment; symmetric pulls halved dz at every step.
# X13: optimised rotation chains turned clockwise.  X14: they rotated about the origin.
# X15: USER_CONF comments broke the route line of the generated inputs.  X16: every-step grids lagged one point.
PROPAGATION = [
    _prop('pull', 'first step', dict(dz=0.5, nz=4, both=False), 'R01', 'X1'),
    _prop('pull', 'first step', dict(dz=0.5, nz=4, both=True), 'R01', 'X1'),
    _prop('pull', 'every step', dict(dz=0.5, nz=4, both=False), 'X15'),
    _prop('pull', 'every step', dict(dz=0.5, nz=4, both=True), 'X15'),
    _prop('grid', 'first step', GRID[0], 'R07'),
    _prop('grid', 'every step', GRID[0], 'X15', 'X16'),
    _prop('rot', 'first step', dict(frm=0, to=45, by=15), 'R06', 'X1', 'X13'),
    _prop('rot', 'every step', dict(frm=0, to=45, by=15), 'X13', 'X15'),
    _prop('rot', 'every step', dict(frm=0, to=45, by=15, mol_xy=(1.0, 0.5)), 'X13', 'X14', 'X15'),
]


@pytest.mark.parametrize('kind,chains,kw', PROPAGATION)
def test_chain_propagates_the_requested_geometries(junction, outdir, opt_data, dft_opts, kind, chains, kw):
    kw = dict(kw)
    mol_xy = kw.pop('mol_xy', None)

    def objs():
        o = junction()
        if mol_xy:
            o[1].spatial_coords = np.array([mol_xy[0], mol_xy[1], 0.0])
        return o

    ref_dir, ref_name = outdir()
    generate(kind, objs(), ref_dir, ref_name, 'none', opt_data, dft_opts, **kw)
    reference = [Gjf(f'{ref_dir}/{f}/input{f}.gjf').coords for f in step_folders(ref_dir, ref_name)]

    d, name = outdir()
    generate(kind, objs(), d, name, chains, opt_data, dft_opts, **kw)
    folders = step_folders(d, name)
    first = f'{d}/{folders[0]}'
    np.testing.assert_allclose(Gjf(f'{first}/input{folders[0]}.gjf').coords, reference[0], atol=1e-4)

    if chains == 'first step':
        fake_gaussian_log(f'{first}/input{folders[0]}.gjf')
        run_link_opt(first)
    else:
        for f in folders[:-1]:
            fake_gaussian_log(f'{d}/{f}/input{f}.gjf')
            run_link_opt(f'{d}/{f}')

    for k, f in enumerate(folders[1:], 1):
        got = Gjf(f'{d}/{f}/input{f}.gjf').coords
        np.testing.assert_allclose(got, reference[k], atol=1e-4, err_msg=f'step {k + 1} ({f})')
