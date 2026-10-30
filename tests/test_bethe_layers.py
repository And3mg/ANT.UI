"""Bethe layers relaxed along z (R04).

Gaussian09 ties each Bethe plane to its own z variable, so the planes relax
independently. link_opt.py must then keep the innermost plane where Gaussian
put it and move the outer planes back to their original distance from it, so
the block attached to the Bethe lattice keeps its geometry.
"""
import os

import numpy as np
import pytest

import _OUT
from antui_testlib import Gjf, fake_gaussian_log, run_link_opt

BETHE_Z = dict(Top_Bethe='z', Bottom_Bethe='z', Molecule='all')


def bethe_planes(gjf):
    """{variable: atom indices} for the Bethe plane variables (zB*, zT*) of an input."""
    planes = {}
    for k, tokens in enumerate(gjf.symbols):
        if tokens[-1].startswith(('zB', 'zT')):
            planes.setdefault(tokens[-1], []).append(k)
    return planes


def innermost(planes, coords, side):
    pick = max if side == 'zB' else min
    return pick((v for v in planes if v.startswith(side)), key=lambda v: coords[planes[v][0], 2])


def shifts_that_change_spacing(planes):
    """Move every plane by a different amount (innermost included)."""
    return {var: 0.05 * (k + 1) * (-1) ** k for k, var in enumerate(sorted(planes))}


def assert_spacing_restored(original, relaxed_shift, result_coords, planes):
    for side in ('zB', 'zT'):
        ref = innermost(planes, original.coords, side)
        z_ref = original.coords[planes[ref][0], 2] + relaxed_shift[ref]
        for var in (v for v in planes if v.startswith(side)):
            offset = original.coords[planes[var][0], 2] - original.coords[planes[ref][0], 2]
            np.testing.assert_allclose(result_coords[planes[var], 2], z_ref + offset, atol=1e-5,
                                       err_msg=f'plane {var} (reference {ref})')


@pytest.mark.parametrize('electrode', ['Point_FCC_[111]', 'Block_FCC_[111]_80'])  # 2 and 3 Bethe planes per side
def test_optimised_input_keeps_bethe_spacing(junction, outdir, opt_data, dft_opts, electrode):
    d, name = outdir()
    _OUT.Create_output(junction(electrode=electrode), d, name, opt_data(**BETHE_Z), dft_opts())
    gjf = f'{d}/input{name}.gjf'
    original = Gjf(gjf)
    planes = bethe_planes(original)
    assert len([v for v in planes if v.startswith('zB')]) >= 2
    shift = shifts_that_change_spacing(planes)
    fake_gaussian_log(gjf, shift)
    run_link_opt(d)
    assert_spacing_restored(original, shift, Gjf(f'{d}/input{name}_opt.gjf').coords, planes)


@pytest.mark.parametrize('chains', ['first step', 'every step'])
def test_pull_chain_keeps_bethe_spacing(junction, outdir, opt_data, dft_opts, chains):
    d, name = outdir()
    dz, nz = 0.5, 3
    _OUT.Create_pull_output(junction(), d, name, dz, nz, True, opt_data(Chains=chains, **BETHE_Z), dft_opts())
    first = f'{d}/{name}_1/input{name}_1.gjf'
    original = Gjf(first)
    planes = bethe_planes(original)
    shift = shifts_that_change_spacing(planes)
    fake_gaussian_log(first, shift)
    run_link_opt(f'{d}/{name}_1')

    step2 = Gjf(f'{d}/{name}_2/input{name}_2.gjf').coords
    # step 2 = restored relaxed geometry, then the top moves +dz/2 and the bottom -dz/2
    expected_shift = dict(shift)
    for var in planes:
        expected_shift[var] += dz / 2 if var.startswith('zT') else -dz / 2
    assert_spacing_restored(original, expected_shift, step2, planes)
