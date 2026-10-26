"""Pre-export summary and checks (Referee 2 suggestion 2 / report items 9-10)."""
import os

import numpy as np
import pytest

import _DICS
import _OUT
from _SUMMARY import describe_export, summary_text
from antui_testlib import Gjf


def text_of(lines, warnings):
    return '\n'.join(lines), '\n'.join(warnings)


def test_summary_matches_the_written_input(junction, outdir, opt_data, dft_opts):
    objs = junction()
    lines, warnings = describe_export(objs, 'ANT', opt_data(), dft_opts())
    d, name = outdir()
    _OUT.Create_output(junction(), d, name, opt_data(), dft_opts())
    g = Gjf(f'{d}/input{name}.gjf')
    text, _ = text_of(lines, warnings)
    assert f'Atoms: {len(g.elements)} = bottom electrode 20' in text
    assert f'charge {g.charge}, multiplicity {g.mult} (automatic)' in text
    assert 'molecule 12 (Benzene, C6H6)' in text
    assert 'C, H: lanl2dz / lanl2dz' in text
    assert 'Au electrode (' in text and 'Au Bethe anchor (' in text


def test_shipped_config_is_flagged(junction, opt_data, dft_opts):
    _, warnings = text_of(*describe_export(junction(), 'ANT', opt_data(), dft_opts()))
    assert 'JOBHEAD' in warnings


def test_mixed_soc_lists_every_ecp(junction, opt_data, dft_opts):
    text, warnings = text_of(*describe_export(junction(top='Cu'), 'SOC', opt_data(), dft_opts()))
    assert 'Au electrode' in text and 'SOC_params/Au_SOC_pseudopot.gbs (60 core electrons)' in text
    assert 'Cu electrode' in text and 'none (all-electron basis)' in text
    assert 'Au SOC parameters' in text and 'Cu SOC parameters' in text


def test_overlapping_parts_are_flagged(junction, opt_data, dft_opts):
    _, warnings = text_of(*describe_export(junction(top_z=2.0), 'ANT', opt_data(), dft_opts()))
    assert 'overlap' in warnings


def test_compression_overlap_is_flagged_at_the_last_step(junction, opt_data, dft_opts):
    params = dict(dz=-2.0, nz=6, both=True)    # each electrode moves 5 A towards the molecule
    text, warnings = text_of(*describe_export(junction(), 'pull', opt_data(), dft_opts(), params))
    assert 'total change -10.000 A' in text
    assert 'Last step contacts' in warnings


def test_impossible_user_spin_is_flagged(junction, opt_data, dft_opts):
    _, warnings = text_of(*describe_export(junction(), 'ANT', opt_data(), dft_opts(Multiplicity=2)))
    assert 'impossible' in warnings


def test_user_charge_reaches_the_input(junction, outdir, opt_data, dft_opts):
    d, name = outdir()
    _OUT.Create_output(junction(), d, name, opt_data(), dft_opts(Shell='Open Shell', Charge=-1, Multiplicity=2))
    g = Gjf(f'{d}/input{name}.gjf')
    assert (g.charge, g.mult) == (-1, 2)


def test_chain_without_optimisation_is_flagged(junction, opt_data, dft_opts):
    params = dict(dz=0.25, nz=4, both=False)
    _, warnings = text_of(*describe_export(junction(), 'pull', opt_data(Chains='first step'), dft_opts(), params))
    assert 'nothing will be optimised' in warnings


def test_rotation_summary_lists_the_angles(junction, opt_data, dft_opts):
    params = dict(angles=np.array([10., 30., 50., 70.]), components=(False, True, False), soc=False)
    text, _ = text_of(*describe_export(junction(), 'rotation', opt_data(), dft_opts(), params))
    assert '4 configurations at 10, 30, 50, 70 degrees' in text and 'rotating: molecule' in text


def test_summary_text_layout():
    out = summary_text(['a', 'b'], [], 'Outputs/X-0001', 'v1.1.0')
    assert out.startswith('ANT.UI v1.1.0 export summary\nDestination: Outputs/X-0001')
    assert out.rstrip().endswith('Warnings: none')
    assert summary_text(['a'], ['w1'], 'X').rstrip().endswith('Warnings:\n  ! w1')


def test_site_config_reload(tmp_path):
    conf = tmp_path / 'USER_CONF.txt'
    conf.write_text('#%%\nJOBHEAD\n#!/bin/bash\n#SBATCH -p jjp\n#%%\n\n#%%\nG09call\ng16 # my call\n#%%\n'
                    '\n#%%\nANTDIRcall\n/opt/ant/bin\n#%%\n\n#%%\nOUTPUTDIR\nD:/calcs/ant runs  # spaces ok\n#%%\n')
    try:
        _DICS.load_user_conf(str(conf))
        assert _DICS.UserInfo['G09call'] == 'g16'
        assert _DICS.UserInfo['JOBHEAD'] == '#!/bin/bash\n#SBATCH -p jjp\n'
        assert _DICS.UserInfo['OUTPUTDIR'] == 'D:/calcs/ant runs'
        assert _DICS.UserInfo['ZMAX'] == 20.0          # default kept
        assert _DICS.site_config_warnings() == []
    finally:
        _DICS.load_user_conf()
    assert _DICS.site_config_warnings()                 # shipped placeholder again
