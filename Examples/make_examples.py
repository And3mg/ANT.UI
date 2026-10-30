"""Regenerate the worked examples of the user manual without the GUI.

    python Examples/make_examples.py

Each example builds the junction exactly as the GUI does (same electrode
templates, molecule files and placement), writes it with the same generators
as the buttons and assistants, and saves the export summary next to it.
Existing example folders are replaced.
"""
import importlib.util
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
spec = importlib.util.spec_from_file_location('antui', os.path.join(ROOT, 'ANT.UI.py'))
antui = importlib.util.module_from_spec(spec)
spec.loader.exec_module(antui)      # also moves the working directory to the program folder
import _OUT
from _SUMMARY import describe_export, summary_text

OPT_NONE = {'Top_Bethe': 'fix', 'Top': 'fix', 'Bottom_Bethe': 'fix', 'Bottom': 'fix', 'Molecule': 'fix',
            'Molecule_z_atoms': [], 'Chains': 'none', 'Opt_steps': 20, 'Opt_in_ANT': False}


def junction(molecule, mol_z, top_z, bottom='Point_FCC_[111]', top='Point_FCC_[111]', metal='Au'):
    """Electrodes and molecule placed as with the GUI sliders: molecule z and top-electrode z, everything else 0."""
    down = antui.place_bottom_electrode(antui.xyz_object(f'ElectrodesDFT/{bottom}.xyz'))
    up = antui.xyz_object(f'ElectrodesDFT/{top}.xyz')
    if metal != 'Au':
        down.change_elem(metal)
        up.change_elem(metal)
    up.spatial_coords = np.array([0., 0., top_z])
    mol = antui.xyz_object(f'Molecules/{molecule}.xyz')
    mol.spatial_coords = np.array([0., 0., mol_z])
    return [down, mol, up]


def export(folder, kind, create, objs, *args, opt=OPT_NONE, dft=None, params=None):
    dft = dft or {'Functional': 'HSE06', 'Shell': 'Close Shell', 'DFTU': False, 'EABS': False}
    target = os.path.join(HERE, folder)
    shutil.rmtree(target, ignore_errors=True)
    os.makedirs(target)
    name = objs[1].name.replace(' ', '_') + '-0001'
    Dir = f'Examples/{folder}/{name}'
    os.makedirs(Dir)
    lines, warnings = describe_export(objs, kind, opt, dft, params)
    create(objs, Dir, name, *args)
    _OUT.verify_output(Dir)
    with open(f'{Dir}/ANTUI_summary.txt', 'w') as out:
        out.write(summary_text(lines, warnings, Dir, f'v{antui.__version__}'))
    print(f'{Dir}: {len(warnings)} warning(s)')


def main():
    hse = {'Functional': 'HSE06', 'Shell': 'Close Shell', 'DFTU': False, 'EABS': False}
    pbe = dict(hse, Functional='PBE')

    # 1. Single junction: flat benzene between two Au(111) tips (apex gap 5.76 A), single-point ANT input
    export('1_single_junction', 'ANT', _OUT.Create_output, junction('Benzene', -3.0, 5.2), OPT_NONE, hse, dft=hse)

    # 2. Fig. 2(a): symmetric pull of the same junction, 20 configurations, gap +0.1 A per step,
    #    molecule relaxed at every step
    opt = dict(OPT_NONE, Molecule='all', Chains='every step')
    export('2_pull_every_step', 'pull', _OUT.Create_pull_output, junction('Benzene', -3.0, 5.2), 0.1, 20, True, opt, hse,
           opt=opt, dft=hse, params={'dz': 0.1, 'nz': 20, 'both': True})

    # 3. Fig. 2(b): glycerol bound to both tips, rotated 0-90 degrees in 15-degree steps, first-step relaxation
    opt = dict(OPT_NONE, Molecule='all', Chains='first step')
    angles = antui.rotation_angles(0, 90, 15)
    export('3_rotation_first_step', 'rotation', _OUT.Create_rot_output, junction('Glycerol', -2.0, 6.6), angles,
           (False, True, False), False, opt, pbe,
           opt=opt, dft=pbe, params={'angles': angles, 'components': (False, True, False), 'soc': False})

    # 4. Fig. 2(c): toluene lying on an Au(111) slab, 10 x 10 tip scan over 8 x 8 A at constant height,
    #    first-step relaxation
    opt = dict(OPT_NONE, Molecule='all', Chains='first step')
    grid = [[-4.0, 4.0, 10], [-4.0, 4.0, 10], 3.7]
    export('4_grid_first_step', 'grid', _OUT.Create_grid_output,
           junction('Toluene', -6.7, 3.7, bottom='Block_FCC_[111]_40'), grid, opt, pbe,
           opt=opt, dft=pbe, params={'grid': grid})

    # 5. Get U: PBE single points of isolated glycerol with N, N-1 and N+1 electrons and the Get_U.py analysis
    export('5_get_u', 'getU', _OUT.Create_getU_output, junction('Glycerol', -2.0, 6.6), pbe, dft=pbe)


if __name__ == '__main__':
    main()
