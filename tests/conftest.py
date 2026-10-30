"""Shared fixtures. ANT.UI reads its data files relative to the working
directory, so the whole session runs from the repository root."""
import gc
import importlib.util
import itertools
import os
import shutil
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault('MPLBACKEND', 'Agg')


def known_bug(bug_id, reason, condition=True):
    """Mark a test that documents an open bug. strict: a fix makes it fail until the mark is removed."""
    return pytest.mark.xfail(condition, reason=f'{bug_id}: {reason}', strict=True)


@pytest.fixture(scope='session')
def antui():
    spec = importlib.util.spec_from_file_location('antui', os.path.join(ROOT, 'ANT.UI.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def junction(antui):
    """Build [bottom, molecule, top] exactly as tkinter_window does at start-up."""
    def build(mol='Benzene', top='Au', bottom='Au', top_z=9.0, electrode='Point_FCC_[111]', drop_top=0):
        up = antui.xyz_object(f'ElectrodesDFT/{electrode}.xyz')
        down = antui.place_bottom_electrode(antui.xyz_object(f'ElectrodesDFT/{electrode}.xyz'))
        if top != 'Au':
            up.change_elem(top)
        if bottom != 'Au':
            down.change_elem(bottom)
        up.spatial_coords = np.array([0., 0., top_z])
        up.euler_angles = np.array([0., 0., 0.])
        for _ in range(drop_top):       # remove the tip apex to change the metal-atom parity
            keep = np.ones(up.N, bool)
            keep[np.argmin(up.poss[:, 2])] = False
            for attr in ('poss', 'poss0', 'elements', 'semi_infinite_layer', 'extremes', 'neigh_count'):
                if hasattr(up, attr):
                    setattr(up, attr, getattr(up, attr)[keep])
            up.N -= 1
        molecule = antui.xyz_object(f'Molecules/{mol}.xyz')
        molecule.euler_angles = np.array([0., 0., 0.])
        return [down, molecule, up]
    return build


_counter = itertools.count(1)


def _sweep_test_outputs():
    """Remove leftovers. Folders of crashed generators stay locked while pytest holds their tracebacks."""
    gc.collect()
    for folder in os.listdir('Outputs'):
        if folder.startswith('zzpytest'):
            shutil.rmtree(f'Outputs/{folder}', ignore_errors=True)


def pytest_sessionstart(session):
    _sweep_test_outputs()


def pytest_unconfigure(config):
    _sweep_test_outputs()


@pytest.fixture
def outdir():
    """Return a factory giving (Dir, name) laid out like the GUI: Outputs/<name>."""
    made = []

    def make():
        name = f'zzpytest{next(_counter)}-0001'
        path = f'Outputs/{name}'
        shutil.rmtree(path, ignore_errors=True)
        os.makedirs(path)
        made.append(path)
        return path, name
    yield make
    gc.collect()    # generators that crash leave files open; Windows cannot delete open files
    for path in made:
        shutil.rmtree(path, ignore_errors=True)


@pytest.fixture
def opt_data():
    def make(**changes):
        data = {'Top_Bethe': 'fix', 'Top': 'fix', 'Bottom_Bethe': 'fix', 'Bottom': 'fix',
                'Molecule': 'fix', 'Molecule_z_atoms': [], 'Chains': 'none', 'Opt_steps': 10,
                'Opt_in_ANT': False}
        data.update(changes)
        return data
    return make


@pytest.fixture
def dft_opts():
    def make(**changes):
        data = {'Functional': 'HSE06', 'Shell': 'Close Shell', 'DFTU': False, 'EABS': False}
        data.update(changes)
        return data
    return make
