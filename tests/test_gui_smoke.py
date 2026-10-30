"""Open the real window, drive it through its own callbacks, close it."""
import os
import shutil
import time
import traceback

import numpy as np
import pytest

pytestmark = pytest.mark.gui


@pytest.fixture(scope='module')
def gui(antui):
    """One window for the module: Tk cannot create a new root after one was destroyed in this process."""
    import tkinter
    import customtkinter
    import matplotlib.pyplot as plt

    try:
        root = customtkinter.CTk()
    except tkinter.TclError as e:
        pytest.skip(f'no display: {e}')
    errors = []
    root.report_callback_exception = lambda *exc: errors.append(''.join(traceback.format_exception(*exc)))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    win = antui.call_back_dft(root, ax, fig)
    root.update()
    win.errors = errors
    win.root = root
    yield win
    if win.after_process:
        root.after_cancel(win.after_process)
    root.destroy()
    plt.close(fig)


@pytest.fixture
def window(gui):
    before = set(os.listdir('Outputs'))
    gui.errors.clear()
    yield gui
    for folder in set(os.listdir('Outputs')) - before:
        shutil.rmtree(f'Outputs/{folder}', ignore_errors=True)


def let_gui_refresh(win, seconds=1.2):
    """The window applies combobox changes in its 500 ms refresh loop."""
    end = time.time() + seconds
    while time.time() < end:
        win.root.update()
        time.sleep(0.02)


def test_window_starts_and_exports(window):
    created = [window.create_xyz(), window.create_ant(), window.create_pull_dft(0.25, 4, False),
               window.create_rot_dft(0, 30, 15, False, True, False, False)]
    assert all(os.path.isdir(d) for d in created)
    assert len(set(created)) == 4, 'every export gets its own folder'
    assert 'configurations' in window.status.cget('text')
    assert window.errors == [], window.errors[0]


def test_reselected_bottom_electrode_is_unchanged(window):
    """X4 through the real combobox: switch the bottom electrode away and back."""
    start = window.down_electrode.true_coordinates().copy()
    window.combo_bottom_type.set('Point_FCC_[100]_4')
    let_gui_refresh(window)
    assert window.down_electrode.name == 'Point_FCC_[100]_4'
    window.combo_bottom_type.set('Point_FCC_[111]')
    let_gui_refresh(window)
    np.testing.assert_allclose(window.down_electrode.true_coordinates(), start, atol=1e-9)
    assert window.errors == [], window.errors[0]


def test_rotation_assistant_leaves_the_junction_as_it_was(window):
    """X11: the assistant used to leave the (int-typed) bottom electrode rotated after export."""
    start = window.down_electrode.true_coordinates().copy()
    window.create_rot_dft(0, 22.5, 7.5, False, False, True, False)
    np.testing.assert_allclose(window.down_electrode.true_coordinates(), start, atol=1e-9)


def test_failed_export_leaves_no_folder(window):
    """X8: a generator error removes the half-written folder and reaches the caller (error dialog)."""
    before = set(os.listdir('Outputs'))

    def broken(objs, Dir, name):
        open(f'{Dir}/partial.gjf', 'w').close()
        raise RuntimeError('boom')

    with pytest.raises(RuntimeError):
        window.export(broken)
    assert set(os.listdir('Outputs')) == before


def test_output_numbering_skips_deleted_folders(antui, tmp_path):
    """X6: counting folders reused a number after a deletion."""
    for n in ('Benzene-0001', 'Benzene-0003', 'Benzene_dimer-0007'):
        (tmp_path / n).mkdir()
    assert antui.next_output_name('Benzene', folder=str(tmp_path)) == 'Benzene-0004'
    assert antui.next_output_name('Toluene', folder=str(tmp_path)) == 'Toluene-0001'
