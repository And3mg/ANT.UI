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
    win.confirm_export.set(False)    # tests write directly; the summary dialog has its own test
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
        window.export(broken, kind='xyz')
    assert set(os.listdir('Outputs')) == before


def test_output_numbering_skips_deleted_folders(antui, tmp_path):
    """X6: counting folders reused a number after a deletion."""
    for n in ('Benzene-0001', 'Benzene-0003', 'Benzene_dimer-0007'):
        (tmp_path / n).mkdir()
    assert antui.next_output_name('Benzene', folder=str(tmp_path)) == 'Benzene-0004'
    assert antui.next_output_name('Toluene', folder=str(tmp_path)) == 'Toluene-0001'


def test_every_export_saves_its_summary(window):
    folder = window.create_ant()
    text = open(f'{folder}/ANTUI_summary.txt').read()
    assert 'Calculation: ANT input' in text and f'Destination: {folder}' in text


def test_output_folder_can_be_anywhere(window, tmp_path):
    target = (tmp_path / 'my calcs').as_posix()      # with a space, created on demand
    previous = window.output_root
    try:
        window.set_output_root(target)
        folder = window.create_pull_dft(0.25, 3, True)
    finally:
        window.set_output_root(previous)
    assert folder.startswith(target) and os.path.isdir(folder)
    assert f'name={os.path.basename(folder)}' in open(f'{folder}/launchANT.sh').read()


def test_charge_fields_reach_the_input(window):
    from antui_testlib import Gjf
    window.entry_charge.insert(0, '-1')
    window.entry_mult.insert(0, '2')
    window.combo_shell.set('Open Shell')
    window.update_DFT(None)
    try:
        folder = window.create_ant()
    finally:
        window.entry_charge.delete(0, 'end')
        window.entry_mult.delete(0, 'end')
        window.combo_shell.set('Close Shell')
        window.update_DFT(None)
    name = os.path.basename(folder)
    g = Gjf(f'{folder}/input{name}.gjf')
    assert (g.charge, g.mult) == (-1, 2)


def test_invalid_charge_is_refused(window):
    window.entry_charge.insert(0, 'one')
    try:
        with pytest.raises(ValueError):
            window.create_ant()
    finally:
        window.entry_charge.delete(0, 'end')


def test_summary_dialog_writes_only_on_generate(window):
    before = set(os.listdir(window.output_root))
    window.confirm_export.set(True)
    try:
        assert window.create_ant() is None
        dialog = window.summary_window
        assert 'Calculation: ANT input' in dialog.textbox.get('0.0', 'end')
        assert set(os.listdir(window.output_root)) == before          # nothing written yet
        folder = dialog.generate()
    finally:
        window.confirm_export.set(False)
    assert os.path.isdir(folder)
    window.confirm_export.set(True)
    try:
        window.create_ant()
        window.summary_window.close()                                  # Cancel
    finally:
        window.confirm_export.set(False)
    assert set(os.listdir(window.output_root)) - before == {os.path.basename(folder)}
    assert window.errors == [], window.errors[0]


def test_site_config_panel_shows_the_placeholder_warning(window, antui):
    panel = antui.toplevel_site_config()
    try:
        text = panel.textbox.get('0.0', 'end')
        assert 'G09call: g09' in text and 'JOBHEAD' in text and '!' in text
        panel.reload()
    finally:
        panel.window.destroy()


def test_assistant_previews(window, antui):
    rot = antui.toplevel_rotate(lambda *a: None)
    pull = antui.toplevel_pull(lambda *a: None)
    try:
        assert rot.label_preview.cget('text').startswith('7 configurations: 0, 15, 30')
        rot.text_by.delete(0, 'end'); rot.text_by.insert(0, '20')
        rot.text_from.delete(0, 'end'); rot.text_from.insert(0, '10')
        rot.text_to.delete(0, 'end'); rot.text_to.insert(0, '70')
        rot.preview()
        assert rot.label_preview.cget('text').startswith('4 configurations: 10, 30, 50, 70')
        pull.check_both.select(); pull.preview()
        assert 'top +0.25 A, bottom -0.25 A' in pull.label_preview.cget('text')
    finally:
        rot.window.destroy()
        pull.window.destroy()
