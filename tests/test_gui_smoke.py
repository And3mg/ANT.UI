"""Open the real window, drive it through its own callbacks, close it."""
import os
import shutil
import subprocess
import sys
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


def run_main(mainloop):
    """Run ANT.UI.main() with this body as the window's mainloop, in its own process: Tk cannot create
    a new root after one was destroyed."""
    code = f'''
import importlib.util, tkinter
import customtkinter, matplotlib
import matplotlib.pyplot as plt

tk_mainloop = customtkinter.CTk.mainloop
def mainloop(self):
{mainloop}

customtkinter.CTk.mainloop = mainloop
spec = importlib.util.spec_from_file_location('antui', 'ANT.UI.py')
antui = importlib.util.module_from_spec(spec)
spec.loader.exec_module(antui)
antui.main()
'''
    return subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=30)


def test_main_uses_the_tk_viewer_and_a_ctk_default_root(gui):
    """With Qt installed matplotlib picked it, and dragging the window aborted Python inside flush_events."""
    out = run_main('''
    print(matplotlib.get_backend(), type(tkinter._default_root).__name__)
    self.destroy()
''')
    assert out.stdout.split()[-2:] == ['TkAgg', 'CTk'], out.stderr


def test_closing_the_viewer_leaves_the_controls_running(gui):
    """The refresh loop drew on the destroyed Tk viewer: segmentation fault."""
    out = run_main('''
    self.after(1000, lambda: plt.close('all'))
    self.after(3000, lambda: (print('still running'), self.destroy()))
    tk_mainloop(self)
''')
    assert out.returncode == 0 and 'still running' in out.stdout, out.stderr


def test_closing_the_controls_ends_the_program_with_the_viewer_open(gui):
    """The Tk loop waits for every Tk window: the process stayed alive with the viewer."""
    out = run_main('''
    self.after(1000, lambda: self.tk.call(self.protocol('WM_DELETE_WINDOW')))
    tk_mainloop(self)
    print('finished')
''')
    assert 'finished' in out.stdout, out.stderr


def test_viewer_redraws_only_when_the_junction_changes(window, antui, monkeypatch):
    """A full 3D redraw every cycle blocked the window: dragging it jumped and could snap back."""
    let_gui_refresh(window)
    draws = []
    monkeypatch.setattr(antui.fig.canvas, 'draw', lambda: draws.append(1))
    let_gui_refresh(window, 1.5)
    assert draws == []
    z = window.mol_controls.slider_z.get()
    try:
        window.mol_controls.slider_z.set(z + 1)
        let_gui_refresh(window)
        assert len(draws) == 1
    finally:
        window.mol_controls.slider_z.set(z)
        let_gui_refresh(window)
    assert window.errors == [], window.errors[0]


def test_plus_minus_buttons_move_one_step_and_the_label_follows(window):
    """The +/- buttons moved the slider one step but labelled it one step further: the label jumped by
    0.2 on the first press and by 0.3 when the direction changed."""
    for scale, step in ((window.mol_controls.slider_x, 0.1), (window.mol_controls.slider_phi, 1)):
        start = scale.get()
        try:
            values = []
            for button in ('+', '+', '+', '-', '-', '+'):
                (scale.slider_plus if button == '+' else scale.slider_minus).invoke()
                values.append(scale.get())
                shown = float(scale.label.cget('text').split('=')[1])
                assert shown == pytest.approx(scale.get(), abs=1e-9), scale.label.cget('text')
            np.testing.assert_allclose(np.diff([start] + values), step * np.array([1, 1, 1, -1, -1, 1]), atol=1e-9)
        finally:
            scale.set(start)
    assert window.errors == [], window.errors[0]


def test_slider_positions_stay_on_the_step_grid(window, antui):
    """(25.3 + 20)/0.1 = 452.999...: int() gave one step fewer, so the positions were 0.1002 A apart."""
    import customtkinter
    frame = customtkinter.CTkFrame(window.root)
    try:
        scale = antui.Scale(frame, 'z', -20, 25.3, 0.1)
        scale.set(0)
        for _ in range(7):
            scale.slider_plus.invoke()
        assert scale.get() == pytest.approx(0.7, abs=1e-9)
        assert scale.label.cget('text') == 'z = 0.7'
        scale.set(-1e-12)
        assert scale.label.cget('text') == 'z = 0.0'      # not "-0.0"
    finally:
        frame.destroy()


@pytest.fixture
def extra_molecule():
    """Write Molecules/zzpytest_<name>.xyz for one test."""
    made = []

    def make(name, text):
        path = f'Molecules/zzpytest_{name}.xyz'
        with open(path, 'w') as f:
            f.write(text)
        made.append(path)
        return f'zzpytest_{name}'
    yield make
    for path in made:
        os.remove(path)


WATER = '3\n\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'


def test_molecule_file_with_blank_lines_loads(antui, extra_molecule):
    name = extra_molecule('blank', WATER.replace('H 0.96', '\nH 0.96') + '\n\n')
    mol = antui.xyz_object(f'Molecules/{name}.xyz')
    assert list(mol.elements) == ['O', 'H', 'H']


def test_unreadable_molecule_keeps_the_window_working(window, antui, extra_molecule):
    """An unreadable molecule file removed the molecule controls and stopped the refresh loop."""
    name = extra_molecule('short', '5\n\nO 0 0 0\nH 0.96 0 0\n')
    start = window.molecule.name
    window.combo.set(name)
    let_gui_refresh(window)
    assert len(window.errors) == 1 and 'says 5 atoms' in window.errors[0]
    assert window.combo.get() == start and window.molecule.name == start
    window.errors.clear()
    window.combo.set('Toluene')
    let_gui_refresh(window)
    window.mol_controls.slider_z.set(1.5)
    let_gui_refresh(window)
    assert window.molecule.name == 'Toluene' and window.molecule.spatial_coords[2] == pytest.approx(1.5)
    window.mol_controls.slider_z.set(0)
    window.combo.set(start)
    let_gui_refresh(window)
    assert window.errors == [], window.errors[0]


def test_soc_needs_factors_for_every_metal(window):
    """Ni has no SOC factors: the button stayed enabled (widget['state'] = ... does nothing in customtkinter)
    and the export failed with KeyError: 'Ni'."""
    try:
        window.combo_top_elem.set('Ni')
        let_gui_refresh(window)
        assert window.soc_btn.cget('state') == 'disabled'
        with pytest.raises(ValueError, match = 'No SOC factors for Ni'):
            window.create_soc()
        with pytest.raises(ValueError, match = 'No SOC factors for Ni'):
            window.create_rot_dft(0, 30, 15, True, False, False, True)
    finally:
        window.combo_top_elem.set('Au')
        let_gui_refresh(window)
    assert window.soc_btn.cget('state') == 'normal'


def test_dft_u_box_is_disabled_without_a_stored_u(window):
    window.combo_funct.set('HSE06')
    window.update_DFT(None)
    has_u = 'HSE06' in window.molecule.U
    assert window.check_U.cget('state') == ('normal' if has_u else 'disabled')


def test_options_dialog_keeps_opt_in_ant(window, antui):
    """The checkbox started unticked, so Ok after any other change switched optimisation in ANT off."""
    window.opt_data['Opt_in_ANT'] = True
    try:
        antui.toplevel_opt(window.opt_data, window.molecule, window.update_opt_options).ok()
        assert window.opt_data['Opt_in_ANT'] is True
    finally:
        window.opt_data['Opt_in_ANT'] = False


def test_grid_assistant_uses_what_is_typed(window, antui):
    """Ok returned the values of the last 0.5 s preview, and silently kept them when an entry was invalid."""
    got = []
    g = antui.toplevel_grid(-1, 1, -1, 1, 5.0, window.down_electrode.true_coordinates(),
                            window.molecule.true_coordinates(), got.append)
    try:
        g.update()
        g.text_nx.delete(0, 'end'); g.text_nx.insert(0, '4')
        g.text_xmax.delete(0, 'end'); g.text_xmax.insert(0, '2,5')
        with pytest.raises(ValueError, match = 'decimal point'):
            g.ok()
        g.text_xmax.delete(0, 'end'); g.text_xmax.insert(0, '2.5')
        g.ok()
    finally:
        if g.window.winfo_exists():
            g.window.destroy()
    assert got == [[[-1.0, 2.5, 4], [-1.0, 1.0, 2], 5.0]]


def test_z_atoms_are_reset_with_the_molecule(window, antui):
    """Indices chosen for one molecule were kept for the next; out of range they constrained electrode atoms."""
    start = window.molecule.name
    window.opt_data = dict(window.opt_data, Molecule = 'all, some only z', Molecule_z_atoms = [0, 1])
    try:
        window.combo.set('Toluene' if start != 'Toluene' else 'Benzene')
        let_gui_refresh(window)
        assert window.opt_data['Molecule_z_atoms'] == []
        window.opt_data['Molecule_z_atoms'] = [window.molecule.N]
        with pytest.raises(ValueError, match = 'not in the molecule'):
            window.create_ant()
    finally:
        window.opt_data = dict(window.opt_data, Molecule = 'fix', Molecule_z_atoms = [])
        window.combo.set(start)
        let_gui_refresh(window)
    assert antui.z_atom_indices('3, 1 1', 5) == [1, 3]
    with pytest.raises(ValueError, match = r'\[7\] not in the molecule'):
        antui.z_atom_indices('7', 5)


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


def test_get_u_button_writes_the_gas_phase_inputs(window):
    window.u_btn.invoke()
    folder = window.status.cget('text').split('Created ')[1].split(' ')[0]
    mol = window.molecule.name
    assert sorted(f for f in os.listdir(folder) if f.endswith('.gjf')) == \
           sorted(f'{mol}_{tag}.gjf' for tag in ('N', 'anion', 'cation'))
    assert 'Calculation: Hubbard U of the isolated molecule' in open(f'{folder}/ANTUI_summary.txt').read()
    assert not window.errors


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
