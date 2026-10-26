"""Open the real window, press its export buttons through their callbacks, close it."""
import os
import shutil
import traceback

import pytest

pytestmark = pytest.mark.gui


def test_window_starts_and_exports(antui):
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
    before = set(os.listdir('Outputs'))
    try:
        window = antui.call_back_dft(root, ax, fig)
        root.update()
        window.create_xyz()
        window.create_ant()
        window.create_pull_dft(0.25, 4, False)
        root.update()
        created = set(os.listdir('Outputs')) - before
        assert len(created) == 3, created
        assert errors == [], errors[0]
    finally:
        for folder in set(os.listdir('Outputs')) - before:
            shutil.rmtree(f'Outputs/{folder}', ignore_errors=True)
        if getattr(window, 'after_process', None):
            root.after_cancel(window.after_process)
        root.destroy()
        plt.close(fig)
