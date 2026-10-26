"""Junction geometry must be a function of the GUI state, not of its history."""
import numpy as np


def _sorted(coords):
    return coords[np.lexsort(np.round(coords, 3).T)]


def test_bottom_electrode_is_placed_the_same_way_every_time(antui):
    """X4: start-up and re-selection used different transforms. The GUI path is in test_gui_smoke."""
    first = antui.place_bottom_electrode(antui.xyz_object('ElectrodesDFT/Point_FCC_[111].xyz'))
    second = antui.place_bottom_electrode(antui.xyz_object('ElectrodesDFT/Point_FCC_[111].xyz'))
    np.testing.assert_allclose(first.true_coordinates(), second.true_coordinates(), atol=1e-9)
    assert first.true_coordinates()[:, 2].max() < 0, 'bottom electrode must point up from below the molecule'


def _stacking(coords, ref, a1, a2):
    """A/B/C site of every (111) layer, in ascending z."""
    letters = ''
    for z in np.unique(np.round(coords[:, 2], 2)):
        atom = coords[np.isclose(np.round(coords[:, 2], 2), z)][0]
        f = np.linalg.solve(np.array([a1, a2]).T, atom[:2] - ref)
        letters += 'ACB'[int(round(3 * (f[0] + f[1]))) % 3]
    return letters


def test_bottom_111_continues_the_fcc_stacking(antui):
    """X4: the [111] bottom tip is the top tip rotated (y180, z60), so ABC stacking continues across the
    gap (one crystal). The mirror used after re-selecting the electrode reversed it (a twin boundary)."""
    top = antui.xyz_object('ElectrodesDFT/Point_FCC_[111].xyz')
    top.spatial_coords = np.array([0., 0., 9.])
    t = top.true_coordinates()
    b = antui.place_bottom_electrode(antui.xyz_object('ElectrodesDFT/Point_FCC_[111].xyz')).true_coordinates()
    nn = 2.84  # in-plane Au-Au distance of the template; a1 along x
    a1, a2 = nn * np.array([1., 0.]), nn * np.array([0.5, np.sqrt(3) / 2])
    ref = t[np.argmin(t[:, 2]), :2]
    bottom, upper = _stacking(b, ref, a1, a2), _stacking(t, ref, a1, a2)
    assert bottom == upper == 'ABCA', (bottom, upper)


def test_orientation_depends_only_on_slider_values(antui):
    direct = antui.xyz_object('Molecules/Benzene.xyz')
    direct.euler_angles = np.array([30., 40., 0.])
    target = direct.true_coordinates().copy()

    stepwise = antui.xyz_object('Molecules/Benzene.xyz')
    stepwise.euler_angles = np.array([0., 40., 0.])
    stepwise.true_coordinates()
    stepwise.euler_angles = np.array([30., 40., 0.])
    np.testing.assert_allclose(stepwise.true_coordinates(), target, atol=1e-6)


def test_element_change_round_trip(antui):
    electrode = antui.xyz_object('ElectrodesDFT/Point_FCC_[111].xyz')
    original = electrode.poss.copy()
    electrode.change_elem('Cu')
    electrode.change_elem('Au')
    np.testing.assert_allclose(electrode.poss, original, atol=1e-9)
