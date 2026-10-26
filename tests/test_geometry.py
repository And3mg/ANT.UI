"""Junction geometry must be a function of the GUI state, not of its history."""
import numpy as np

from conftest import known_bug


def _sorted(coords):
    return coords[np.lexsort(np.round(coords, 3).T)]


@known_bug('X4', 'start-up bottom electrode is rotated (y180, z60); after re-selection it is mirrored')
def test_bottom_electrode_same_after_reselection(antui):
    startup = antui.xyz_object('ElectrodesDFT/Point_FCC_[111].xyz')
    startup.euler_angles = np.array([0, 180, 60])
    startup.spatial_coords = np.array([0, 0, -11])
    reselected = antui.xyz_object('ElectrodesDFT/Point_FCC_[111].xyz')
    reselected.reflection_z()
    reselected.spatial_coords = np.array([0, 0, -11])
    np.testing.assert_allclose(_sorted(startup.true_coordinates()), _sorted(reselected.true_coordinates()), atol=1e-3)


@known_bug('X5', 'true_coordinates applies rotation increments, so the result depends on slider order')
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
