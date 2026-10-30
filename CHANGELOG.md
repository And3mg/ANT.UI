# Changelog

## 1.1.0 (2026-11)

Revision for the SoftwareX review. Inputs written by 1.0.0 for SOC, mixed-metal, Pt/Ni or optimisation-chain
calculations should be regenerated.

### Fixed
- SOC with two different metals: the pseudopotential of the active atoms of the first metal was missing.
- SOC with molecular optimisation: the generated `link_opt.py` was not valid Python.
- Pull, Rotation and Grid assistants with first-step optimisation: the connector scripts failed (unquoted folder
  names; `np.float64(...)` values under numpy 2).
- First-step chains moved every propagated configuration one increment too far, and symmetric pulls halved the
  step at every configuration.
- Rotation chains turned clockwise (the viewer turns counter-clockwise) and about the origin instead of the
  component's own centre. Rotation angles are now relative to the orientation shown in the viewer.
- Every-step grid scans lagged one point behind.
- Symmetric pulls with decimal steps (e.g. dz = 0.1, nz = 3) stopped with an `IndexError`.
- Bethe-lattice layers in z mode: the second plane got one z variable per atom. Gaussian now relaxes one z per
  plane and `link_opt.py` restores the original spacing between the Bethe planes.
- Charge and multiplicity are computed from the electron count after the pseudopotential cores (Pt and Ni
  junctions with an odd number of atoms were rejected by Gaussian).
- Atoms of elements missing from the element table were silently dropped from the Gaussian input.
- Inline `#` comments in `USER_CONF.txt` broke the route line of optimisation chains.
- The bottom electrode differed after re-selecting it; the orientation depended on the order in which the sliders
  were moved; output folder numbers could collide after a deletion; settings leaked between exports.
- The program can be started from any folder.

### Added
- Pre-export summary of the scientific setup (atoms per region, contacts, electrons, charge and multiplicity, basis
  set and pseudopotential per element and region, SOC, DFT+U, constraints, sequence settings, warnings), shown
  before writing and saved as `ANTUI_summary.txt` in every output folder.
- User-selectable output folder (`OUTPUTDIR` in `USER_CONF.txt` for the default).
- Charge and Multiplicity fields (empty = automatic).
- Site configuration panel and a warning for the placeholder job header.
- Live previews in the Pull and Rotation assistants ("step (deg)").
- Error dialog; generated scripts are syntax-checked and incomplete output folders are removed.
- Automated test suite (`python -m pytest`) and GitHub Actions on Windows and Linux, numpy 1.x and 2.x.
- Worked examples (`Examples/`) and a revised user manual.
- `ANT.UI.bat` launcher for Windows.

## 1.0.0 (2026)

First public version.
