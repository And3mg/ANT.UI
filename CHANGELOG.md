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
- Stored PBE U of glycerol: 11.86 eV is replaced by 5.88 eV, the value given by (IP - EA) - (e_LUMO - e_HOMO) from
  its gas-phase runs. Inputs written with the old value should be regenerated.
- The program can be started from any folder.
- Dragging the control window could abort Python on machines with Qt bindings installed, because matplotlib then
  opened the 3D viewer with Qt. The viewer now always uses Tk, like the controls. Closing the viewer leaves the
  controls working.
- The 3D viewer redraws only when the junction changes. It used to redraw every cycle, which blocked the controls
  and made a dragged window jump or snap back to where it started.
- Molecule comment lines: only "FUNCTIONAL U" pairs are read; other text (e.g. leftover `Energy:` lines in
  Benzene, Toluene and Cyclohexane, now removed) is ignored.

### Changed
- SOC factors of Au, Pt and Ag recalibrated against two-component CRYSTAL23 bands computed with the same basis
  sets (band error 3–10 meV). The previous fit counted part of the SOC twice and gave factors that were too small.
  RCUT is now 2.0 bohr for all three. Al and Cu keep 1, 1, 0. `SOC_params/antg-basis-sets-and-socfacs.pdf`
  documents the method and every element.

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
- Warning for mixed-metal junctions that ANT.Gaussian 2.8.2 and earlier read wrongly: the top electrode's
  Bethe-lattice parameters inherit the energy unit (and ISPIN) of the bottom ones, so for example a Pt bottom and an
  Au top halve every Au parameter. Swap the electrodes or use an ANT build with the fix.
- Get U button: Gaussian single points of the isolated molecule with N, N-1 and N+1 electrons and a `Get_U.py`
  script that computes U = (IP - EA) - (e_LUMO - e_HOMO) and writes the molecule file with the U of the chosen
  functional in its comment line.

## 1.0.0 (2026)

First public version.
