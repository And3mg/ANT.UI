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
  component's own centre. They now turn about the axis through the centre the component has in the viewer, also
  when the optimisation moves it sideways (the relaxed centroid would drag the contact atoms off their apices).
  Rotation angles are now relative to the orientation shown in the viewer.
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
- The + and - buttons of the position and angle sliders labelled the new value one step too far, so the label
  seemed to jump by 0.2 (or 0.3 when changing direction) although the geometry moved by 0.1. With a `ZMAX` that is
  not a whole number of steps from -20, the slider positions were also off the 0.1 A grid.
- A molecule file with a blank line (e.g. at the end) could not be read, and choosing it removed the molecule
  controls and stopped the window from refreshing until restart. Blank lines are now skipped; an unreadable file
  gives one error message and the previous molecule stays selected.
- The SOC button was never disabled (customtkinter ignores `widget['state'] = ...`), so Ni junctions failed with
  `KeyError: 'Ni'`. It is now disabled when a metal has no SOC factors, and SOC exports explain why.
  The DFT+U box is now also disabled for functionals without a stored U.
- Reopening Optimization Options and pressing Ok switched "Opt. in ANT" off.
- Grid Assistant: Ok used the values of the last preview refresh (every 0.5 s) and silently kept the previous
  value when an entry was invalid (e.g. a decimal comma). It now reads the entries and reports invalid ones.
- z-only atoms were kept when another molecule was chosen; indices beyond the new molecule constrained electrode
  atoms instead. They are now cleared with the molecule, and invalid indices are refused.
- The 3D viewer redraws only when the junction changes. It used to redraw every cycle, which blocked the controls
  and made a dragged window jump or snap back to where it started.
- Molecule comment lines: only "FUNCTIONAL U" pairs are read; other text (e.g. leftover `Energy:` lines in
  Benzene, Toluene and Cyclohexane, now removed) is ignored.

### Changed
- SOC factors of Au, Pt and Ag recalibrated against two-component CRYSTAL23 bands computed with the same basis
  sets and pseudopotentials (RMS band error 3–10 meV). The previous factors gave Au spin-orbit constants about 20 %
  too weak, and Pt and Ag constants too strong (d shells by 20–30 %, p shells by 40–110 %). RCUT is now 2.0 bohr for
  all three (it was 2.5 for Pt and Ag). Al and Cu keep 1, 1, 0 (no rescaling). `SOC_params/antg-basis-sets-and-socfacs.pdf`
  documents the method and every element.
- Au electrodes use the `Papacon` Bethe-lattice set instead of `lsda-crenbs`. ANT reads `lsda-crenbs` as an
  orthogonal tight-binding set, which gives the Au electrode 2.4 times the bulk DOS at E_F (0.44 vs 0.19 states/eV per
  spin), a band about twice as wide, a dip in T(E) at E_F + 3.5 eV, a Fermi level 0.35-0.50 eV higher, and changes
  T(E_F) by +1.5 % (atomic nanocontact) to -32 % (5 A tunnelling gap). Au inputs written by 1.0.0 should be regenerated.

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
- Get U button: Gaussian single points of the isolated molecule with N, N-1 and N+1 electrons and a `Get_U.py`
  script that computes U = (IP - EA) - (e_LUMO - e_HOMO) and writes the molecule file with the U of the chosen
  functional in its comment line.

## 1.0.0 (2026)

First public version.
