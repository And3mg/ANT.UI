# ANT.UI

**ANT.UI** is a Python graphical interface for building molecular junction geometries and generating ready-to-run input files for [ANT.Gaussian](https://github.com/juanjosepalacios/ANT.Gaussian) (NEGF-DFT transport), [Gaussian 09/16](https://gaussian.com), and XYZ coordinate files compatible with external visualizers such as OVITO.

It turns what was previously hours of manual scripting into a point-and-click workflow: load your electrodes and molecule, position them interactively in 3D, set your DFT parameters, and click to generate all required input files.

---

## Requirements

- Python ≥ 3.8
- [ANT.Gaussian](https://github.com/juanjosepalacios/ANT.Gaussian) installed on your HPC cluster

### Desktop dependencies (for the GUI)

```bash
pip install -r requirements.txt
```

This installs `numpy`, `matplotlib`, and `customtkinter`. `tkinter` is included with the standard Python installer on Windows and macOS; on Linux install it via your package manager (e.g. `sudo apt install python3-tk`).

### Cluster-side dependency

The Python helper scripts that ANT.UI generates for sequential jobs (optimization chains, pull sequences, grid scans) only require `numpy`:

```bash
pip install -r requirements-hpc.txt
```

---

## Quick start

1. **Configure `USER_CONF.txt`** — set `ANTDIRcall` to the path of your ANT installation on the cluster, and fill in `JOBHEAD` and `GAUSSIANHEAD` to match your cluster's submission system. See the manual for all available keywords.

2. **Launch the application:**
   ```bash
   python ANT.UI.py
   ```
   It can be started from any folder. On Windows you can also double-click `ANT.UI.bat` (or send it to the desktop as a shortcut).

3. **Build your junction** — select electrodes and molecule, position them with the sliders, set DFT options, and press **Create ANT** (or **Create xyz** / **SOC ANT**) to generate the output folder.

4. **Transfer outputs to your cluster** and run.

---

## Running the tests

The test suite checks the generated inputs without Gaussian or a cluster: every atom has a basis set (and an ECP when its basis needs one), the charge and multiplicity agree with the electron count, the generated helper scripts are valid Python, and the sequence assistants produce every step. For optimisation chains it writes the Gaussian log an unchanged optimisation would leave, runs the generated `link_opt.py`, and checks that the propagated geometries match the requested scan.

```bash
pip install pytest
python -m pytest
```

Run it from the repository root. The same tests run on every push on Windows and Linux (GitHub Actions).

---

## File structure

| Path | Contents |
|---|---|
| `ANT.UI.py` | Main application entry point |
| `ANT.UI.bat` | Windows launcher (double-click; no console window) |
| `_OUT.py` | Output file generation (Gaussian `.gjf`, `.ant`, `.sh`, `.xyz`) |
| `_DFT.py` | DFT file writers and Python connector scripts |
| `_SOC.py` | Spin-orbit coupling output generation |
| `_SUMMARY.py` | Pre-export summary and checks (saved as `ANTUI_summary.txt` with every output) |
| `CHANGELOG.md` | Changes between versions |
| `_DICS.py` | Element dictionaries, functional list, user config loader |
| `Examples/` | Inputs of the manual's worked examples (the three cases of Fig. 2 of the paper) and the script that regenerates them |
| `tests/` | Automated tests of the generated inputs (`python -m pytest`) |
| `ElectrodesDFT/` | Pre-built electrode geometries (XYZ, top-electrode convention) |
| `Molecules/` | Molecule library (XYZ format) |
| `DFT_basis_pseudopot/` | Basis sets for electrode atoms shown in the viewer |
| `DFT_min_basis_pseudopot/` | Minimal basis sets for Bethe-lattice support atoms |
| `SOC_params/` | Basis sets and SOC parameters for supported elements |
| `PythonRutines/` | Template scripts embedded in sequential outputs |
| `USER_CONF.txt` | User configuration (cluster paths, Gaussian header, etc.) |

---

## Documentation

A full user manual is available in `ANT.UI_UserManual.pdf` in the repository.

---

## Citation

If you use ANT.UI in your research, please cite:

> A. Martinez-Garcia, J. J. Palacios, C. Sabater, *ANT.UI: An interactive 3D tool for preparing ANT.Gaussian molecular junction geometries*, SoftwareX (2026). *(DOI to be added upon publication)*

The software itself is archived on Zenodo; `CITATION.cff` gives the citation of each release.

---

## License

MIT — see [LICENSE](LICENSE).

## Authors

- Andrés Martinez-Garcia — [andres.martinez@ua.es](mailto:andres.martinez@ua.es)
- Carlos Sabater — [carlos.sabater@ua.es](mailto:carlos.sabater@ua.es)

Departamento de Física & Instituto Universitario de Materiales de Alicante (IUMA), Universidad de Alicante, Spain.
