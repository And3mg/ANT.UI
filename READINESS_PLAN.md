# ANT.UI — GitHub & SoftwareX Readiness Plan

Generated: 2026-06-18

---

## 1. Manual vs Code Discrepancies

These are concrete errors or gaps between `ANT.UI_Manual/main.tex` and the actual source code that must be fixed before publication.

### 1.1 L501 → L502 (Critical error)

**Manual (§1, Introduction):** says the `ANTDIRcall` folder contains "the L501 and L101 files needed to call ANT."

**Code (`_DFT.py`, lines 161–162, 249–251, 288–291, 311–313, 345–347):** the actual link written to every `.gjf` and Python connector is `%Subst L502`, not L501.

**Fix:** Replace every mention of L501 with L502 in the manual.

---

### 1.2 REAX button — described but disabled (Medium)

**Manual (§2, last paragraph):** "the REAX input will generate the geometry input for LAMMPS" — written as an active feature.

**Code (`ANT.UI.py`, lines 957–958):** the button is commented out:
```python
#self.reax_btn = Button(self.frame_btn,  text = 'REAX Input', command = self.create_reax)
#self.reax_btn.pack(side = 'right', expand = True)
```
The backend functions `Create_reax_output` and `create_reax` still exist in `_OUT.py` and `ANT.UI.py` but are unreachable from the UI.

**Fix (choose one):**
- Option A: Update the manual to say REAX support is experimental / not available in the current release, and remove the button from the output section description.
- Option B: Re-enable the button for v1.0.0 and test it.

---

### 1.3 Undocumented USER_CONF.txt keywords (Medium)

**Manual (§1):** lists exactly 6 keywords: `ANTDIRcall`, `JOBHEAD`, `GAUSSIANHEAD`, `G09call`, `PYTHONcall`, `ANTpredet`.

**Code (`_DICS.py`, UserInfo dict and parsing loop) + actual `USER_CONF.txt`:** two additional keywords are parsed and used:
- `ZMAX` — controls the upper limit of the z-slider in the UI (`object_controls.__init__`, line 315). Defaults to 20.0 Å. Already present in `USER_CONF.txt` but not documented.
- `U_ELECTRON_SCREENING` — used to scale the DFT+U value written to the `.ant` file. Defaults to 1.0. Not documented anywhere.

**Fix:** Add both keywords to the manual's §1 keyword list, explaining their purpose and defaults.

---

### 1.4 "Abs. Energy" (EABS) checkbox — completely undocumented (Medium)

**Manual:** nowhere mentions an absolute energy option.

**Code (`ANT.UI.py`, lines 928–931):** there is a `Abs. Energy` checkbox in the DFT options row, stored as `DFT_opts['EABS']`. It is passed to `Create_output`, `Create_grid_output`, `Create_pull_output`, and `Create_rot_output`, and controls whether `create_python_get_EABS` is called (which generates an extra bookkeeping Python script on the cluster).

**Fix:** Add a paragraph to §5 (DFT section) or create a new subsection explaining what absolute energy bookkeeping does and when to use it.

---

### 1.5 Button names differ between manual and UI (Minor)

| Manual label | Actual UI button label (`ANT.UI.py`) |
|---|---|
| "Scan tool" | **Grid Assistant** |
| "Pull tool" | **Pull Assistant** |
| "Rotation tool" | **Rotation Assistant** |

The manual uses "tool" throughout §6 (Advance functionalities), but the UI says "Assistant".

**Fix:** Either rename the buttons to match the manual or update the manual. Renaming the manual is easier and more consistent since §6 headings are "Pull tool", "Scan tool", "Rotation tool" — consider also updating section titles to match the UI labels.

---

### 1.6 Old filename in PyInstaller comment (Cosmetic)

**Code (`ANT.UI.py`, line 1):**
```python
#pyinstaller --windowed --icon="Icon.ico" --paths="." DFT_editor.py
```
The file is now `ANT.UI.py`, not `DFT_editor.py`.

**Fix:** Update the comment to `ANT.UI.py`.

---

### 1.7 Typo in "Specular Reflection" button (Cosmetic)

**Code (`ANT.UI.py`, line 340):**
```python
text = 'Specular Refelction'
```
"Refelction" should be "Reflection". This appears on screen in the UI.

**Fix:** Correct the spelling in the source.

---

### 1.8 Molecule basis set claim (Verify)

**Manual (§4, Basis sets):** "For atoms in the molecule, ANTUI will always use lanl2dz."

This needs verification against `_OUT.py` / `_DFT.py` `opt_gjf` and `create_gjf` to confirm no per-element lookup is done for molecule atoms. The current audit suggests this is correct, but verify before finalising the paper.

---

## 2. GitHub Readiness Checklist

### 2.1 Fix the LICENSE placeholder (Blocking)

`LICENSE` currently reads:
```
Copyright (c) 2026 [Andrés <FULL NAME / INSTITUTION>]
```
You must fill in your full legal name and institution before making the repo public. MIT license itself is fine — it is on SoftwareX's approved list.

---

### 2.2 Write a README.md (Blocking for GitHub)

A `README.md` at the repo root is the first thing anyone sees. It should cover:
- What ANT.UI is and what problem it solves (2–3 sentences)
- Prerequisites: Python ≥ 3.8, ANT.Gaussian installation
- Installation: `pip install -r requirements.txt`
- Quick start: how to launch, how to configure `USER_CONF.txt`
- Pointer to the manual PDF / paper DOI once published
- License badge

---

### 2.3 Add a .gitignore (Blocking for GitHub)

At minimum exclude:
```
__pycache__/
*.pyc
Outputs/
*.log
```
The `__pycache__` directory is already in the folder and must not be committed.

---

### 2.4 Remove stray files (Good hygiene)

- `SOC_params/Nuevo Documento de texto.txt` — an empty "New Text Document" left over from Windows. Delete it before the first commit.

---

### 2.5 Pin or clarify version number (Good hygiene)

`ANT.UI.py` defines `__version__ = '1.0.0'`. Ensure this matches the GitHub release tag you create when submitting to SoftwareX.

---

### 2.6 Create a GitHub release + Zenodo DOI (Required by SoftwareX)

After the repo is public:
1. Tag the commit: `git tag v1.0.0`
2. Create a GitHub Release from that tag.
3. Link Zenodo to the repo so it mints a DOI automatically on each release.
4. Include the Zenodo DOI in the SoftwareX manuscript metadata table.

---

### 2.7 Consider a `CITATION.cff` file

A machine-readable citation file makes it easy for users to cite the software correctly once the paper is published. Add one after the DOI is known.

---

## 3. SoftwareX Paper Requirements

Source: [Guide for authors](https://www.sciencedirect.com/journal/softwarex/publish/guide-for-authors)

| Requirement | Status |
|---|---|
| Open-source licence (approved list) | ✅ MIT — approved |
| Public GitHub repository | ❌ Not yet created |
| Max 3,000 words (excl. title, authors, refs, metadata tables) | ❌ Manual is longer than this; paper needs to be written fresh using SoftwareX template |
| Max 6 figures | ❌ Manual has 6 figures — fits, but one is a composite (Fig. 6) |
| Must use SoftwareX Word (.docx) template | ❌ Current manual is LaTeX; the paper must be submitted as .docx |
| Code DOI (Zenodo or similar) | ❌ Requires GitHub release first |
| APC | USD 1,560 (check if your institution has a read-and-publish deal with Elsevier) |

### 3.1 Paper structure (SoftwareX template)

The SoftwareX article format is not a traditional paper. It must include:

1. **Metadata table** — software name, version, repository URL, code DOI, language, dependencies, OS, licence, link to paper.
2. **Motivation and significance** — why this software fills a gap (educational value, reducing manual input effort, ANT.Gaussian accessibility).
3. **Software description** — architecture overview (4 Python modules), key UI components, input/output file formats.
4. **Illustrative examples** — reuse the benzene/fullerene scan and pull examples from the manual's §6.4.
5. **Impact** — who uses it, what research it enables.
6. **Conclusions**.
7. **Conflict of interest**.
8. **References** — cite ANT.Gaussian, Gaussian 09/16, the key NEGF-DFT papers.

The current LaTeX manual is an excellent source of content but needs to be condensed to ≤3,000 words and reformatted. The biggest rewriting work is §1–§3.

---

## 4. Suggested Order of Work

1. **Fix code first**
   - Correct the `Specular Refelction` typo (§1.7)
   - Update the PyInstaller comment filename (§1.6)
   - Decide on REAX: remove dead code or re-enable and test (§1.2)

2. **Prepare the repository**
   - Fill in LICENSE (§2.1)
   - Write README.md (§2.2)
   - Add .gitignore (§2.3)
   - Delete `SOC_params/Nuevo Documento de texto.txt` (§2.4)
   - Create the GitHub repo, push, tag v1.0.0, link Zenodo

3. **Update the manual / write the paper**
   - Fix L501 → L502 (§1.1)
   - Add ZMAX and U_ELECTRON_SCREENING to keyword list (§1.3)
   - Add EABS/Abs. Energy section (§1.4)
   - Align tool/assistant naming (§1.5)
   - Rewrite as SoftwareX .docx using their template (§3)

4. **Submit**
   - Article type: "Original Software Publication"
   - Attach repo URL + Zenodo DOI in metadata table
