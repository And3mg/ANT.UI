"""Helpers for the ANT.UI test suite.

Nothing here imports ANT.UI code: the checks are independent oracles for what
the generated files must contain (atom list, basis/ECP coverage, electron
count) and a stand-in for Gaussian so the generated link_opt.py can be run.
"""
import os
import re
import subprocess
import sys

import numpy as np

Z = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
     'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
     'K': 19, 'Ca': 20, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Br': 35,
     'Ag': 47, 'I': 53, 'Pt': 78, 'Au': 79}

# Core electrons replaced by the Hay-Wadt LANL2DZ ECP, by atomic number range.
LANL2DZ_CORE = [(1, 10, 0), (11, 30, 10), (31, 48, 28), (49, 57, 46), (72, 80, 60), (81, 83, 78)]


def lanl2dz_core(z):
    for lo, hi, core in LANL2DZ_CORE:
        if lo <= z <= hi:
            return core
    raise KeyError(z)


# ----------------------------------------------------------------- .gjf parsing
_CENTERS = re.compile(r'^\s*(\d+\s+)*\d+\s*$')
_ELEM_CENTER = re.compile(r'^\s*([A-Z][a-z]?)\s+0\s*$')


class Gjf:
    """Parsed Gaussian input as written by ANT.UI (gen basis + pseudo=read)."""

    def __init__(self, path):
        self.path = path
        self.text = open(path).read()
        lines = self.text.split('\n')
        self.route = next(l for l in lines if l.startswith('#'))
        i = next(k for k, l in enumerate(lines) if re.match(r'^\s*-?\d+\s*,\s*\d+\s*$', l))
        self.charge, self.mult = (int(v) for v in lines[i].split(','))

        atoms, k = [], i + 1
        while lines[k].strip() and not lines[k].startswith('Variables'):
            atoms.append(lines[k].split())
            k += 1
        self.variables = {}
        if lines[k].startswith('Variables'):
            k += 1
            while lines[k].strip():
                name, value = lines[k].split('=')
                self.variables[name.strip()] = float(value)
                k += 1
        self.symbols = [a[-3:] for a in atoms]
        self.elements = [a[0] for a in atoms]
        self.coords = np.array([[self.variables[t] if t in self.variables else float(t)
                                 for t in a[-3:]] for a in atoms])
        while not lines[k].strip():
            k += 1
        self.basis, k = self._parse_basis(lines, k)
        self.ecp = self._parse_ecp(lines, k)

    def _centers(self, line):
        """Atom indices (1-based) a basis/ECP block applies to."""
        m = _ELEM_CENTER.match(line)
        if m:
            return [i + 1 for i, e in enumerate(self.elements) if e == m.group(1)]
        if _CENTERS.match(line):
            return [int(x) for x in line.split()]
        raise ValueError(f'{self.path}: not a centre line: {line!r}')

    def _parse_basis(self, lines, k):
        basis = {}
        while lines[k].strip():
            centers = self._centers(lines[k])
            body = []
            k += 1
            while '****' not in lines[k]:
                body.append(lines[k])
                k += 1
            k += 1
            for c in centers:
                basis.setdefault(c, []).append(body)
        return basis, k

    def _parse_ecp(self, lines, k):
        ecp = {}
        n = len(lines)
        while k < n:
            if not lines[k].strip():
                k += 1
                continue
            centers = self._centers(lines[k])
            k += 1
            head = lines[k].split()
            k += 1
            if len(head) == 1:                     # keyword, e.g. lanl2dz
                entry = ('keyword', head[0].lower())
            else:                                  # NAME lmax ncore + lmax+1 blocks
                lmax, core = int(head[1]), int(head[2])
                for _ in range(lmax + 1):
                    k += 1                         # block title
                    k += 1 + int(lines[k].split()[0])
                entry = ('explicit', core)
            for c in centers:
                ecp.setdefault(c, []).append(entry)
        return ecp


def basis_is_all_electron(body):
    """True when the basis carries core functions (tightest exponent > 500)."""
    if len(body) == 1 and body[0].strip().lower() == 'lanl2dz':
        return False
    exps = []
    for line in body:
        parts = line.split()
        if len(parts) >= 2 and re.match(r'^[SPDFL]+$', parts[0]):
            continue
        try:
            exps.append(float(parts[0].replace('D', 'E')))
        except (ValueError, IndexError):
            pass
    return bool(exps) and max(exps) > 500


def basis_ecp_problems(gjf):
    """List of human-readable problems with basis/ECP coverage of every atom."""
    problems = []
    for idx, elem in enumerate(gjf.elements, 1):
        b = gjf.basis.get(idx, [])
        e = gjf.ecp.get(idx, [])
        if len(b) != 1:
            problems.append(f'atom {idx} ({elem}): {len(b)} basis assignments')
            continue
        if len(e) > 1:
            problems.append(f'atom {idx} ({elem}): {len(e)} ECP assignments')
        needs_ecp = Z[elem] > 18 and not basis_is_all_electron(b[0])
        if needs_ecp and not e:
            problems.append(f'atom {idx} ({elem}): basis needs an ECP but none is assigned')
    return problems


def electron_count(gjf):
    """Electrons of the neutral system after removing ECP cores."""
    total = 0
    for idx, elem in enumerate(gjf.elements, 1):
        core = 0
        for kind, value in gjf.ecp.get(idx, []):
            core = lanl2dz_core(Z[elem]) if kind == 'keyword' else value
        total += Z[elem] - core
    return total


def spin_consistent(gjf):
    electrons = electron_count(gjf) - gjf.charge
    return electrons % 2 == (gjf.mult - 1) % 2


# --------------------------------------------------------------- fake Gaussian
def fake_gaussian_log(gjf_path):
    """Write the .log Gaussian would leave if the optimisation did not move any atom.

    Only the sections link_opt.py parses are written: the 'Input orientation'
    table and the final 'Optimized Parameters' table of z-matrix variables.
    """
    g = Gjf(gjf_path)
    out = [' Input orientation:', ' ' + '-' * 69,
           ' Center     Atomic      Atomic             Coordinates (Angstroms)',
           ' Number     Number       Type             X           Y           Z',
           ' ' + '-' * 69]
    for k, (e, r) in enumerate(zip(g.elements, g.coords), 1):
        out.append(f' {k:6d} {Z[e]:10d} {0:11d}  {r[0]:12.6f}{r[1]:12.6f}{r[2]:12.6f}')
    out += [' ' + '-' * 69, ' Rotational constants (GHZ):      0.0000000      0.0000000      0.0000000', '',
            '                       !   Optimized Parameters   !',
            '                       ! (Angstroms and Degrees)  !',
            ' ! Name  Definition              Value          Derivative Info.                !']
    for name, value in g.variables.items():
        out.append(f' ! {name:10s} {value:14.6f}   -DE/DX =    0.0                              !')
    out.append(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
    with open(gjf_path[:-4] + '.log', 'w') as f:
        f.write('\n'.join(out) + '\n')


def run_link_opt(folder):
    r = subprocess.run([sys.executable, 'link_opt.py'], cwd=folder, capture_output=True, text=True)
    assert r.returncode == 0, f'link_opt.py failed in {folder}:\n{r.stderr}'


def python_syntax_error(path):
    import ast
    try:
        ast.parse(open(path).read(), filename=path)
        return None
    except SyntaxError as e:
        return f'{os.path.basename(path)} line {e.lineno}: {e.text!r}'
