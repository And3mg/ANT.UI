'''Plain-text description of an export: shown before the files are written and saved with them
(ANTUI_summary.txt), with warnings for settings that would give a broken or surprising calculation.'''
import numpy as np

from _DICS import *
from _DFT import optimization_protocol, charge_multiplicity, count_electrons, gas_states
from _OUT import assemble
from _SOC import SOC_PARAMS, SOC_basis, SOC_pseudopot

KIND_TITLES = {'xyz': 'XYZ geometry', 'ANT': 'ANT input', 'SOC': 'ANT input with spin-orbit coupling',
               'pull': 'Pull Assistant sequence', 'grid': 'Grid Assistant sequence',
               'rotation': 'Rotation Assistant sequence', 'getU': 'Hubbard U of the isolated molecule'}
SEQUENCES = ['pull', 'grid', 'rotation']
CONTACT_WARNING = 1.0   # Angstrom: closer atoms almost certainly overlap

def formula(elements):
    '''Hill order: C, H, then alphabetical.'''
    count = {e: int(np.sum(np.asarray(elements) == e)) for e in np.unique(elements)}
    order = [e for e in ('C', 'H') if e in count] + sorted(e for e in count if e not in ('C', 'H'))
    return ''.join(e + (str(count[e]) if count[e] > 1 else '') for e in order)

def min_distance(a, b):
    return float(np.min(np.linalg.norm(a[:, None, :] - b[None, :, :], axis = 2)))

def contact_lines(Poss, counts, label = 'Closest contacts'):
    Nb, Nm, _ = counts
    bottom, molecule, top = Poss[:Nb], Poss[Nb:Nb+Nm], Poss[Nb+Nm:]
    d = {'molecule-bottom': min_distance(molecule, bottom), 'molecule-top': min_distance(molecule, top),
         'bottom-top': min_distance(bottom, top)}
    line = f'  {label}: ' + ', '.join(f'{k} {v:.2f} A' for k, v in d.items())
    warnings = [f'{label}: {k} atoms are {v:.2f} A apart (< {CONTACT_WARNING} A): the parts overlap.'
                for k, v in d.items() if v < CONTACT_WARNING]
    return line, warnings

def ecp_label(text, elem, source):
    if not text.strip():
        return 'none (all-electron basis)'
    return f'{source} ({ecp_core(text, elem)} core electrons)'

def describe_get_u(molecule, DFT_opts):
    '''Summary of the gas-phase calculations behind a molecular Hubbard U.'''
    functional = DFT_opts['Functional']
    route = Functionals[functional]
    electrons = count_electrons(molecule.elements, np.zeros(molecule.N, dtype = 'bool'))
    lines = [f'Calculation: {KIND_TITLES["getU"]}',
             f'Molecule: {molecule.name} ({formula(molecule.elements)}, {molecule.N} atoms),'
             f' coordinates as in Molecules/{molecule.name}.xyz',
             f'Functional: {functional} (Gaussian route {route}); basis lanl2dz, ECP lanl2dz',
             f'Electrons after ECP cores: {electrons}']
    names = {'N': 'neutral, N electrons', 'cation': 'cation, N-1', 'anion': 'anion, N+1'}
    for tag, charge, mult, shell in gas_states(electrons, DFT_opts):
        lines.append(f'  {molecule.name}_{tag}.gjf: {names[tag]}; charge {charge}, multiplicity {mult} ({shell}{route})')
    lines += ['U = (IP - EA) - (e_LUMO - e_HOMO): the gas-phase correction that ANT splits as -U/2 on occupied'
              ' and +U/2 on empty molecular levels.',
              f'Get_U.py writes U.{molecule.name}.dat and {molecule.name}.xyz with "{functional} U" in its comment line;'
              ' copy that file to Molecules/ to use it.']
    if functional in molecule.U:
        lines.append(f'  The molecule file already has U = {molecule.U[functional]} eV for {functional}; the new file replaces it.')
    lines.append(f"In the junction ANT uses UPLUS = U x U_ELECTRON_SCREENING (now {UserInfo['U_ELECTRON_SCREENING']}),"
                 ' which the researcher chooses for the molecule and contact geometry.')
    return lines, site_config_warnings()

def describe_export(objs, kind, opt_data = None, DFT_opts = None, params = None):
    '''Return (lines, warnings) describing what an export of this kind will write.'''
    if kind == 'getU':
        return describe_get_u(objs[1], DFT_opts)
    params = params or {}
    lines, warnings = [], []
    bottom, molecule, top = objs
    Poss, Elem, Semi, Extremes, z_order, counts = assemble(objs)
    Nb, Nm, Nt = counts
    soc = kind == 'SOC' or (kind == 'rotation' and params.get('soc', False))

    lines.append(f'Calculation: {KIND_TITLES[kind]}')
    lines.append(f'Atoms: {len(Elem)} = bottom electrode {Nb} ({bottom.elements[0]}, {bottom.name})'
                 f' + molecule {Nm} ({molecule.name}, {formula(molecule.elements)})'
                 f' + top electrode {Nt} ({top.elements[0]}, {top.name})')
    lines.append(f'  Bethe region (outer layers): bottom {int(Semi[:Nb].sum())} atoms, top {int(Semi[Nb+Nm:].sum())};'
                 f' anchor layer with minimal basis: bottom {int(Extremes[:Nb].sum())}, top {int(Extremes[Nb+Nm:].sum())}')
    line, w = contact_lines(Poss, counts)
    lines.append(line)
    warnings += w
    if kind == 'xyz':
        return lines, warnings

    # Electronic model
    shell = DFT_opts['Shell']
    lines.append(f"Functional: {DFT_opts['Functional']} (Gaussian route {Shell[shell]}{Functionals[DFT_opts['Functional']]}), {shell.lower()}")
    active_ecp = SOC_pseudopot if soc else DFT_pseudopotentials
    try:
        electrons = count_electrons(Elem, Extremes, active_ecp)
    except KeyError as e:
        electrons = None
        warnings.append(f'No ECP core size known for {e}: cannot count electrons (is lanl2dz defined for it?).')
    if electrons is not None:
        charge, mult = charge_multiplicity(electrons, DFT_opts)
        user_set = DFT_opts.get('Charge') is not None or DFT_opts.get('Multiplicity') is not None
        lines.append(f"Electrons after ECP cores: {electrons}; charge {charge}, multiplicity {mult}"
                     f" ({'set by the user' if user_set else 'automatic'})")
        if not user_set and electrons % 2 and shell == 'Close Shell':
            lines.append('  Odd electron count: one electron removed (charge +1) to keep a closed-shell singlet.')
        if (electrons - charge) % 2 != (mult - 1) % 2:
            warnings.append(f'Charge {charge} and multiplicity {mult} are impossible with {electrons} electrons: Gaussian will stop.')
        if mult > 1 and shell == 'Close Shell':
            warnings.append(f'Multiplicity {mult} needs an open-shell calculation, but Close Shell is selected.')

    # Basis sets and ECPs, per element and region
    lines.append('Basis set / ECP:')
    others = [e for e in np.unique(Elem) if e not in Electrode_elem]
    metals = [e for e in np.unique(Elem) if e in Electrode_elem]
    if others:
        lines.append(f"  {', '.join(others)}: lanl2dz / lanl2dz")
    for met in metals:
        n_active = int(np.sum((Elem == met) & ~Extremes))
        n_anchor = int(np.sum((Elem == met) & Extremes))
        if n_active:
            if soc:
                if met not in SOC_basis:
                    warnings.append(f'No SOC basis for {met} (SOC_params/{met}_SOC_basis.gbs).')
                    continue
                basis = f'SOC_params/{met}_SOC_basis.gbs'
                ecp = ecp_label(SOC_pseudopot.get(met, ''), met, f'SOC_params/{met}_SOC_pseudopot.gbs')
            else:
                basis = DFT_basis_source[met]
                ecp = ecp_label(DFT_pseudopotentials[met], met, 'lanl2dz' if basis.startswith('lanl2dz') else basis)
            lines.append(f'  {met} electrode ({n_active} atoms): basis {basis}; ECP {ecp}')
        if n_anchor:
            lines.append(f'  {met} Bethe anchor ({n_anchor} atoms): {DFT_min_basis_source[met]};'
                         f' ECP {ecp_label(DFT_min_pseudopotentials[met], met, "same file")}')
    lines.append(f'  Bethe lattice parameters: BLPAR1 = {BLPAR[Elem[0]]}, BLPAR2 = {BLPAR[Elem[-1]]}')

    if soc:
        n_soc = int(np.sum([(e in Electrode_elem) and not s for e, s in zip(Elem, Semi)]))
        lines.append(f'SOC: SOCFACATOM on {n_soc} electrode atoms outside the Bethe region')
        for met in metals:
            if met in SOC_PARAMS:
                lines.append(f"  {met} SOC parameters: {', '.join(SOC_PARAMS[met])}")
            else:
                warnings.append(f'No SOC parameters for {met}: the SOC input cannot be written.')

    if DFT_opts.get('DFTU'):
        U = float(molecule.U[DFT_opts['Functional']])
        s = UserInfo['U_ELECTRON_SCREENING']
        lines.append(f'DFT+U: molecule atoms {Nb+1}-{Nb+Nm}, U = {U} eV x screening {s} = {U*s:.2f} eV')

    # Optimisation
    chains = opt_data['Chains']
    if kind in SEQUENCES and chains == 'none':
        modes = np.full(len(Elem), 'fix')
    else:
        modes = optimization_protocol(opt_data, Semi, z_order, counts)
    if (modes != 'fix').any():
        mol_mode = opt_data['Molecule']
        if mol_mode not in ['fix', 'all', 'z']:
            mol_mode = f"all, z only for atoms {opt_data['Molecule_z_atoms']}"
        lines.append(f"Optimisation: Gaussian, at most {opt_data['Opt_steps']} steps"
                     f"{', inside ANT' if opt_data['Opt_in_ANT'] else ''}")
        lines.append(f"  bottom Bethe {opt_data['Bottom_Bethe']}, bottom electrode {opt_data['Bottom']}, molecule {mol_mode},"
                     f" top electrode {opt_data['Top']}, top Bethe {opt_data['Top_Bethe']}")
        if 'z' in (opt_data['Bottom_Bethe'], opt_data['Top_Bethe']):
            lines.append('  Bethe planes relax along z; link_opt.py then restores their original spacing from the innermost plane.')
        if opt_data['Opt_steps'] < 1:
            warnings.append('Number of optimisation steps is below 1.')
    else:
        lines.append('Optimisation: none (single-point calculations)')
        if kind in SEQUENCES and chains != 'none':
            warnings.append(f'Chains is "{chains}" but every region is fixed: nothing will be optimised.')
    if kind in SEQUENCES:
        chain_text = {'first step': ' (only the first folder is relaxed in Gaussian; its link_opt.py writes the inputs'
                                    ' of the other folders from the relaxed geometry, which are not re-optimised)',
                      'every step': ' (every folder is relaxed, starting from the relaxed geometry of the previous one)'}
        lines.append(f'  Chain: {chains}{chain_text.get(chains, "") if (modes != "fix").any() else ""}')

    # Sequence settings
    if kind == 'pull':
        dz, nz, both = params['dz'], params['nz'], params['both']
        lines.append(f"Pull: {nz} configurations; the gap grows by {dz} A per step"
                     f" ({'each electrode moves dz/2' if both else 'only the top electrode moves'});"
                     f" total change {dz*(nz-1):.3f} A")
        last = Poss.copy()
        last[Nb+Nm:, 2] += dz*(nz-1)/2 if both else dz*(nz-1)
        if both:
            last[:Nb, 2] -= dz*(nz-1)/2
        line, w = contact_lines(last, counts, 'Last step contacts')
        lines.append(line)
        warnings += w
    elif kind == 'grid':
        (x0, x1, nx), (y0, y1, ny), z = params['grid']
        lines.append(f'Grid: {nx} x {ny} = {nx*ny} tip positions, x {x0} to {x1} A, y {y0} to {y1} A,'
                     f' tip centre at z = {float(z):.3f} A (LDOS on atom {len(Elem)-Nt+1})')
    elif kind == 'rotation':
        angles = [float(a) for a in params['angles']]
        parts = [name for name, flag in zip(['bottom electrode', 'molecule', 'top electrode'], params['components']) if flag]
        lines.append(f"Rotation: {len(angles)} configurations at {', '.join(f'{a:g}' for a in angles)} degrees"
                     f" about z, relative to the current orientation; rotating: {', '.join(parts) or 'nothing'}")
        if not parts:
            warnings.append('Rotation Assistant: no component selected, every configuration will be identical.')

    lines.append(f"Site configuration (USER_CONF.txt): Gaussian call '{UserInfo['G09call']}', "
                 f"Python call '{UserInfo['PYTHONcall']}', ANT links '{UserInfo['ANTDIRcall']}'")
    warnings += site_config_warnings()
    return lines, warnings

def summary_text(lines, warnings, destination, version = ''):
    text = [f'ANT.UI {version} export summary'.replace('  ', ' '), f'Destination: {destination}', ''] + lines + ['']
    if warnings:
        text += ['Warnings:'] + [f'  ! {w}' for w in warnings]
    else:
        text += ['Warnings: none']
    return '\n'.join(text) + '\n'
