import numpy as np
from os import path

#%%
from _DICS import *

p_line = '#p HSEH1PBE/gen pseudo=read pop=min Scf=(conver=5,maxcycle=1000,dsymm,nodamp) nosymm'
Flag_opt = False
Flag_EABS = False
Flag_opt_in_ANT = False
bethe_layers = []

def begin_export():
    '''Reset the module state shared by the writers, so no setting leaks from one export to the next.'''
    global Flag_opt, Flag_EABS, Flag_opt_in_ANT, bethe_layers
    Flag_opt = False
    Flag_EABS = False
    Flag_opt_in_ANT = False
    bethe_layers = []

def read_template(name):
    with open(f'PythonRutines/{name}', 'r') as read:
        return read.read()

#Launcher SA
def create_unisa(name, Dir):
    global Flag_opt, Flag_EABS
    h = open(Dir+'/launchANT.sh','w', newline='\n')
    h.write(UserInfo['JOBHEAD'])
    bar = '/'
    h.write(f"{UserInfo['G09call']} input{Dir.split(bar)[-1]}.gjf")
    if Flag_opt:
        h.write('\n')
        h.write(f"{UserInfo['PYTHONcall']} link_opt.py \n")
        if Flag_opt_in_ANT:
            h.write(f"find input{Dir.split(bar)[-1]}.ant | xargs sed -i 's/EW1 = 3.0/EW1 = -3.0/g' \n")
        h.write(f"{UserInfo['G09call']} input{Dir.split(bar)[-1]}_opt.gjf")
    if Flag_EABS:
        h.write('\n')
        if not Flag_opt:
            h.write(f"{UserInfo['PYTHONcall']} link_opt.py \n")
            if Flag_opt_in_ANT:
                h.write(f"find input{Dir.split(bar)[-1]}.ant | xargs sed -i 's/EW1 = 3.0/EW1 = -3.0/g' \n")
        h.write(f"cp input{Dir.split(bar)[-1]}.ant input{Dir.split(bar)[-1]}_EABS.ant \n")
        h.write(f"{UserInfo['G09call']} input{Dir.split(bar)[-1]}_elec_opt.gjf \n")
        h.write(f"{UserInfo['G09call']} input{Dir.split(bar)[-1]}_mol_opt.gjf \n")
        create_python_get_EABS(Dir, name)
        h.write(f"{UserInfo['PYTHONcall']} Get_EABS.py \n")
    h.close()
    return

#Multiexe functions
def create_unisa_multi(N, Dir, Opt_chain = 'none'):
    h = open(Dir+'/launchANT.sh','w', newline='\n')
    h.write(UserInfo['JOBHEAD'])
    h.write(f'name={path.basename(path.normpath(Dir))}\n')
    h.write(f'printf -v first "%0{len(str(N))}d" 1 \n')
    h.write('# Launch program \n')
    h.write(f'for i in $(seq -f "%0{len(str(N))}g" 1 {N})\n')
    h.write('do \n')
    h.write('  if [[ $i != $first ]] \n')
    h.write('  then \n')
    h.write('    cp P.* ../$name\\_$i/P.input$name\\_$i.dat \n')
    h.write('    cd ../$name\\_$i \n')
    if Flag_opt_in_ANT and Opt_chain == 'first step':
        h.write(f"    find input$name\\_$i.ant | xargs sed -i 's/EW1 = 3.0/EW1 = -3.0/g' \n")
    h.write(f"    {UserInfo['G09call']} input$name\\_$i.gjf \n")
    if Opt_chain == 'every step':
        h.write(f"    {UserInfo['PYTHONcall']} link_opt.py \n")
        if Flag_opt_in_ANT:
            h.write(f"    find input$name\\_$i.ant | xargs sed -i 's/EW1 = 3.0/EW1 = -3.0/g' \n")
        h.write(f"    {UserInfo['G09call']} input$name\\_$i\\_opt.gjf \n")
    h.write('  else \n')
    h.write('    cd $name\\_$i \n')
    h.write(f"    {UserInfo['G09call']} input$name\\_$i.gjf \n")
    if Opt_chain != 'none':
        h.write(f"    {UserInfo['PYTHONcall']} link_opt.py \n")
        if Flag_opt_in_ANT:
            h.write(f"    find input$name\\_$i.ant | xargs sed -i 's/EW1 = 3.0/EW1 = -3.0/g' \n")
        h.write(f"    {UserInfo['G09call']} input$name\\_$i\\_opt.gjf \n")

    h.write('  fi\n')
    h.write('done \n')
    h.write('cd .. \n')
    h.write('mkdir T; \n')
    h.write('mkdir xyz; \n')
    h.write('mkdir Tlog; \n')
    h.write('find . -type f -name "*opt.log" -exec cp {} Tlog/ \\; \n')
    h.write('find . -type f -name "T.*" -exec cp {} T/ \\; \n')
    h.write('find . -type f -name "input*.xyz" -exec cp {} xyz/ \\; \n')
    h.close()
    return

def create_ant(name, Dir, bl1, bl2, DFTU, DOSAT = None):
    global Flag_opt_in_ANT
    out_name = Dir + '/input' + name+'.ant'
    h = open(out_name, "w")
    h.write(UserInfo['ANTpredet'])
    if DOSAT != None:
        h.write(f"LDOS_BEG = {DOSAT}  \n")
        h.write(f"LDOS_END = {DOSAT}  \n")
    h.write(f"BLPAR1 = {BLPAR[bl1]} \n")
    h.write(f"BLPAR2 = {BLPAR[bl2]} \n")
    if Flag_opt_in_ANT:
        h.write('EW1 = 3.0 \n')
    if DFTU:
        h.write('DFT+U \n')
        h.write('MOLBLOCKS \n')
        h.write('1 \n')
        h.write(f'{DFTU[0]},{DFTU[1]} \n')
        h.write(f'  UPLUS = {UserInfo["U_ELECTRON_SCREENING"]*float(DFTU[-1]):.2f} \n')
    h.close()
    return

def optimization_protocol(opt_data, Semi_infinite_layer, z_order, counts):
    '''Optimisation mode of every atom, by region: bottom electrode, molecule, top electrode.
    Bethe layers relaxed along z get 'zB'/'zT' so each layer is written as one rigid plane.'''
    Nbtm, Nmol, Ntop = counts
    opt_resume = np.zeros(Nbtm+Nmol+Ntop, dtype = '<U20')
    for i in range(Nbtm+Nmol+Ntop):
        top = i >= Nbtm+Nmol
        if Nbtm <= i < Nbtm+Nmol:
            opt_resume[i] = opt_data['Molecule'] if opt_data['Molecule'] in ['all', 'fix', 'z'] else 'all'
        elif Semi_infinite_layer[i]:
            mode = opt_data['Top_Bethe'] if top else opt_data['Bottom_Bethe']
            opt_resume[i] = mode + ('T' if top else 'B') if mode == 'z' else mode
        else:
            opt_resume[i] = opt_data['Top'] if top else opt_data['Bottom']
    if opt_data['Molecule'] not in ['all', 'fix', 'z']:
        # Molecule_z_atoms index the molecule file; z_order maps output position -> input atom
        for x in opt_data['Molecule_z_atoms']:
            opt_resume[np.where(z_order == Nbtm + x)[0][0]] = 'z'
    return opt_resume

def bethe_layer_offsets(layers):
    '''Gaussian09 can only tie each Bethe plane to its own z variable, so the planes relax independently.
    For each side with several planes, record every plane's original offset from the innermost one
    (the layer touching the active electrode); link_opt.py restores these offsets after the optimisation.'''
    sides = []
    for side in ['zB', 'zT']:
        planes = {var: layer for var, layer in layers.items() if var.startswith(side)}
        if len(planes) < 2:
            continue
        innermost = max if side == 'zB' else min
        reference = innermost(planes, key = lambda var: planes[var][0])
        sides.append({'reference': reference,
                      'layers': {var: [round(z0 - planes[reference][0], 6), atoms] for var, (z0, atoms) in planes.items()}})
    return sides

def opt_gjf(Poss, Elem, Opt_resume, Opt_data, DFT_opts):
    global p_line, Flag_opt, bethe_layers
    #Juega con esto para añadir HOMO/LUMO representation, cambiar funcional...
    Flag_opt = False
    P_line = ''
    P_line += f"#p {Shell[DFT_opts['Shell']]}{Functionals[DFT_opts['Functional']]}/gen pseudo=read pop=min"
    P_line += " Scf=(conver=6,maxcycle=1000,dsymm,nodamp) nosymm IOP(3/124=3)"
    #P_line += " Scf=(conver=6,maxcycle=500,dsymm,nodamp) nosymm EmpiricalDispersion=GD3BJ"
    p_line = P_line
    character = ''
    Variable_lines = ''
    if (Opt_resume != 'fix').any():
        P_line += f" opt(z-matrix,maxcycles={Opt_data[r'Opt_steps']},loose)"
        Variable_lines = 'Variables: \n'
        character = '0'
        Flag_opt = True
    P_line += ' \n'

    Poss_lines = ''
    k_opt = 1
    layer_vars = {}
    layers = {}   # Bethe plane variable -> [original z, [atom indices]]
    for idx, (elem, poss, opt) in enumerate(zip(Elem, Poss, Opt_resume)):
        if opt == 'fix':
            Poss_lines += f'{elem} {character}  {poss[0]:.6f}    {poss[1]:.6f}    {poss[2]:.6f}\n'
        elif opt == 'all':
            Poss_lines += f'{elem} {character}  x{k_opt}    y{k_opt}    z{k_opt}\n'
            Variable_lines += f'x{k_opt} = {poss[0]:.6f} \n'
            Variable_lines += f'y{k_opt} = {poss[1]:.6f} \n'
            Variable_lines += f'z{k_opt} = {poss[2]:.6f} \n'
            k_opt +=1
        elif opt == 'z':
            Poss_lines += f'{elem} {character}  {poss[0]:.6f}    {poss[1]:.6f}    z{k_opt}\n'
            Variable_lines += f'z{k_opt} = {poss[2]:.6f} \n'
            k_opt += 1
        elif opt in ['zB', 'zT']:
            layer = (opt, f'{poss[2]:.4f}')
            if layer not in layer_vars:
                layer_vars[layer] = f'{opt}{sum(side == opt for side, _ in layer_vars)+1}'
                Variable_lines += f'{layer_vars[layer]} = {poss[2]:.6f} \n'
            Poss_lines += f'{elem} {character}  {poss[0]:.6f}    {poss[1]:.6f}    {layer_vars[layer]}\n'
            layers.setdefault(layer_vars[layer], [float(poss[2]), []])[1].append(idx)
        else:
            raise ValueError(f'Unknown optimisation mode {opt!r} for {elem}')
    Variable_lines += '\n'
    bethe_layers = bethe_layer_offsets(layers)
    return P_line, Poss_lines, Variable_lines

def count_electrons(Elem, Semi, active_ecp = None):
    '''Electrons of the neutral system once the ECP cores are removed: lanl2dz for molecule atoms,
    the minimal-basis ECP for Bethe anchor atoms (Semi), active_ecp for the other electrode atoms.'''
    if active_ecp is None:
        active_ecp = DFT_pseudopotentials
    electrons = 0
    for i, e in enumerate(Elem):
        if e not in Electrode_elem:
            core = lanl2dz_core(e)
        elif Semi[i]:
            core = ecp_core(DFT_min_pseudopotentials[e], e)
        else:
            core = ecp_core(active_ecp.get(e, ''), e)
        electrons += atomic_number[e] - core
    return electrons

def charge_multiplicity(electrons, DFT_opts):
    '''Auto: closed shell removes one electron from an odd count (charge +1, singlet);
    open shell keeps the system neutral (doublet if odd). DFT_opts Charge/Multiplicity override.'''
    if electrons % 2 == 0:
        charge, mult = 0, 1
    elif DFT_opts['Shell'] == 'Close Shell':
        charge, mult = 1, 1
    else:
        charge, mult = 0, 2
    if DFT_opts.get('Charge') is not None:
        charge = int(DFT_opts['Charge'])
    if DFT_opts.get('Multiplicity') is not None:
        mult = int(DFT_opts['Multiplicity'])
    return charge, mult

def create_gjf(name, Dir, Poss, Elem, Semi, Opt_resume, Opt_data, DFT_opts, active_basis = None, active_ecp = None):
    '''Semi marks the Bethe anchor atoms (minimal basis). active_basis/active_ecp give the basis and ECP
    of the other electrode atoms (default: DFT_basis_pseudopot, SOC files for SOC inputs).'''
    global Flag_opt, Flag_opt_in_ANT
    if active_basis is None:
        active_basis, active_ecp = DFT_basis, DFT_pseudopotentials
    P_line, Poss_lines, Variable_lines = opt_gjf(Poss, Elem, Opt_resume, Opt_data, DFT_opts)
    Flag_opt_in_ANT = Opt_data['Opt_in_ANT']

    Metals = [x for x in np.unique(Elem) if x in Electrode_elem]
    Others = [x for x in np.unique(Elem) if x not in Electrode_elem]
    active = {met: [i+1 for i, e in enumerate(Elem) if e == met and not Semi[i]] for met in Metals}
    anchor = {met: [i+1 for i, e in enumerate(Elem) if e == met and Semi[i]] for met in Metals}

    charge, mult = charge_multiplicity(count_electrons(Elem, Semi, active_ecp), DFT_opts)

    with open(f'{Dir}/input{name}.gjf', "w") as g:
        g.write(UserInfo['GAUSSIANHEAD'])
        if not Flag_opt or Flag_opt_in_ANT:
            g.write(f"%Subst L502 {UserInfo['ANTDIRcall']} \n")
            g.write(f"%Subst L101 {UserInfo['ANTDIRcall']} \n")
        g.write(P_line)
        g.write("\n")
        g.write(f'input{name}')
        g.write("\n")
        g.write("\n")
        g.write(f'{charge},{mult} \n')

        g.write(Poss_lines)
        g.write(Variable_lines)

        for elem in Others:
            g.write(f'{elem} 0 \n')
            g.write('lanl2dz \n')
            g.write('**** \n')
        for met in Metals:
            if active[met]:
                g.write(' '.join(map(str, active[met])) + ' \n')
                g.write(active_basis[met].rstrip('\n') + '\n')
                g.write('**** \n')
        for met in Metals:
            if anchor[met]:
                g.write(' '.join(map(str, anchor[met])) + ' \n')
                g.write(DFT_min_basis[met].rstrip('\n') + '\n')
                g.write('**** \n')
        g.write('\n')

        for elem in Others:
            g.write(f'{elem} 0 \n')
            g.write('lanl2dz \n')
        for met in Metals:  # an all-electron basis (e.g. the Cu SOC basis) has no ECP block
            if active[met] and active_ecp.get(met, '').strip():
                g.write(' '.join(map(str, active[met])) + ' \n')
                g.write(active_ecp[met].rstrip('\n') + '\n')
        for met in Metals:
            if anchor[met] and DFT_min_pseudopotentials[met].strip():
                g.write(' '.join(map(str, anchor[met])) + ' \n')
                g.write(DFT_min_pseudopotentials[met].rstrip('\n') + '\n')
        g.write(4*'\n ')
    return

def create_python_conector(Dir, name, FlagEABS, text):
    global Flag_opt, Flag_EABS, Flag_opt_in_ANT, bethe_layers
    Flag_EABS = FlagEABS
    if not Flag_EABS and not Flag_opt:
        return
    global p_line
    common = read_template('CommonLink.txt')
    main = ''
    if Flag_opt:
        main += read_template('Main.txt')
    if Flag_EABS:
        main += read_template('Main_EABS.txt')

    with open(f'{Dir}/link_opt.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        out.write(f"optgjf_path = 'input{name}.gjf' \n")
        out.write(f"log_file_path = 'input{name}.log' \n")
        out.write(f"P_line = f'{p_line}")
        out.write(" {intro}' \n")
        out.write(f"Bethe_layers = {bethe_layers!r} \n")
        out.write(f"{text} \n")
        if not Flag_opt or Flag_opt_in_ANT:
            out.write(f"Ant1 = ''\n")
            out.write(f"Ant2 = ''\n")
        else:
            out.write(f"Ant1 = f'%Subst L502 {UserInfo['ANTDIRcall']}")
            out.write(" {intro}' \n")
            out.write(f"Ant2 = f'%Subst L101 {UserInfo['ANTDIRcall']}")
            out.write(" {intro}' \n")
        out.write(common)
        out.write(main)

def create_python_get_EABS(Dir, name):
    global Flag_opt, Flag_EABS
    main = read_template('Get_EABS.txt')
    with open(f'{Dir}/Get_EABS.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        if Flag_opt:
            out.write(f"Ee_file = 'input{name}_opt.log'")
        else:
            out.write(f"Ee_file = 'input{name}.log'")
        out.write('\n')
        out.write(f"Eelec_file = 'input{name}_elec_opt.log' \n")
        out.write(f"Emol_file = 'input{name}_mol_opt.log' \n")
        out.write(main)

def create_python_OptFirstStep(Dir, name, Folder_name, N_folder, text_func):
    global Flag_opt, Flag_opt_in_ANT
    global p_line
    common = read_template('CommonLink.txt')
    main = read_template('Main_OptFirstStep.txt')

    with open(f'{Dir}/link_opt.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        out.write(f"optgjf_path = 'input{name}.gjf' \n")
        out.write(f"log_file_path = 'input{name}.log' \n")
        out.write(f"P_line = f'{p_line}")
        out.write(" {intro}' \n")
        out.write(f"Bethe_layers = {bethe_layers!r} \n")
        if not Flag_opt_in_ANT:
            out.write(f"Ant1 = f'%Subst L502 {UserInfo['ANTDIRcall']}")
            out.write(" {intro}' \n")
            out.write(f"Ant2 = f'%Subst L101 {UserInfo['ANTDIRcall']}")
            out.write(" {intro}' \n")
        else:
            out.write(f"Ant1 = '' \n")
            out.write(f"Ant2 = '' \n")
        out.write('\n')
        out.write(f"N_folder = '{N_folder}' \n")
        out.write(f"Folder_name = '{Folder_name}' \n")
        out.write(text_func)
        out.write('\n')
        out.write(common)
        out.write('\n')
        out.write(main)



def create_python_OptAllSteps(Dir, name, nxt, text_func):
    global Flag_opt, Flag_opt_in_ANT
    global p_line
    common = read_template('CommonLink.txt')
    main = read_template('Main_OptAllSteps.txt')

    with open(f'{Dir}/link_opt.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        out.write(f"optgjf_path = 'input{name}.gjf' \n")
        out.write(f"log_file_path = 'input{name}.log' \n")
        out.write(f"P_line = f'{p_line}")
        out.write(" {intro}' \n")
        out.write(f"Bethe_layers = {bethe_layers!r} \n")
        if not Flag_opt_in_ANT:
            out.write(f"Ant1 = f'%Subst L502 {UserInfo['ANTDIRcall']}")
            out.write(" {intro}' \n")
            out.write(f"Ant2 = f'%Subst L101 {UserInfo['ANTDIRcall']}")
            out.write(" {intro}' \n")
        else:
            out.write(f"Ant1 = '' \n")
            out.write(f"Ant2 = '' \n")
        out.write('\n')
        out.write(f"Next_folder = '{nxt}' \n")
        out.write(text_func)
        out.write('\n')
        out.write(common)
        out.write('\n')
        out.write(main)
