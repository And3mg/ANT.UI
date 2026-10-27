from _DICS import *
from _DFT import *
from _SOC import *
import ast
import numpy as np
from os import listdir, mkdir, walk, path


def create_xyz(Dir, name, Poss, Elem, N):
    out = open(f'{Dir}/{name}.xyz','w')
    out.write(f'{N}\n')
    out.write('\n')
    for i in range(N):
        out.write(f'{Elem[i]}  {Poss[i][0]}    {Poss[i][1]}    {Poss[i][2]} \n')
    out.close()
    return

def create_xyz_film(Dir, name, POSS, ELEM, N):
    out = open(f'{Dir}/{name}.xyz','w')
    t = 0
    for Poss, Elem in zip(POSS, ELEM):
        out.write(f'{N}\n')
        out.write(f'Atoms. Timestep: {t}\n')
        for i in range(N):
            out.write(f'{Elem[i]}  {Poss[i][0]}    {Poss[i][1]}    {Poss[i][2]} \n')
        t+=1
    out.close()
    return


def create_reax(Dir, name, poss, Elem, N, types):

    with open(f'{Dir}/{name}.dat','w') as file:

        file.write('#Data file for LAMMPS \n')
        file.write('\n')
        file.write(f'{N} atoms \n')
        file.write(f'{len(types)} atom types \n')
        file.write('\n')
        file.write(f'{poss[:,0].min()-5}   {poss[:,0].max()+5}  xlo xhi \n')
        file.write(f'{poss[:,1].min()-5}   {poss[:,1].max()+5}  ylo yhi \n')
        file.write(f'{poss[:,2].min()-5}   {poss[:,2].max()+5}  zlo zhi \n \n')
        file.write('Masses \n \n')

        for i,tp in enumerate(types):
            file.write(f'{i+1} {masses[tp]} \n')

        file.write('\n')
        file.write('Atoms \n \n')

        index = 0

        for x, elem in zip(poss, Elem):
            index += 1
            file.write(f'    {index}   {types.searchsorted(elem)+1}   0   {np.round(x[0], 5)}    {np.round(x[1], 5)}    {np.round(x[2], 5)}\n')

    return

def assemble(objs):
    '''Concatenate bottom electrode, molecule and top electrode, each sorted by z.
    z_order maps output position -> index in the concatenated input.'''
    counts = tuple(obj.N for obj in objs)
    N = sum(counts)
    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    Extremes = np.zeros(N, dtype = 'bool')
    z_order = np.zeros(N, dtype = 'int64')
    n = 0
    for obj in objs:
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        Extremes[n:obj.N+n] = obj.extremes.copy()
        z_order[n:obj.N+n] = obj.poss[:,2].argsort()+n
        n += obj.N
    return Poss[z_order], Elem[z_order], Semi_infinite_layer[z_order], Extremes[z_order], z_order, counts

def dftu_block(objs, DFT_opts):
    if not DFT_opts['DFTU']:
        return None
    return [objs[0].N+1, objs[0].N+objs[1].N, objs[1].U[DFT_opts['Functional']]]

def verify_output(Dir):
    '''Raise if any generated helper script is not valid Python. Returns the number of step folders.'''
    problems = []
    for folder, _, files in walk(Dir):
        for f in files:
            if f.endswith('.py'):
                try:
                    with open(f'{folder}/{f}') as src:
                        ast.parse(src.read())
                except SyntaxError as e:
                    problems.append(f'{folder}/{f}, line {e.lineno}: {e.msg}')
    if problems:
        raise RuntimeError('Generated scripts are not valid Python:\n' + '\n'.join(problems))
    return len([x for x in listdir(Dir) if path.isdir(f'{Dir}/{x}')])

def Create_output(objs, Dir, name, opt_data, DFT_opts):
    begin_export()
    Poss, Elem, Semi_infinite_layer, Extremes, z_order, counts = assemble(objs)
    opt_resume = optimization_protocol(opt_data, Semi_infinite_layer, z_order, counts)

    text = f'Nbtm, Nmol, Ntop = {counts[0]}, {counts[1]}, {counts[2]} \n'

    create_gjf(name, Dir, Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)
    create_python_conector(Dir, name, DFT_opts['EABS'] ,text)
    create_unisa(name, Dir)
    create_ant(name, Dir, Elem[0], Elem[-1], dftu_block(objs, DFT_opts))
    create_xyz(Dir, name, Poss, Elem, len(Elem))
    return


def Create_reax_output(objs, Dir, name):
    Poss, Elem, _, _, _, _ = assemble(objs)
    types = np.unique(Elem)
    create_reax(Dir, name, Poss, Elem, len(Elem), types)
    create_xyz(Dir, name, Poss, Elem, len(Elem))
    return

def Create_xyz_output(objs, Dir, name):
    Poss, Elem, _, _, _, _ = assemble(objs)
    create_xyz(Dir, name, Poss, Elem, len(Elem))
    return

def molecule_file_atoms(molecule):
    '''Atom lines of the molecule exactly as in its file (not centred, not rotated).'''
    with open(f'Molecules/{molecule.name}.xyz') as text:
        lines = text.read().split('\n')
    return '\n'.join(lines[2:2+molecule.N]) + '\n'

def Create_getU_output(objs, Dir, name, DFT_opts):
    '''Gaussian single points of the isolated molecule with N, N-1 and N+1 electrons, and Get_U.py,
    which turns them into U and writes the molecule file with the U of this functional added.'''
    molecule = objs[1]
    Elem = molecule.elements
    atom_lines = molecule_file_atoms(molecule)
    electrons = count_electrons(Elem, np.zeros(len(Elem), dtype = 'bool'))
    logs = {}
    for tag, charge, mult, shell in gas_states(electrons, DFT_opts):
        create_gjf_gas(f'{molecule.name}_{tag}', Dir, Elem, atom_lines, charge, mult, shell, DFT_opts)
        logs[tag] = f'{molecule.name}_{tag}.log'
    create_launch_get_u(Dir, [x[:-4] for x in logs.values()])
    create_python_get_u(Dir, molecule.name, DFT_opts['Functional'], logs, dict(molecule.U), atom_lines)
    return

def Create_output_soc(objs, Dir, name, opt_data, DFT_opts): #Crear SOC dft
    begin_export()
    Poss, Elem, Semi_infinite_layer, Extremes, z_order, counts = assemble(objs)
    opt_resume = optimization_protocol(opt_data, Semi_infinite_layer, z_order, counts)

    text = f'Nbtm, Nmol, Ntop = {counts[0]}, {counts[1]}, {counts[2]} \n'

    create_gjf_soc(name, Dir, Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)
    create_python_conector(Dir, name, DFT_opts['EABS'], text)
    create_unisa(name, Dir)
    create_ant_soc(name, Dir, Elem[0], Elem[-1], dftu_block(objs, DFT_opts), Elem, Semi_infinite_layer) #Can I apply SOC to the second layer of the semi-infinite electrode? I think yes, but I have to check it
    create_xyz(Dir, name, Poss, Elem, len(Elem))
    return

def sequence_opt_resume(opt_data, Semi_infinite_layer, z_order, counts):
    if opt_data['Chains'] == 'none':
        return np.full(sum(counts), 'fix', dtype = '<U20'), False
    return optimization_protocol(opt_data, Semi_infinite_layer, z_order, counts), True

def Create_grid_output(objs, Dir, name, grid, opt_data, DFT_opts):
    begin_export()
    Poss, Elem, Semi_infinite_layer, Extremes, z_order, counts = assemble(objs)
    N = len(Elem)
    N_up = counts[2]

    x = np.linspace(grid[0][0], grid[0][1], grid[0][2])
    y = np.linspace(grid[1][0], grid[1][1], grid[1][2])
    n_points = grid[0][2]*grid[1][2]

    up_poss = objs[-1].poss.copy()[z_order[N-N_up:]-np.min(z_order[N-N_up:])]

    DFTU = dftu_block(objs, DFT_opts)
    opt_resume, opt_flag = sequence_opt_resume(opt_data, Semi_infinite_layer, z_order, counts)
    name_number = [str(i+1).zfill(len(str(n_points))) for i in range(n_points)]

    # Snake path over the grid: y runs forwards, then backwards on the next x
    Grid = []
    for i in x:
        for j in y:
            Grid.append((float(i), float(j)))
        y = y[::-1]

    POSS = []
    ELEM = []
    with open(f'{Dir}/grid_data.dat', 'w') as grid_dat:
        grid_dat.write(f'z = {float(grid[2])} \n \n')
        for count, (i, j) in enumerate(Grid):
            step = f'{name}_{name_number[count]}'
            step_dir = f'{Dir}/{step}'
            Poss[N-N_up:] = up_poss + np.array([i,j,grid[2]])
            mkdir(step_dir)

            if not opt_flag or count == 0:
                create_gjf(step, step_dir, Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)

            text = f'Nbtm, Nmol, Ntop = {counts[0]}, {counts[1]}, {N_up} \n'
            text += f'Grid = {Grid} \n'
            if opt_data['Chains'] == 'first step' and count == 0:
                text += f'I = 0 \n \n'
                text += f"Next_folder = '{name}_{'0'*len(str(n_points))}' \n"
                text += read_template('Move_grid.txt') + '\n'
                create_python_OptFirstStep(step_dir, step, name, n_points, text)
            elif opt_data['Chains'] == 'every step':
                text += f'I = {count+1} \n \n'  # move from this grid point to the next one
                text += read_template('Move_grid.txt') + '\n'
                nxt = f'{name}_{name_number[count+1]}' if count+1 < n_points else ''
                create_python_OptAllSteps(step_dir, step, nxt, text)

            create_ant(step, step_dir, Elem[0], Elem[-1], DFTU, DOSAT= N-N_up+1)
            POSS.append(Poss.copy())
            ELEM.append(Elem)
            grid_dat.write(f'{name_number[count]}: x = {i} , y = {j} \n')

    create_unisa_multi(n_points, Dir, opt_data['Chains'])
    create_xyz_film(Dir, name, POSS, ELEM, N)
    return

def Create_pull_output(objs, Dir, name, dz, nz, both, opt_data, DFT_opts):
    '''Step k moves the top electrode by k*dz, or each electrode by k*dz/2 in opposite directions if both.
    Either way the gap grows by dz per step.'''
    begin_export()
    Poss, Elem, Semi_infinite_layer, Extremes, z_order, counts = assemble(objs)
    N = len(Elem)
    N_down, N_mol, N_up = counts

    steps = np.arange(nz)
    if both:
        z_up = steps*dz/2
        z_down = -steps*dz/2
    else:
        z_up = steps*dz
        z_down = np.zeros(nz)

    up_poss = Poss[N-N_up:].copy()
    down_poss = Poss[:N_down].copy()

    DFTU = dftu_block(objs, DFT_opts)
    opt_resume, opt_flag = sequence_opt_resume(opt_data, Semi_infinite_layer, z_order, counts)
    name_number = [str(i+1).zfill(len(str(nz))) for i in range(nz)]

    POSS = []
    ELEM = []
    with open(f'{Dir}/pull_data.dat', 'w') as pull_dat:
        pull_dat.write(f'{dz} {nz} {both} \n \n')
        for count, (zu, zd) in enumerate(zip(z_up, z_down)):
            step = f'{name}_{name_number[count]}'
            step_dir = f'{Dir}/{step}'
            Poss[N-N_up:] = up_poss + np.array([0,0,zu])
            Poss[:N_down] = down_poss + np.array([0,0,zd])
            mkdir(step_dir)

            if not opt_flag or count == 0:
                create_gjf(step, step_dir, Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)

            text = f'Nbtm, Nmol, Ntop = {N_down}, {N_mol}, {N_up} \n'
            text += f'dz = {float(dz)} \n'
            text += f'both = {both} \n \n'
            if opt_data['Chains'] == 'first step' and count == 0:
                text += f"Next_folder = '{name}_{'0'*len(str(nz))}' \n"
                text += read_template('Move_pull.txt') + '\n'
                create_python_OptFirstStep(step_dir, step, name, nz, text)
            elif opt_data['Chains'] == 'every step':
                text += read_template('Move_pull.txt') + '\n'
                nxt = f'{name}_{name_number[count+1]}' if count+1 < nz else ''
                create_python_OptAllSteps(step_dir, step, nxt, text)

            create_ant(step, step_dir, Elem[0], Elem[-1], DFTU)
            POSS.append(Poss.copy())
            ELEM.append(Elem)
            pull_dat.write(f'{name_number[count]}: top dz = {zu:.4f} , bottom dz = {zd:.4f} , gap change = {zu-zd:.4f} \n')

    create_unisa_multi(nz, Dir, opt_data['Chains'])
    create_xyz_film(Dir, name, POSS, ELEM, N)
    return

def Create_rot_output(objs, Dir, name, angles, boolean_data, SOC, opt_data, DFT_opts):
    '''Rotate the selected components about the z axis through their centres, by each angle in
    angles, measured from the orientation they have in the viewer.'''
    begin_export()
    Poss, Elem, Semi_infinite_layer, Extremes, z_order, counts = assemble(objs)
    N = len(Elem)
    ns = (0, counts[0], counts[0]+counts[1], N)
    angles = [float(a) for a in angles]   # plain floats: numpy>=2 scalar reprs would leak into link_opt.py
    dtheta = angles[1]-angles[0] if len(angles) > 1 else 0.0
    width = len(str(len(angles)))

    DFTU = dftu_block(objs, DFT_opts)
    opt_resume, opt_flag = sequence_opt_resume(opt_data, Semi_infinite_layer, z_order, counts)
    write_gjf = create_gjf_soc if SOC else create_gjf

    start_angles = [obj.euler_angles.astype(float).copy() for obj in objs]
    POSS = []
    ELEM = []
    try:
        with open(f'{Dir}/rotation_data.dat', 'w') as rot_dat:
            rot_dat.write(f'From {angles[0]} to {angles[-1]} by steps of {dtheta} (degrees, relative to the starting orientation) \n \n')
            for i,psi in enumerate(angles):
                for j,flag in enumerate(boolean_data):
                    if flag:
                        objs[j].euler_angles = start_angles[j] + np.array([0., 0., psi])
                        Poss[ns[j]:ns[j+1]] = objs[j].true_coordinates()[z_order[ns[j]:ns[j+1]]-ns[j]]

                step = f'{name}_{str(i+1).zfill(width)}'
                step_dir = f'{Dir}/{step}'
                mkdir(step_dir)

                if not opt_flag or i == 0:
                    write_gjf(step, step_dir, Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)

                text = f'Nbtm, Nmol, Ntop = {counts[0]}, {counts[1]}, {counts[2]} \n'
                text += f'dtheta = {dtheta} \n'
                text += f'Flag_bottom = {boolean_data[0]} \n'
                text += f'Flag_mol = {boolean_data[1]} \n'
                text += f'Flag_top = {boolean_data[2]} \n'
                if opt_data['Chains'] == 'first step' and i == 0:
                    text += f"Next_folder = '{name}_{'0'*width}' \n"
                    text += read_template('Move_rot.txt') + '\n'
                    create_python_OptFirstStep(step_dir, step, name, len(angles), text)
                elif opt_data['Chains'] == 'every step':
                    text += read_template('Move_rot.txt') + '\n'
                    nxt = f'{name}_{str(i+2).zfill(width)}' if i+1 < len(angles) else ''
                    create_python_OptAllSteps(step_dir, step, nxt, text)

                if SOC:
                    create_ant_soc(step, step_dir, Elem[0], Elem[-1], DFTU, Elem, Semi_infinite_layer)
                else:
                    create_ant(step, step_dir, Elem[0], Elem[-1], DFTU)

                POSS.append(Poss.copy())
                ELEM.append(Elem)
                rot_dat.write(f'{str(i+1).zfill(width)}: {psi} \n')
    finally:
        for obj, start in zip(objs, start_angles):
            obj.euler_angles = start
            obj.true_coordinates()

    create_unisa_multi(len(angles), Dir, opt_data['Chains'])
    create_xyz_film(Dir, name, POSS, ELEM, N)
    return
