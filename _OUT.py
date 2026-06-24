from _DICS import *
from _DFT import *
from _SOC import *
import numpy as np
from os import listdir, mkdir


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

def Create_output(objs, Dir, name, opt_data, DFT_opts):
    N = 0
    for obj in objs:
        N += obj.N
        
    Molecular_atoms = objs[1].Molecular_atoms.copy()
    
    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    Extremes = np.zeros(N, dtype = 'bool')
    z_order = np.zeros(N, dtype = 'int64')
    n = 0
    for k,obj in enumerate(objs):
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        Extremes[n:obj.N+n] = obj.extremes.copy()
        if n == 0:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)
        else:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)+1
        n += obj.N
    #z_order = Poss[:,2].argsort()
    
    Poss = Poss[z_order]
    Elem = Elem[z_order]
    Semi_infinite_layer = Semi_infinite_layer[z_order]
    Extremes = Extremes[z_order]
    opt_resume = optimization_protocol(opt_data, Elem, Semi_infinite_layer, z_order, N)
    
    DFTU = None
    if DFT_opts['DFTU']:
        DFTU = [objs[0].N+1, objs[0].N+objs[1].N, objs[1].U[DFT_opts['Functional']]]
    
    text = f'Nbtm, Nmol, Ntop = {objs[0].N}, {objs[1].N}, {objs[-1].N} \n'

    create_gjf(name, Dir, Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)
    create_python_conector(Dir, name, DFT_opts['EABS'] ,text)
    create_unisa(name, Dir)
    create_ant(name, Dir, Elem[0], Elem[-1], DFTU)
    create_xyz(Dir, name, Poss, Elem, N)
    return

def optimization_protocol(opt_data, Elem, Semi_infinite_layer, z_order, N):
    opt_resume = np.zeros(N, dtype = '<U20')
    flag_top = False
    first_mol = True

    for i in range(N):
        if Semi_infinite_layer[i]:
            if flag_top: opt_resume[i] = opt_data['Top_Bethe']
            else: opt_resume[i] = opt_data['Bottom_Bethe']
        elif Elem[i] in Electrode_elem:
            if flag_top: opt_resume[i] = opt_data['Top']
            else: opt_resume[i] = opt_data['Bottom']
        elif Elem[i] in Molecular_atoms:
            flag_top = True
            if opt_data['Molecule'] in ['all', 'fix', 'z']: opt_resume[i] = opt_data['Molecule']
            else:
                if first_mol: 
                    n = i
                    Mol_z_atoms = [np.where( x+n == z_order)[0][0] for x in opt_data['Molecule_z_atoms']]
                    first_mol = False
                if i in Mol_z_atoms:
                    opt_resume[i] = 'z'
                else: opt_resume[i] = 'all'
    return opt_resume


def Create_reax_output(objs, Dir, name):
    N = 0
    for obj in objs:
        N += obj.N
    
    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    n = 0
    for obj in objs:
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        n += obj.N
    
    z_order = Poss[:,2].argsort()
    
    Poss = Poss[z_order]
    Elem = Elem[z_order]
    
    types = np.unique(Elem)
    
    create_reax(Dir, name, Poss, Elem, N, types)
    create_xyz(Dir, name, Poss, Elem, N)
    return

def Create_xyz_output(objs, Dir, name):
    N = 0
    for obj in objs:
        N += obj.N
    
    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    n = 0
    for obj in objs:
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        n += obj.N
    
    z_order = Poss[:,2].argsort()
    
    Poss = Poss[z_order]
    Elem = Elem[z_order]
    
    create_xyz(Dir, name, Poss, Elem, N)
    return

def Create_output_soc(objs, Dir, name, opt_data, DFT_opts): #Crear SOC dft
    N = 0
    for obj in objs:
        N += obj.N
    
    Molecular_atoms = objs[1].Molecular_atoms.copy()

    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    Extremes = np.zeros(N, dtype = 'bool')
    z_order = np.zeros(N, dtype = 'int64')

    n = 0
    for k,obj in enumerate(objs):
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        Extremes[n:obj.N+n] = obj.extremes.copy()
        if n == 0:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)
        else:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)+1
        n += obj.N
        
    #z_order = Poss[:,2].argsort()
        
    Poss = Poss[z_order]
    Elem = Elem[z_order]
    Semi_infinite_layer = Semi_infinite_layer[z_order]
    Extremes = Extremes[z_order]
    opt_resume = optimization_protocol(opt_data, Elem, Semi_infinite_layer, z_order, N)
    
    DFTU = None
    if DFT_opts['DFTU']:
        DFTU = [objs[0].N+1, objs[0].N+objs[1].N, objs[1].U[DFT_opts['Functional']]]
                
    create_gjf_soc(name, Dir, Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts) 
    create_unisa(name, Dir)
    create_ant_soc(name, Dir, Elem[0], Elem[-1], DFTU, Elem, Semi_infinite_layer) #Can I apply SOC to the second layer of the semi-infinite electrode? I think yes, but I have to check it
    create_xyz(Dir, name, Poss, Elem, N)
    create_python_opt_conector(Dir,name)
    return

def Create_grid_output(objs, Dir, name, grid, opt_data, DFT_opts):
    N = 0
    for obj in objs:
        N += obj.N
    
    Molecular_atoms = objs[1].Molecular_atoms.copy()

    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    Extremes = np.zeros(N, dtype = 'bool')
    z_order = np.zeros(N, dtype = 'int64')
    n = 0
    
    N_up = objs[-1].N

    x = np.linspace(grid[0][0], grid[0][1], grid[0][2])
    y = np.linspace(grid[1][0], grid[1][1], grid[1][2])
        
    for obj in objs:
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        Extremes[n:obj.N+n] = obj.extremes.copy()
        if n == 0:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)
        else:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)+1
        n += obj.N
    
    #z_order = Poss[:,2].argsort()
    
    Poss = Poss[z_order]
    Elem = Elem[z_order]
    Semi_infinite_layer = Semi_infinite_layer[z_order]
    Extremes = Extremes[z_order]
    
    up_poss = objs[-1].poss.copy()[z_order[N-N_up:]-np.min(z_order[N-N_up:])]
    
    count = 0
    Names = []
    POSS = []
    ELEM = []
    grid_dat = open(f'{Dir}/grid_data.dat', 'w')
    grid_dat.write(f'z = {grid[2]} \n \n')
    
    DFTU = None
    if DFT_opts['DFTU']:
        DFTU = [objs[0].N+1, objs[0].N+objs[1].N, objs[1].U[DFT_opts['Functional']]]
                
    if opt_data['Chains'] == 'none':
        opt_resume = np.zeros(N, dtype = '<U20')
        opt_resume[:] = 'fix'
        opt_flag = False
    else:
        opt_resume = optimization_protocol(opt_data, Elem, Semi_infinite_layer, z_order, N)
        opt_flag = True
        first = True

    
    name_number = []
    for i in range(grid[0][2]*grid[1][2]):
        i = ('0'*(len(str(grid[0][2]*grid[1][2]))-len(str(i+1)))) + str(i+1)
        name_number.append(i)

    Grid = []
    for i in x:
        for j in y:
            Grid.append((i,j))
        y = y[::-1]

    y = np.linspace(grid[1][0], grid[1][1], grid[1][2]) 

    for i in x:
        for j in y:
            Poss[N-N_up:] = up_poss + np.array([i,j,grid[2]])
            mkdir(Dir+f'/{name}_{name_number[count]}')
            Names.append(f'{name}_{name_number[count]}')
            
            if opt_flag:
                if first:
                    first = False
                    create_gjf(Names[-1], f'{Dir}/{name}_{name_number[count]}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)

                    if opt_data['Chains'] == 'first step':
                        text = f'Nbtm, Nmol, Ntop = {objs[0].N}, {objs[1].N}, {N_up} \n'
                        text += f'Grid = {Grid} \n'
                        text += f'I = 0 \n \n'
                        text += f"Next_folder = '{name}_{"0"*(len(str(grid[0][2]*grid[1][2])))}' \n"                

                        with open('PythonRutines/Move_grid.txt', 'r') as read:
                            text += read.read()   
                        text += '\n'
                        create_python_OptFirstStep(f'{Dir}/{name}_{name_number[count]}', Names[-1], name, grid[0][2]*grid[1][2], text)
                    
                if opt_data['Chains'] == 'every step':
                    text = f'Nbtm, Nmol, Ntop = {objs[0].N}, {objs[1].N}, {N_up} \n'
                    text += f'Grid = {Grid} \n'
                    text += f'I = {count} \n \n'
                    with open('PythonRutines/Move_grid.txt', 'r') as read:
                        text += read.read()   
                    text += '\n'
                    try:
                        nxt = f'{name}_{name_number[count+1]}'
                    except:
                        nxt = ''
                    create_python_OptAllSteps(f'{Dir}/{name}_{name_number[count]}', Names[-1], nxt, text)
            else:
                create_gjf(Names[-1], f'{Dir}/{name}_{name_number[count]}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)

            create_ant(Names[-1], f'{Dir}/{name}_{name_number[count]}', Elem[0], Elem[-1], DFTU, n-obj.N +1)
            POSS.append(Poss.copy())
            ELEM.append(Elem)
            grid_dat.write(f'{name_number[count]}: x = {i} , y = {j} \n')

            count+=1
        y = y[::-1]
    grid_dat.close()
    create_unisa_multi(count, Dir, opt_data['Chains'])
    create_xyz_film(Dir, name, POSS, ELEM, N)
    return

def Create_pull_output(objs, Dir, name, dz, nz, both, opt_data, DFT_opts):
    N = 0
    for obj in objs:
        N += obj.N
    
    Molecular_atoms = objs[1].Molecular_atoms.copy()

    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    Extremes = np.zeros(N, dtype = 'bool')
    z_order = np.zeros(N, dtype = 'int64')
    n = 0
    
    N_up = objs[-1].N
    N_down = objs[0].N

    if both:
        z_up = np.arange(0, dz*nz/2, dz/2)
        z_down = -np.arange(0, dz*nz/2, dz/2)
    else:
        z_up = np.arange(0, dz*nz, dz)
        z_down = np.zeros(nz)

            
    for obj in objs:
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        Extremes[n:obj.N+n] = obj.extremes.copy()
        if n == 0:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)
        else:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)+1
        n += obj.N
    
    #z_order = Poss[:,2].argsort()
    
    Poss = Poss[z_order]
    Elem = Elem[z_order]
    Semi_infinite_layer = Semi_infinite_layer[z_order]
    Extremes = Extremes[z_order]
    
    up_poss = objs[-1].true_coordinates().copy()[z_order[N-N_up:]-np.min(z_order[N-N_up:])]
    down_poss = objs[0].true_coordinates().copy()[z_order[:N_down]-np.min(z_order[:N_down])]
    
    count = 0
    Names = []
    POSS = []
    ELEM = []
    pull_dat = open(f'{Dir}/pull_data.dat', 'w')
    pull_dat.write(f'{dz} {nz} {both} \n \n')
    
    DFTU = None
    if DFT_opts['DFTU']:
        DFTU = [objs[0].N+1, objs[0].N+objs[1].N, objs[1].U[DFT_opts['Functional']]]
                
    if opt_data['Chains'] == 'none':
        opt_resume = np.zeros(N, dtype = '<U20')
        opt_resume[:] = 'fix'
        opt_flag = False
    else:
        opt_resume = optimization_protocol(opt_data, Elem, Semi_infinite_layer, z_order, N)
        opt_flag = True
        first = True

    name_number = []
    for i in range(nz):
        i = ('0'*(len(str(nz))-len(str(i+1)))) + str(i+1)
        name_number.append(i)

    for zu,zd in zip(z_up, z_down):
        
        Poss[N-N_up:] = up_poss + np.array([0,0,zu])
        Poss[:N_down] = down_poss + np.array([0,0,zd])
        mkdir(Dir+f'/{name}_{name_number[count]}')
        Names.append(f'{name}_{name_number[count]}')
        
        if opt_flag:
            if first:
                first = False
                create_gjf(Names[-1], f'{Dir}/{name}_{name_number[count]}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)

                if opt_data['Chains'] == 'first step':
                    text = f'Nbtm, Nmol, Ntop = {N_down}, {objs[1].N}, {N_up} \n'
                    text += f'dz = {dz} \n'
                    text += f'both = {both} \n \n'
                    text += f'Next_folder = {name}_{str(i).zfill(len(str(len(nz))))} \n'                 

                    with open('PythonRutines/Move_pull.txt', 'r') as read:
                        text += read.read()   
                    text += '\n'
                    create_python_OptFirstStep(f'{Dir}/{name}_{name_number[count]}', Names[-1], name, nz, text)

            if opt_data['Chains'] == 'every step':
                text = f'Nbtm, Nmol, Ntop = {N_down}, {objs[1].N}, {N_up} \n'
                text += f'dz = {dz} \n'
                text += f'both = {both} \n \n'
                with open('PythonRutines/Move_pull.txt', 'r') as read:
                    text += read.read()   
                text += '\n'
                try:
                    nxt = f'{name}_{name_number[count+1]}'
                except:
                    nxt = ''
                create_python_OptAllSteps(f'{Dir}/{name}_{name_number[count]}', Names[-1], nxt, text)

        else:
            create_gjf(Names[-1], f'{Dir}/{name}_{name_number[count]}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)
        create_ant(Names[-1], f'{Dir}/{name}_{name_number[count]}', Elem[0], Elem[-1], DFTU)
        POSS.append(Poss.copy())
        ELEM.append(Elem)
        count+=1
            
    pull_dat.close()
    create_unisa_multi(count, Dir, opt_data['Chains'])
    create_xyz_film(Dir, name, POSS, ELEM, N)
    return

def Create_rot_output(objs, Dir, name, angles, boolean_data, SOC, opt_data, DFT_opts):
    N = 0
    for obj in objs:
        N += obj.N
 
    Molecular_atoms = objs[1].Molecular_atoms.copy()        

    Poss = np.zeros((N,3))
    Elem = np.zeros(N, dtype = '<U2')
    Semi_infinite_layer = np.zeros(N, dtype = 'bool')
    Extremes = np.zeros(N, dtype = 'bool')
    z_order = np.zeros(N, dtype = 'int64')
    n = 0
    
    ns = (0, objs[0].N, objs[0].N+objs[1].N, N)


    for obj in objs:
        Poss[n:obj.N+n] = obj.true_coordinates()
        Elem[n:obj.N+n] = obj.elements.copy()
        Semi_infinite_layer[n:obj.N+n] = obj.semi_infinite_layer.copy()
        Extremes[n:obj.N+n] = obj.extremes.copy()
        if n == 0:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)
        else:
            z_order[n:obj.N+n] = obj.poss[:,2].argsort()+np.max(z_order)+1
        n += obj.N
        
    Poss = Poss[z_order]
    Elem = Elem[z_order]
    Semi_infinite_layer = Semi_infinite_layer[z_order]
    Extremes = Extremes[z_order]

    opt_resume = optimization_protocol(opt_data, Elem, Semi_infinite_layer, z_order, N)

    Names = []
    POSS = []
    ELEM = []
    rot_dat = open(f'{Dir}/rotation_data.dat', 'w')
    rot_dat.write(f'From {angles[0]} to {angles[-1]} by steps of {angles[1]-angles[0]} \n \n')
    
    DFTU = None
    if DFT_opts['DFTU']:
        DFTU = [objs[0].N+1, objs[0].N+objs[1].N, objs[1].U[DFT_opts['Functional']]]
    
    if opt_data['Chains'] == 'none':
        opt_resume = np.zeros(N, dtype = '<U20')
        opt_resume[:] = 'fix'
        opt_flag = False
    else:
        opt_resume = optimization_protocol(opt_data, Elem, Semi_infinite_layer, z_order, N)
        opt_flag = True
        first = True

    for i,psi in enumerate(angles):
        for j,flag in enumerate(boolean_data):
            if flag:
                objs[j].euler_angles[2] = psi  #Si metes rotar en x o y o z, es cambiar ese número a 0, 1 o 2
                Poss[ns[j]:ns[j+1]] = objs[j].true_coordinates()[z_order[ns[j]:ns[j+1]]-np.min(z_order[ns[j]:ns[j+1]])]

        mkdir(Dir+f'/{name}_{str(i+1).zfill(len(str(len(angles))))}')
        Names.append(f'{name}_{str(i+1).zfill(len(str(len(angles))))}')

        if opt_flag:
            if first:
                first = False
                if SOC:
                    create_gjf_soc(Names[-1], f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)
                else:
                    create_gjf(Names[-1], f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)

                if opt_data['Chains'] == 'first step':
                    text = f'Nbtm, Nmol, Ntop = {objs[0].N}, {objs[1].N}, {objs[-1].N} \n'
                    text += f'dtheta = {angles[1]-angles[0]} \n'
                    text += f'Flag_bottom = {boolean_data[0]} \n'
                    text += f'Flag_mol = {boolean_data[1]} \n'
                    text += f'Flag_top = {boolean_data[2]} \n'   
                    text += f'Next_folder = {name}_{str(i).zfill(len(str(len(angles))))} \n'                 
                    with open('PythonRutines/Move_rot.txt', 'r') as read:
                        text += read.read()   
                    text += '\n'
                    create_python_OptFirstStep(f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}', Names[-1], name, len(angles), text)
                    
            if opt_data['Chains'] == 'every step':
                text = f'Nbtm, Nmol, Ntop = {objs[0].N}, {objs[1].N}, {objs[-1].N} \n'
                text += f'dtheta = {angles[1]-angles[0]} \n'
                text += f'Flag_bottom = {boolean_data[0]} \n'
                text += f'Flag_mol = {boolean_data[1]} \n'
                text += f'Flag_top = {boolean_data[2]} \n'
                with open('PythonRutines/Move_rot.txt', 'r') as read:
                    text += read.read()   
                text += '\n'
                try:
                    nxt = f'{name}_{str(i+2).zfill(len(str(len(angles))))}'
                except:
                    nxt = ''
                create_python_OptAllSteps(f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}', Names[-1], nxt, text)

        else:
            if SOC:
                create_gjf_soc(Names[-1], f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)
            else:
                create_gjf(Names[-1], f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}', Poss, Elem, Extremes, opt_resume, opt_data, DFT_opts)
        if SOC:
            create_ant_soc(Names[-1], f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}',  Elem[0], Elem[-1], DFTU, Elem, Semi_infinite_layer)
        else:
            create_ant(Names[-1], f'{Dir}/{name}_{str(i+1).zfill(len(str(len(angles))))}', Elem[0], Elem[-1], DFTU)

        POSS.append(Poss.copy())
        ELEM.append(Elem)
            
    rot_dat.close()
    create_unisa_multi(len(angles), Dir, opt_data['Chains'])
    create_xyz_film(Dir, name, POSS, ELEM, N)
    return

