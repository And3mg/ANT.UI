from enum import Flag
import numpy as np
from os import listdir

#%%
from _DICS import *
from _DFT import opt_gjf, Flag_opt

SOC_PARAMS = {'Au': ['150.0', '35.0',  '5.0',  '2.0'],
              'Pt': ['280.0',  '40.0',  '10.0',  '2.5'],
              'Ag': ['80.0',  '10.0',  '5.0',  '2.5'],
              'Al': ['1.0',  '1.0',  '0.0',  '2.5'],
              'Cu': ['1.0',  '1.0',  '0.0',  '2.5']}
SOC_basis = {}
SOC_pseudopot = {}

for item in listdir('SOC_params'):
    if item.split('.')[-1] == 'gbs':
        _elem = item.split('_')[0]
        with open(f'SOC_params/{item}', 'r') as dat_base:
            _basis = ''
            line = dat_base.readline()
            while line != '':
                _basis += line
                line = dat_base.readline()  
        if item.split('.')[0].split('_')[-1] == 'basis':
            SOC_basis.update({_elem: _basis})
        elif item.split('.')[0].split('_')[-1] == 'pseudopot':
            SOC_pseudopot.update({_elem: _basis})



def create_ant_soc(name, Dir, bl1, bl2, DFTU, Elem, Semi, DOSAT = None):
    out_name = Dir + '/input' + name+'.ant'
    h = open(out_name, "w")
    h.write(UserInfo['ANTpredet'])
    if DOSAT != None:
        h.write(f"LDOS_BEG = {DOSAT}  \n")
        h.write(f"LDOS_END = {DOSAT}  \n")
    h.write(f"BLPAR1 = {BLPAR[bl1]} \n")
    h.write(f"BLPAR2 = {BLPAR[bl2]} \n")

    if DFTU:
        h.write('DFT+U \n')
        h.write('MOLBLOCKS \n')
        h.write('1 \n')
        h.write(f'{DFTU[0]},{DFTU[1]} \n')
        h.write(f'  UPLUS = {UserInfo["U_ELECTRON_SCREENING"]*float(DFTU[-1]):.2f} \n')

    h.write('SOC \n')
    h.write('SOCFACATOM \n')
    
    N = np.sum([x in Electrode_elem for x in Elem]) - np.sum(Semi)
    h.write(f'{N} \n')
    
    for i, val in enumerate(Semi):
        if not val and Elem[i] in Electrode_elem: 
            #h.write(f'{i+1}, {SOC_PARAMS[Elec_elem][0]}, {SOC_PARAMS[Elec_elem][1]}, {SOC_PARAMS[Elec_elem][2]}, {SOC_PARAMS[Elec_elem][3]}')
            h.write(f'{i+1}, {SOC_PARAMS[Elem[i]][0]}, {SOC_PARAMS[Elem[i]][1]}, {SOC_PARAMS[Elem[i]][2]}, {SOC_PARAMS[Elem[i]][3]} \n') #For lanl2dz

    
    h.close()
    return



def create_gjf_soc(name, Dir, Poss, Elem, Semi, Opt_resume, Opt_data, DFT_opts):  
    g = open(f'{Dir}/input{name}.gjf',"w")
    #CABECERA DEL PROGRAMA
    P_line, Poss_lines, Variable_lines = opt_gjf(Poss, Elem, Semi, Opt_resume, Opt_data, DFT_opts)
    from _DFT import Flag_opt
    g.write(UserInfo['GAUSSIANHEAD'])
    if not Flag_opt:
        g.write(f"%Subst L502 {UserInfo['ANTDIRcall']} \n")
        g.write(f"%Subst L101 {UserInfo['ANTDIRcall']} \n")
    g.write(P_line)

    g.write("\n")
    g.write(f'input{name}')
    g.write("\n")
    g.write("\n")
    
    #NAu = np.sum(Elem == Elec_elem)
    NAu = np.sum([x in set(Electrode_elem).difference(set(Molecular_atoms)) for x in Elem])
    g.write(f'{NAu%2},1 \n')
    
    Metals = [x in set(Electrode_elem).difference(set(Molecular_atoms)) for x in np.unique(Elem)]
    Metals = np.unique(Elem)[Metals]
    
    '''
    for elem, poss in zip(Elem, Poss):
        g.write(f'{elem}   {poss[0]}    {poss[1]}    {poss[2]}\n')
    g.write('\n')
    '''
    g.write(Poss_lines)
    g.write(Variable_lines)
    
    for mol_at in Molecular_atoms:
        if (Elem == mol_at).any():
            g.write(f'{mol_at} 0 \n')
            g.write('lanl2dz \n')
            g.write('**** \n')

    for met in Metals:    
        for i, val in enumerate(Semi):
            if not val and Elem[i] == met: g.write(f'{i+1} ')
        g.write('\n'+f'{SOC_basis[met]} \n')
        g.write('**** \n')
    
    for met in Metals:
        for i, val in enumerate(Semi):
            if val and Elem[i] == met: g.write(f'{i+1} ')
        g.write('\n')
        g.write(DFT_min_basis[met])
        g.write('**** \n')
    
    g.write('\n')
    
    
    for mol_at in Molecular_atoms:
        if (Elem == mol_at).any():
            g.write(f'{mol_at} 0 \n')
            g.write('lanl2dz \n')

    try:
        SOC_pseudopot[met]
        for met in Metals:    
            for i, val in enumerate(Semi):
                if not val and Elem[i] == met: g.write(f'{i+1} ')
            g.write('\n'+f'{SOC_pseudopot[met]} \n')
    except:
        print('No SOC pseudopotentials found.')
    
    for met in Metals:
        for i, val in enumerate(Semi):
            if val and Elem[i] == met: g.write(f'{i+1} ')
        g.write('\n')
        g.write(DFT_min_pseudopotentials[met])
    
    g.write(4*'\n ')
    return
