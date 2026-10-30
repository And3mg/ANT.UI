import numpy as np
from os import listdir

#%%
from _DICS import *
from _DFT import create_gjf

# SOCFACATOM factors P, D, F and RCUT (bohr). Au, Pt, Ag: fitted (2026) to two-component CRYSTAL23 bands with the
# SOC_params basis sets; Al, Cu (all-electron bases) keep 1, 1, 0. Method and per-element fits in
# SOC_params/antg-basis-sets-and-socfacs.pdf. The factors belong to these basis sets.
SOC_PARAMS = {'Au': ['201.3', '43.2', '12.6', '2.0'],
              'Pt': ['229.4', '41.6',  '9.3', '2.0'],
              'Ag': ['49.1',   '8.9',  '0.1', '2.0'],
              'Al': ['1.0',    '1.0',  '0.0', '2.5'],
              'Cu': ['1.0',    '1.0',  '0.0', '2.5']}
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
    '''Same input as create_gjf, with the SOC basis (and SOC ECP, when the element has one) on the active electrode atoms.'''
    create_gjf(name, Dir, Poss, Elem, Semi, Opt_resume, Opt_data, DFT_opts,
               active_basis = SOC_basis, active_ecp = SOC_pseudopot)
