import numpy as np

#%%
from _DICS import *

p_line = '#p HSEH1PBE/gen pseudo=read pop=min Scf=(conver=5,maxcycle=1000,dsymm,nodamp) nosymm IOP(3/124=3)'
Flag_opt = False
Flag_EABS = False
Flag_opt_in_ANT = False
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
    barra = '/'
    h = open(Dir+'/launchANT.sh','w', newline='\n')
    h.write(UserInfo['JOBHEAD'])
    h.write(f'name={Dir.split(barra)[1]}\n')
    h.write(f'printf -v first "%0{len(str(N))}d" 1 \n')
    h.write('# Launch program \n') 
    h.write(f'for i in $(seq -f "%0{len(str(N))}g" 1 {N})\n')
    h.write('do \n')
    h.write('  if [[ $i != $first ]] \n')
    h.write('  then \n')
    h.write('    cp P.* ../$name\_$i/P.input$name\_$i.dat \n')
    h.write('    cd ../$name\_$i \n')
    h.write(f"    {UserInfo['G09call']} input$name\_$i.gjf \n")
    if Opt_chain == 'every step':
        h.write(f"    {UserInfo['PYTHONcall']} link_opt.py \n")
        if Flag_opt_in_ANT:
            h.write(f"    find input$name\_$i.ant | xargs sed -i 's/EW1 = 3.0/EW1 = -3.0/g' \n")
        h.write(f"    {UserInfo['G09call']} input$name\_$i\_opt.gjf \n")      
    h.write('  else \n')
    h.write('    cd $name\_$i \n')
    h.write(f"    {UserInfo['G09call']} input$name\_$i.gjf \n")
    if Opt_chain != 'none':
        h.write(f"    {UserInfo['PYTHONcall']} link_opt.py \n")
        if Flag_opt_in_ANT:
            h.write(f"    find input$name\_$i.ant | xargs sed -i 's/EW1 = 3.0/EW1 = -3.0/g' \n")
        h.write(f"    {UserInfo['G09call']} input$name\_$i\_opt.gjf \n")

    h.write('  fi\n')
    h.write('done \n')
    h.write('cd .. \n')
    h.write('mkdir T; \n')
    h.write('mkdir xyz; \n')
    h.write('mkdir Tlog; \n')
    h.write('find . -type f -name "*opt.log" -exec cp {} Tlog/ \; \n')
    h.write('find . -type f -name "T.*" -exec cp {} T/ \; \n')
    h.write('find . -type f -name "input*.xyz" -exec cp {} xyz/ \; \n')
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

def opt_gjf(Poss, Elem, Semi, Opt_resume, Opt_data, DFT_opts):
    global p_line, Flag_opt
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
    flag = True
    iB,iT = 0,0
    Bs,Ts = [],[]
    for elem, poss, semi, opt in zip(Elem, Poss, Semi, Opt_resume):
        if opt == 'fix':
            Poss_lines += f'{elem} {character}  {poss[0]:.6f}    {poss[1]:.6f}    {poss[2]:.6f}\n'
            flag = False
        elif opt == 'all':
            Poss_lines += f'{elem} {character}  x{k_opt}    y{k_opt}    z{k_opt}\n'
            Variable_lines += f'x{k_opt} = {poss[0]:.6f} \n'
            Variable_lines += f'y{k_opt} = {poss[1]:.6f} \n'
            Variable_lines += f'z{k_opt} = {poss[2]:.6f} \n'
            k_opt +=1
            flag = False
        elif opt == 'z':
            if semi and flag:
                if not f"{poss[2]:.6f}" in Bs:
                    iB += 1
                    Variable_lines+= f'zB{iB} = {poss[2]:.6f} \n'
                    Bs.append(f"{poss[2]:.6f}")
                Poss_lines += f'{elem} {character}  {poss[0]:.6f}    {poss[1]:.6f}    zB{iB}\n'
            elif semi:
                if not f"{poss[2]:.6f}" in Ts:
                    iT += 1
                    Variable_lines+= f'zT{iT} = {poss[2]:.6f} \n'
                    Ts.append(f"{poss[2]:.6f}")
                Poss_lines += f'{elem} {character}  {poss[0]:.6f}    {poss[1]:.6f}    zT{iT}\n'
            else:
                Poss_lines += f'{elem} {character}  {poss[0]:.6f}    {poss[1]:.6f}    z{k_opt}\n'
                Variable_lines += f'z{k_opt} = {poss[2]:.6f} \n'
                k_opt += 1
    Variable_lines += '\n'
    return P_line, Poss_lines, Variable_lines

def create_gjf(name, Dir, Poss, Elem, Semi, Opt_resume, Opt_data, DFT_opts):  
    global Flag_opt, Flag_opt_in_ANT
    g = open(f'{Dir}/input{name}.gjf',"w")
    P_line, Poss_lines, Variable_lines = opt_gjf(Poss, Elem, Semi, Opt_resume, Opt_data, DFT_opts)
    g.write(UserInfo['GAUSSIANHEAD'])
    Flag_opt_in_ANT = Opt_data['Opt_in_ANT']
    if not Flag_opt or Flag_opt_in_ANT:
        g.write(f"%Subst L502 {UserInfo['ANTDIRcall']} \n")
        g.write(f"%Subst L101 {UserInfo['ANTDIRcall']} \n")
    g.write(P_line)
    g.write("\n")
    g.write(f'input{name}')
    g.write("\n")
    g.write("\n")
    
    NAu = np.sum([x in set(Electrode_elem).difference(set(Molecular_atoms)) for x in Elem])
    g.write(f'{NAu%2},1 \n')
    
    Metals = [x in set(Electrode_elem).difference(set(Molecular_atoms)) for x in np.unique(Elem)]
    Metals = np.unique(Elem)[Metals]

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
        g.write('\n'+f'{DFT_basis[met]}')
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

        
    for met in Metals:    
        for i, val in enumerate(Semi):
            if not val and Elem[i] == met: g.write(f'{i+1} ')
        g.write('\n'+f'{DFT_pseudopotentials[met]}')
    
    for met in Metals:
        for i, val in enumerate(Semi):
            if val and Elem[i] == met: g.write(f'{i+1} ')
        g.write('\n')
        g.write(DFT_min_pseudopotentials[met])
    
    g.write(4*'\n ')
    return 

def create_python_conector(Dir, name, FlagEABS, text):
    global Flag_opt, Flag_EABS, Flag_opt_in_ANT
    Flag_EABS = FlagEABS
    if not Flag_EABS and not Flag_opt:
        return
    global p_line
    with open('PythonRutines/CommonLink.txt', 'r') as read:
        common = read.read() 
    main = ''
    if Flag_opt:
        with open('PythonRutines/Main.txt', 'r') as read:
            main += read.read()
    if Flag_EABS:
        with open('PythonRutines/Main_EABS.txt', 'r') as read:
            main += read.read()
   
    with open(f'{Dir}/link_opt.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        out.write(f"optgjf_path = 'input{name}.gjf' \n")
        out.write(f"log_file_path = 'input{name}.log' \n")
        out.write(f"P_line = f'{p_line}")
        out.write(" {intro}' \n")
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
    main = ''
    with open('PythonRutines/Get_EABS.txt', 'r') as read:
        main += read.read()
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

 # Python script in script
def create_python_opt_conector(Dir, name):
    global p_line
    with open('PythonRutines/CommonLink.txt', 'r') as read:
        common = read.read()
    with open('PythonRutines/Main.txt', 'r') as read:
        main = read.read()

    with open(f'{Dir}/link_opt.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        out.write(f"optgjf_path = 'input{name}.gjf' \n")
        out.write(f"log_file_path = 'input{name}.log' \n")
        out.write(f"P_line = f'{p_line}")
        out.write(" {intro}' \n")
        out.write(f"Ant1 = %Subst L502 {UserInfo['ANTDIRcall']} \n")
        out.write(" {intro}' \n")
        out.write(f"Ant2 = %Subst L101 {UserInfo['ANTDIRcall']} \n")
        out.write("{intro}' \n")
        out.write(common)
        out.write(main)

def create_python_OptFirstStep(Dir, name, Folder_name, N_folder, text_func):
    global Flag_opt, Flag_opt_in_ANT
    global p_line
    with open('PythonRutines/CommonLink.txt', 'r') as read:
        common = read.read()
    with open('PythonRutines/Main_OptFirstStep.txt', 'r') as read:
        main = read.read()

    with open(f'{Dir}/link_opt.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        out.write(f"optgjf_path = 'input{name}.gjf' \n")
        out.write(f"log_file_path = 'input{name}.log' \n")
        out.write(f"P_line = f'{p_line}")
        out.write(" {intro}' \n")
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
    with open('PythonRutines/CommonLink.txt', 'r') as read:
        common = read.read()
    with open('PythonRutines/Main_OptAllSteps.txt', 'r') as read:
        main = read.read()

    with open(f'{Dir}/link_opt.py','w') as out:
        out.write(r"intro = '\n'")
        out.write('\n')
        out.write(f"optgjf_path = 'input{name}.gjf' \n")
        out.write(f"log_file_path = 'input{name}.log' \n")
        out.write(f"P_line = f'{p_line}")
        out.write(" {intro}' \n")
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
