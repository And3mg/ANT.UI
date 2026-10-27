molecule = 'Glycerol' 
functional = 'PBE' 
logs = {'N': 'Glycerol_N.log', 'cation': 'Glycerol_cation.log', 'anion': 'Glycerol_anion.log'} 
U_pairs = {'PBE': '11.86'} 
xyz_atoms = 'C    3.3246    1.7612    0.9530\nC    2.2775    1.3270    1.9749\nC    1.2844    2.4108    2.3852 \nO    1.5546    0.1848    1.4626 \nO    4.1811    0.6876    0.6145 \nO    0.3392    1.9107    3.3112 \nH    2.8136    2.1659    0.0582 \nH    3.9507    2.5626    1.3605 \nH    2.7812    0.9421    2.8673 \nH    1.8044    3.2457    2.8677 \nH    0.7880    2.8114    1.4804 \nH    1.0866    0.4775    0.6601 \nH    3.6052   -0.0933    0.5441 \nH    0.1107    1.0193    2.9953 \n' 

# Hubbard U of the isolated molecule from three single points (N, N-1 and N+1 electrons):
#   IP = E(N-1) - E(N),  EA = E(N) - E(N+1)
#   Sigma_occ   = IP + e_HOMO        (how far -e_HOMO is from IP, Koopmans' theorem)
#   Sigma_unocc = -(EA + e_LUMO)     (how far -e_LUMO is from EA)
#   U = Sigma_occ + Sigma_unocc = (IP - EA) - (e_LUMO - e_HOMO)
# ANT shifts the molecular levels by -U(n - 1/2): occupied down U/2, empty up U/2, so the
# gap grows by U and the Kohn-Sham gap becomes the gas-phase fundamental gap IP - EA.
import re

HARTREE = 27.211386

def read_log(file):
    '''Last SCF energy (Hartree) and frontier orbitals (HOMO, LUMO in Hartree; alpha and beta together).'''
    E, occ, virt, normal = None, [], [], False
    with open(file, 'r') as data:
        for line in data:
            if 'SCF Done' in line:
                E = float(line.split('=')[1].split()[0])
            elif 'Population analysis' in line:
                occ, virt = [], []
            elif 'occ. eigenvalues' in line:
                occ += [float(x) for x in re.findall(r'-?\d+\.\d+', line.split('--')[1])]
            elif 'virt. eigenvalues' in line:
                virt += [float(x) for x in re.findall(r'-?\d+\.\d+', line.split('--')[1])]
            elif 'Normal termination' in line:
                normal = True
    if not normal or E is None:
        raise SystemExit(f'{file}: Gaussian did not finish normally, U not computed.')
    return E, (max(occ) if occ else None), (min(virt) if virt else None)

if __name__ == '__main__':
    E_N, homo, lumo = read_log(logs['N'])
    E_cat = read_log(logs['cation'])[0]
    E_an = read_log(logs['anion'])[0]
    IP = (E_cat - E_N)*HARTREE
    EA = (E_N - E_an)*HARTREE
    homo, lumo = homo*HARTREE, lumo*HARTREE
    S_occ = IP + homo
    S_unocc = -(EA + lumo)
    U = S_occ + S_unocc

    notes = []
    if S_occ < 0 or S_unocc < 0:
        notes.append('One correction is negative: the symmetric shift applied by ANT does not describe this molecule well.')
    if abs(S_occ - S_unocc) > 0.25*abs(U):
        notes.append('The two corrections differ by more than 25 %: ANT applies the same U/2 to occupied and empty levels.')

    with open(f'U.{molecule}.dat', 'w') as out:
        out.write(f'Hubbard U of {molecule} ({functional}, isolated molecule) \n')
        out.write(f'E(N) = {E_N:.8f}  E(N-1) = {E_cat:.8f}  E(N+1) = {E_an:.8f} Hartree \n')
        out.write(f'IP = {IP:.3f} eV   EA = {EA:.3f} eV   IP - EA = {IP-EA:.3f} eV \n')
        out.write(f'e_HOMO = {homo:.3f} eV   e_LUMO = {lumo:.3f} eV   Kohn-Sham gap = {lumo-homo:.3f} eV \n')
        out.write(f'Sigma_occ = IP + e_HOMO = {S_occ:.3f} eV   Sigma_unocc = -(EA + e_LUMO) = {S_unocc:.3f} eV \n')
        out.write(f'U = Sigma_occ + Sigma_unocc = {U:.3f} eV \n')
        for note in notes:
            out.write(f'Note: {note} \n')
        out.write(f'{molecule}.xyz holds this U; copy it to the Molecules folder of ANT.UI. \n')
        out.write('In the junction ANT uses UPLUS = U x U_ELECTRON_SCREENING (USER_CONF.txt): choose the screening \n')
        out.write('for the molecule and contact geometry (image-charge renormalisation by the electrodes). \n')

    U_pairs[functional] = f'{U:.2f}'
    atoms = [line for line in xyz_atoms.split('\n') if line.strip()]
    with open(f'{molecule}.xyz', 'w') as out:
        out.write(f'{len(atoms)}\n')
        out.write('|'.join(f'{f} {u}' for f, u in U_pairs.items()) + '\n')
        out.write('\n'.join(atoms) + '\n')
    print(f'U({molecule}, {functional}) = {U:.2f} eV')
