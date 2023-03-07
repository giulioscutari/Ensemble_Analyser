from ensemble_analyser import Conformer, get_calculator
from openbabel import openbabel as ob
import numpy as np
from ase.optimize import LBFGS


protocol = {
    0: {
        'func': 'pbe',
        'basis': 'def2-svp', 
        'opt': False,
        'freq': False,
    },
    1: {
        'func': 'pbe',
        'basis': 'def2-svp', 
        'opt': True,
        'freq': False,
    }
}


def _parse_xyz_str(fl):
    fl = fl[2:]
    atoms, geom = [] , []
    for line in fl: 
        a, *g = line.split()
        atoms.append(a)
        geom.append(g)
    return np.array(atoms),np.array(geom, dtype=float)

def convert_file(file):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats(file.split('.')[-1].lower(), "xyz")
    mol = ob.OBMol()
    conv.ReadFile(mol, file)

    # Convert Mol to RDKit Mol
    return conv.WriteString(mol)


def read_ensemble(file):

    confs = []

    if not file.endswith('.xyz'):
        fl = convert_file(file)
    else:
        with open(file) as f:
            fl = f.read()

    fl = fl.splitlines()
    n_atoms = int(fl[0])
    old_idx = 0
    for i in range(0, len(fl)+1, n_atoms+2):
        if i==old_idx: continue
        atoms, geom = _parse_xyz_str(fl[old_idx:i])
        confs.append(Conformer(i//n_atoms, geom=geom, atoms=atoms))
        old_idx = i

    return confs






if __name__ == '__main__':
    confs = read_ensemble('tests/struct/cluster_water.xyz')

    print(confs)

    for idx, p in protocol.items():
        calc = get_calculator(func=p['func'], basis=p['basis'], opt=p['opt'], freq=p['freq'])
        print(calc)
        for i in confs:
            atm = i.get_ase_atoms(calc)
            if p['opt']:
                opt = LBFGS(atm)
                opt.run()
                
            i.energies[idx] = atm.get_potential_energy()

    for i in confs:
        print(i.energies)