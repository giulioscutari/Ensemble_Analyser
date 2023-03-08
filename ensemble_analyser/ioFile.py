from .conformer import Conformer
import os
import numpy as np

def _parse_xyz_str(fl):
    fl = fl[2:]
    atoms, geom = [] , []
    for line in fl: 
        a, *g = line.split()
        atoms.append(a)
        geom.append(g)
    return np.array(atoms),np.array(geom, dtype=float)

def convert_file(file):
    output = '_'.join(file.split('.')[:-1])+'.xyz'
    os.system(f'obabel {file} -O{output}')
    return output

def read_ensemble(file):

    confs = []

    if not file.endswith('.xyz'):
        file = convert_file(file)

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