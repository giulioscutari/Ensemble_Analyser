#!/usr/bin/python3

from ensemble_analyser.IOsystem import _parse_xyz_str
from ensemble_analyser.conformer import Conformer
import os

def convert_file(file):
    output = '_'.join(file.split('.')[:-1])+'.xyz'
    os.system(f'obabel {file} -O{output}')
    return output

def read_ensemble(file, charge, multiplicity):

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
        confs.append(Conformer(i//n_atoms, geom=geom, atoms=atoms, charge=charge, mult=multiplicity))
        old_idx = i

    return confs


