

from ensemble_analyser.IOsystem import _parse_xyz_str
from ensemble_analyser.conformer import Conformer

import os

def convert_file(file) -> str:
    """
    Convert the input file into xyz multigeometry XYZ file.
    OPENBABEL is required
    
    file | str : input filename

    return | str : input converted filename
    """
    output = '_'.join(file.split('.')[:-1])+'.xyz'
    os.system(f'obabel {file} -O{output}')
    return output


def read_ensemble(file, charge, multiplicity, log) -> list:
    """
    Read the initial ensemble and return the ensemble list
    Not only XYZ file is supported. OBABEL is required
    
    file | str : initial ensemble file 
    charge | int : charge of the molecule
    multiplicity | int : multiplicity of the molecule
    log : logger instance

    return | list : whole ensemble list as Conformer instances
    """

    confs = []

    if not file.endswith('.xyz'):
        file = convert_file(file)

    with open(file) as f:
        fl = f.readlines()

    n_atoms = int(fl[0])
    old_idx = 0
    for i in range(0, len(fl)+1, n_atoms+2):
        if i==old_idx: continue
        atoms, geom = _parse_xyz_str(fl[old_idx:i])
        confs.append(Conformer(i//n_atoms, geom=geom, atoms=atoms, charge=charge, mult=multiplicity))
        old_idx = i

    return confs


