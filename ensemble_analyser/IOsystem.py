import numpy as np
import shutil
import os


def _parse_xyz_str(fl):
    fl = fl[2:]
    atoms, geom = [] , []
    for line in fl:
        a, *g = line.split()
        atoms.append(a)
        geom.append(g)
    return np.array(atoms),np.array(geom, dtype=float)



def mkdir(directory:str):
    if os.path.exists(directory):
        raise IOError(f'Directory {directory} already exists. Going to exit!')
        shutil.rmtree(directory)
    os.mkdir(directory)
    return None

