import numpy as np
import shutil, json
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



class SerialiseEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj.__dict__