from ensemble_analyser import Conformer, get_calculator
from openbabel import openbabel as ob
import numpy as np
import os, cclib, subprocess
from ase.optimize import LBFGS
import re
from rmsd import quaternion_rmsd as rmsd


protocol = {
    0: {
        'func': 'pbe',
        'basis': '3-21g', 
        'opt': False,
        'freq': False,
        'solv' : {
            'solv': 'toluene',
            'smd' : False,
        },
    },
    1: {
        'func': 'pbe',
        'basis': 'def2-svp', 
        'opt': True,
        'freq': False,
        'solv' : {
            'solv': 'toluene',
            'smd' : False,
        },
    }
}


def parse_rotational_const(x):
    if re.search('Rotational constants in cm-1', x):
        return x

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


def get_conf_parameters(conf, number):


    data = cclib.io.ccread('ORCA.out')
    e, g = data.__dict__.get('scfenergies', None)[-1], data.__dict__.get('freeenergy', None)
    with open('ORCA.out') as f:
        fl = f.readlines()

    B = np.linalg.norm(np.array(list(filter(parse_rotational_const, fl))[-1].strip().split(':')[-1].split(), dtype=float))

    conf.energies[number] = {
        'E' : e*23.060541945329 if e else e,
        'G' : g*627.50960803059 if g else g,
        'B' : B if B else None,
    }
    return None


def launch(number, conf, protocol):

    if protocol['opt']:
        calculator = get_calculator(protocol, opt=True)
        atm = conf.get_ase_atoms(calculator)
        opt = LBFGS(atm, logfile='ORCA_ase.log', trajectory='ORCA_ase.trj')
        opt.run()
        conf.last_geometry = opt.atoms.get_positions()

    if protocol['freq']:
        calculator = get_calculator(protocol)
        atm = conf.get_ase_atoms(calculator)
        atm.get_potential_energy()
    
    if not (protocol['freq'] and protocol['opt']):
        calculator = get_calculator(protocol)
        atm = conf.get_ase_atoms(calculator)
        atm.get_potential_energy()

    get_conf_parameters(conf, number)
    

def check_ensemble(confs):
    for idx, i in enumerate(confs):
        for j in range(0, idx):
            print(rmsd(i.last_geometry, confs[j].last_geometry))
        print(i.number, i.energies)


    return confs


def save_snapshot(output, confs):
    with open(output, 'w') as f:
        for i in confs:
            f.write(f'{i}\n')


def main(file):
    confs = read_ensemble(file)

    for idx, p in protocol.items():
        for i in confs:
            launch(idx, i, p)

        confs = sorted(confs)
        confs = check_ensemble(confs)
        save_snapshot(f'ensemble_after_{idx}.xyz', confs)

    save_snapshot('final_ensemble.xyz', confs)

    return None


if __name__ == '__main__':
    
    main('tests/struct/cluster_water.xyz') 


    # os.system('rm ORCA*')