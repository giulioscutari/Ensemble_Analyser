#!/usr/bin/python3

import cclib, re
import numpy as np

def parse_rotational_const(x):
    if re.search('Rotational constants in cm-1', x):
        return x
    
def parse_dipole_moment(x):
    if re.search('Total Dipole Moment', x):
        return x

def get_conf_parameters(conf, number, time):

    data = cclib.io.ccread('ORCA.out')
    e, g = data.__dict__.get('scfenergies', None)[-1], data.__dict__.get('freeenergy', None)

    with open('ORCA.out') as f:
        fl = f.readlines()

    B = np.linalg.norm(np.array(list(filter(parse_rotational_const, fl))[-1].strip().split(':')[-1].split(), dtype=float))

    M = np.linalg.norm(np.array(list(filter(parse_dipole_moment, fl))[-1].strip().split(':')[-1].split(), dtype=float))

    conf.energies[int(number)] = {
        'E' : e*23.060541945329 if e else e, # Electronic Energy [kcal/mol]
        'G' : g*627.50960803059 if g else g, # Free Gibbs Energy [kcal/mol]
        'B' : B if B else None,              # Rotatory Constant [cm-1]
        'm' : M if M else None,              # dipole momenti [Debye]
        'time' : time,                       # elapsed time [sec] 
    }

    return None