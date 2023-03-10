#!/usr/bin/python3

import cclib, re
import numpy as np

def parse_rotational_const(x):
    if re.search('Rotational constants in cm-1', x):
        return x
    
def parse_dipole_moment(x):
    if re.search('Total Dipole Moment', x):
        return x

def get_conf_parameters(conf, number):

    # log.debug(f'# Reading output file of CONF{conf.number}')
    data = cclib.io.ccread('ORCA.out')
    e, g = data.__dict__.get('scfenergies', None)[-1], data.__dict__.get('freeenergy', None)

    with open('ORCA.out') as f:
        fl = f.readlines()

    B = np.linalg.norm(np.array(list(filter(parse_rotational_const, fl))[-1].strip().split(':')[-1].split(), dtype=float))

    # print(list(filter(parse_dipole_moment, fl)))
    M = np.linalg.norm(np.array(list(filter(parse_dipole_moment, fl))[-1].strip().split(':')[-1].split(), dtype=float))
    # print(M)

    conf.energies[int(number)] = {
        'E' : e*23.060541945329 if e else e, # Electronic Energy in kcal/mol
        'G' : g*627.50960803059 if g else g, # Free Gibbs Energy in kcal/mol
        'B' : B if B else None,              # Rotatory Constant in cm-1
        'm' : M if M else None,              # dipole momenti in Debye
    }

    # log.debug(f'{conf.number} - {conf.energies}')

    return None