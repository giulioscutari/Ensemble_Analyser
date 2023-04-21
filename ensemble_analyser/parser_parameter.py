#!/data/bin/python_env/bin/python3

# import cclib
import re, os, sys
import numpy as np
from ensemble_analyser.rrho import free_gibbs_energy

EH_TO_KCAL = 627.5096080305927



def parse_rotational_const(x):
    """
    ORCA parsing for Rotational Constant
    """
    if re.search('Rotational constants in cm-1', x):
        return x
    
def parse_dipole_moment(x):
    """
    ORCA parsing for Dipole Moment
    """
    if re.search('Total Dipole Moment', x):
        return x
    
def parse_single_point_energy(x):
    """
    ORCA parsing for single point energy (SCF + Dispersion + Solvation)
    """
    if re.search('FINAL SINGLE POINT ENERGY', x):
        return x


def get_freq(fl):
    """
    ORCA parsing for frequencies
    """

    fl = '\n'.join(''.join(fl).split('VIBRATIONAL FREQUENCIES')[-1].strip().splitlines()[4:]).split('------------')[0].strip().splitlines()
    freq = np.array([float(i.split()[1]) for i in fl if float(i.split()[1]) != 0.0 ])
    return freq


def get_conf_parameters(conf, number, p, time, temp, log):

    # data = cclib.io.ccread(os.path.join(conf.folder, f'protocol_{number}.out'))

    with open(os.path.join(conf.folder, f'protocol_{number}.out')) as f:
        fl = f.readlines() 

    e = float(list(filter(parse_single_point_energy, fl))[-1].strip().split()[-1])    

    freq = np.array([])
    if p.freq:
        freq = get_freq(fl)
        if freq.size == 0:
            log.error(('\n'.join(fl[-6:])).strip())
            raise RuntimeError('Some sort of error have been encountered during the calculation of the calcultor.')
    
    B = np.array(list(filter(parse_rotational_const, fl))[-1].strip().split(':')[-1].split(), dtype=float)
    b = np.linalg.norm(B)

    M = np.linalg.norm(np.array(list(filter(parse_dipole_moment, fl))[-1].strip().split(':')[-1].split(), dtype=float))

    g = ''
    if freq.size > 0:
        g = free_gibbs_energy(
            SCF = e, T = temp, freq=freq, 
            mw=conf.weight_mass, 
            B = B,
            m=conf.mult
            )
        
        print(g)


    conf.energies[int(number)] = {
        'E' : e * EH_TO_KCAL if e else e,    #   Electronic Energy [kcal/mol]
        'G' : g * EH_TO_KCAL if g else None, #   Free Gibbs Energy [kcal/mol]
        'B' : b if b else None,                 #   Rotatory Constant [cm-1]
        'm' : M if M else None,                 #   dipole momenti [Debye]
        'time' : time,                          #   elapsed time [sec] 
    }

    return None




if __name__ == '__main__':

    # from logger import create_log

    # class Conf:
    #     def __init__(self, number, mult, folder):
    #         self.number = number
    #         self.mult = mult
    #         self.folder = folder
    #         self.energies = {}

    # c = Conf('1', 1, 'conf_1')
    # log = create_log('test.out')

    # get_conf_parameters(c, 0, 1, 298.15, log)
    # print(c.energies)

    with open('files/opt.out') as f:
        fl = f.readlines()
    get_freq(fl)