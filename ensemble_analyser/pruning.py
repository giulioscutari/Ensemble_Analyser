#!/usr/bin/python3

from ase.build import minimize_rotation_and_translation
from scipy.constants import R
from tabulate import tabulate
import numpy as np


from .logger import save_snapshot, DEBUG





def cut_over_thr_max(confs: list, thrGMAX: float, log) -> list:


    ens = np.array([i.get_energy for i in confs if i.active])
    ens = ens-min(ens)

    remov_confs = np.array([i for i in confs if i.active])[ens > thrGMAX]
    log.info(f'\nGetting number of conformers lying out of the energy windows (over {thrGMAX} kcal/mol) - {len(remov_confs)}')
    for i in list(remov_confs):
        i.active = False
        log.info(f'{i.number} - {ens[confs.index(i)]}')
    log.info('\n')



def rmsd(check, ref):
    ref_pos, check_pos = ref.copy(), check.copy()
    minimize_rotation_and_translation(ref_pos, check_pos)
    return np.sqrt(1/len(ref.get_positions())) * np.linalg.norm(np.array(ref_pos.get_positions())-np.array(check_pos.get_positions()))




def check(check, conf_ref, protocol, controller) -> bool:

    if not conf_ref.active:
        return False

    l = len(controller)
    controller[l] = {
        'Check'         : check.number,
        'Ref'           : conf_ref.number,
        '∆E [kcal/mol]' : check.get_energy - conf_ref.get_energy,
        '∆B [e-3 cm-1]' : np.abs(check.rotatory - conf_ref.rotatory)*10**3,
        '∆m [Debye]'    : np.abs(check.moment - conf_ref.moment),
        'RMSD [Å]'      : rmsd(check.get_ase_atoms(), conf_ref.get_ase_atoms()),
        'Deactivate'    : False
    }

    if ( controller[l]['∆E [kcal/mol]'] < protocol.thrG  and controller[l]['∆B [e-3 cm-1]']*10**-3 < protocol.thrB):
        check.active = False
        check.diactivated_by = conf_ref.number
        controller[l] = True
        return True

    return False




def refactor_dict(controller):
    
    if not controller: return {}

    keys = list(controller[0].keys())
    d = {i : [] for i in keys}

    for i in controller:
        for j in d:
            d[j].append(controller[i][j])
    return d


def check_ensemble(confs, protocol, log) -> list:

    cut_over_thr_max(confs, protocol.thrGMAX, log)

    if DEBUG: save_snapshot(f'after_protocol_{protocol.number}_before_check.xyz', confs, log)

    controller = {}

    for idx, i in enumerate(confs):
        if not i.active: continue # Not check the non active conformers
        for j in range(0, idx):
            print('check ', j, idx)
            if check(i, confs[j], protocol, controller): 
                break
    
    
    controller = refactor_dict(controller)

    log.info('')
    log.info(tabulate(controller, headers="keys", floatfmt=".3f"))
    log.info('')

    return confs

def calculate_rel_energies(conformers, T):
    c = [i for i in conformers if i.active]
    ens = np.array([i.get_energy for i in conformers if i.active])
    ens -= min(ens)
    bolz = np.exp((-ens*4186)/(R*T))
    pop = (bolz/np.sum(bolz))*100
    for idx, i in enumerate(list(ens)):
        c[idx]._last_energy['Erel'] = i
        c[idx]._last_energy['Pop'] = pop[idx]




