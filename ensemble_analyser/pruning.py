#!/data/bin/python_env/bin/python3

from ase.build import minimize_rotation_and_translation
from scipy.constants import R
from tabulate import tabulate
import numpy as np


from .logger import save_snapshot, DEBUG





def cut_over_thr_max(confs: list, thrGMAX: float, log) -> list:


    ens = np.array([(i, i.get_energy) for i in confs if i.active])
    ens[:,1] = ens[:,1]-min(ens[:,1])

    log.info(f'\nGetting number of conformers lying out of the energy windows (over {thrGMAX} kcal/mol)')
    for i, en in list(ens):
        if en > thrGMAX:
            i.active = False
            log.info(f'{i.number} - {en:.3f}')
    log.info('\n')



def rmsd(check, ref):
    ref_pos, check_pos = ref.copy(), check.copy()
    minimize_rotation_and_translation(ref_pos, check_pos)
    return np.sqrt(1/len(ref.get_positions())) * np.linalg.norm(np.array(ref_pos.get_positions())-np.array(check_pos.get_positions()))


def dict_compare(check, conf_ref, deactivate=True):
    return {
        'Check'         : check.number,
        'Ref'           : conf_ref.number,
        '∆E [kcal/mol]' : check.get_energy - conf_ref.get_energy,
        '∆B [e-3 cm-1]' : np.abs(check.rotatory - conf_ref.rotatory)*10**3,
        '∆m [Debye]'    : np.abs(check.moment - conf_ref.moment),
        'RMSD [Å]'      : rmsd(check.get_ase_atoms(), conf_ref.get_ase_atoms()),
        'Deactivate'    : deactivate
    }



def check(check, conf_ref, protocol, controller, log) -> bool:

    if not conf_ref.active:
        return False

    l = len(controller)

    controller[l] = dict_compare(check, conf_ref, deactivate=False)

    if ( controller[l]['∆E [kcal/mol]'] < protocol.thrG  and controller[l]['∆B [e-3 cm-1]']*10**-3 < protocol.thrB):
        check.active = False
        check.diactivated_by = conf_ref.number
        controller[l] = dict_compare(check, conf_ref)
        return True

    controller.pop(l)

    return False




def refactor_dict(controller, log):
    
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
            if check(i, confs[j], protocol, controller, log): 
                break
    
    controller = refactor_dict(controller, log)

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




