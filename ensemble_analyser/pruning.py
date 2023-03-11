#!/usr/bin/python3

# from rmsd import quaternion_rmsd as rmsd
import numpy as np

def cut_over_thr_max(confs: list, thrGMAX: float, log) -> list:


    ens = np.array([i.get_energy for i in confs if i.active])
    ens = ens-min(ens)

    remov_confs = np.array([i for i in confs if i.active])[ens > thrGMAX]
    log.info(f'\nGetting number of conformers lying out of the energy windows (over {thrGMAX} kcal/mol) - {len(remov_confs)}')
    for idx, i in enumerate(list(remov_confs)):
        i.active = False
        log.info(f'{i.number} - {ens[confs.index(i)]}')
    log.info('\n\n')

def check(check, conf_ref, protocol, log) -> None:


    log.debug(f'{check.number} VS {conf_ref.number}: ∆Energy = {check.get_energy - conf_ref.get_energy} - ∆B = {check.rotatory - conf_ref.rotatory}')
    if (check.get_energy - conf_ref.get_energy < protocol.thrG  and check.rotatory - conf_ref.rotatory < protocol.thrB):
        check.active = False
        check.diactivated_by = conf_ref.number
        log.info(f'{check.number} deactivated by {conf_ref.number}. ENERGY\t{check.get_energy:.4f} - {conf_ref.get_energy:.4f}\tDIPOL MOMENTS\t{check.moment:.4f} - {conf_ref.moment:.4f}')
        return True

    return False



def check_ensemble(confs, protocol, log) -> list:

    cut_over_thr_max(confs, protocol.thrGMAX, log)

    for idx, i in enumerate(confs):
        if not i.active: continue # Not check the non active conformers
        for j in range(0, idx):
            print('check ', j, idx)
            if check(i, confs[j], protocol, log): 
                break
            # rmsd(i.last_geometry, confs[j].last_geometry)

    log.debug('\n'.join([f'{i.number}: {i._last_energy} -- {i.active}' for i in confs]))


    return confs



