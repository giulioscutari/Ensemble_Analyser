from rmsd import quaternion_rmsd as rmsd
import numpy as np
from .logger import log

def cut_over_thr_max(confs: list, thrGMAX: float) -> list:


    ens = np.array([i.get_energy for i in confs])
    ens = (ens-min(ens))*627.51
    print(ens)
    remov_confs = np.array(confs)[ens > thrGMAX]
    log.info(f'\nGetting number of conformers lying out of the energy windows (over {thrGMAX} kcal/mol) - {len(remov_confs)}')
    for idx, i in enumerate(list(remov_confs)):
        i.active = False
        log.info(f'{i.number} - {ens[confs.index(i)]}')
    log.info('\n\n')

def check(check, conf_ref, protocol) -> None:

    if (check.get_energy - conf_ref.get_energy < (protocol.thrG / 627.51) and check.rotatory - conf_ref.rotatory < protocol.thrB):
        check.active = False
        check.diactivated_by = conf_ref.number
        log.info(f'{check.number} deactivated by {conf_ref.number}. {check.moment:.4f} - {conf_ref.moment:.4f}')

    return None



def check_ensemble(confs, protocol) -> list:

    cut_over_thr_max(confs, protocol.thrGMAX)

    for idx, i in enumerate(confs):
        if not i.active: continue # Not check the non active conformers

        for j in range(0, idx):
            check(i, confs[j], protocol)
            # rmsd(i.last_geometry, confs[j].last_geometry)

    log.debug('\n'.join([f'{i.number}: {i.energies} -- {i.active}' for i in confs]))


    return confs



