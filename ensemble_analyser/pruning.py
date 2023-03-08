from rmsd import quaternion_rmsd as rmsd
import numpy as np


def cut_over_thr_max(confs: list, thrGMAX: float) -> list:
    """
    confs is a sorted list of conformers

    """

    ens = np.array([i.get_energy for i in confs])
    ens = (ens-min(ens))*627.51
    remov_confs = np.array(confs)[ens > thrGMAX]
    for i in remov_confs:
        i.active = False
        print(i.number, i.active, 'cut_over')


def check(check, conf_ref, thrs) -> None:
    if (check.get_energy - conf_ref.get_energy < thrs['thrG'] and check.rotatory - conf_ref.rotatory < thrs['thrB']):
        check.active = False
        print(check.number, check.active, 'check')

    return None



def check_ensemble(confs, thrs) -> list:

    cut_over_thr_max(confs, thrs['thrGMAX'])

    for idx, i in enumerate(confs):
        if not i.active: continue # Not check the non active conformers

        for j in range(0, idx):
            check(i, confs[j], thrs)
            # print(rmsd(i.last_geometry, confs[j].last_geometry))

        print(i.number, i.energies)


    return confs



