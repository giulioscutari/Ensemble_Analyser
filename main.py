from openbabel import openbabel as ob
from ase.optimize import LBFGS

from ensemble_analyser.pruning import check_ensemble
from ensemble_analyser.calculator import get_calculator
from ensemble_analyser.ioFile import read_ensemble
from ensemble_analyser.parser_parameter import get_conf_parameters
from ensemble_analyser.logger import save_snapshot
from ensemble_analyser.protocol import load_protocol, load_threshold



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
    


def main(ensemble: str, protocol_file: str | None, threshold_file: str | None):

    protocol = load_protocol(protocol_file)
    thrs = load_threshold(threshold_file)

    confs = read_ensemble(ensemble)

    for idx, p in protocol.items():
        for i in confs:
            if not i.active: continue
            launch(idx, i, p)

        confs = sorted(confs)
        confs = check_ensemble(confs, thrs)
        save_snapshot(f'ensemble_after_{idx}.xyz', confs)

    save_snapshot('final_ensemble.xyz', confs)

    return None


if __name__ == '__main__':
    
    main(ensemble = 'tests/struct/cluster_water.xyz', protocol_file = None, threshold_file = None) 


    # os.system('rm ORCA*')