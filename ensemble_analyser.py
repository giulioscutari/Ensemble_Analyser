#!/usr/bin/python3

import ase
from ase.optimize import LBFGS


import sys, os, cclib
from ensemble_analyser.protocol import create_protocol

from ensemble_analyser.pruning import check_ensemble
from ensemble_analyser.ioFile import read_ensemble
from ensemble_analyser.parser_parameter import get_conf_parameters
from ensemble_analyser.logger import save_snapshot, create_log, ordinal
from ensemble_analyser.protocol import load_threshold, load_protocol
from ensemble_analyser.parser_arguments import parser_arguments




def launch(conf, protocol, cpu, log):
    
    log.info(f'Running PROTOCOL {ordinal(int(protocol.number))} CONF{conf.number}')
    try:
        if protocol.opt:
            calculator = protocol.get_calculator(cpu=cpu, opt=True)
            atm = conf.get_ase_atoms(calculator)
            opt = LBFGS(atm, logfile='ORCA_ase.log', trajectory='ORCA_ase.trj')
            opt.run()
            conf.last_geometry = opt.atoms.get_positions()

        if protocol.freq:
            calculator = protocol.get_calculator(cpu=cpu)
            atm = conf.get_ase_atoms(calculator)
            atm.get_potential_energy()

        if not (protocol.opt and protocol.freq):
            calculator = protocol.get_calculator(cpu=cpu)
            atm = conf.get_ase_atoms(calculator)
            atm.get_potential_energy()

    except ase.calculators.calculator.CalculationFailed:
        log.error(f'Calulator error.')
        with open('ORCA.out') as f:
            fl = f.read()
        log.error('\n'.join(fl.splitlines()[-6:-3]))
        raise RuntimeError('Some sort of error have been encountered during the calculation of the calcultor.')


    get_conf_parameters(conf, protocol.number)
    

def main(ensemble: str, protocol_file: str | None, threshold_file: str | None, cpu:int, output):

    log = create_log(output)

    log.debug('# Reading protocol and threshold files')
    thrs = load_threshold(threshold_file)
    protocol = load_protocol(protocol_file)
    protocol = create_protocol(protocol, thrs, log)    


    log.debug('# Reading ensemble file')
    confs = read_ensemble(ensemble)

    for p in protocol:

        log.info(f'STARTING PROTOCOL {p.number}')

        for i in confs:
            if not i.active: continue
            launch(i, p, cpu, log)

        confs = sorted(confs)

        log.info(f'\n{"CONFS":<10s}\t{"E[kcal/mol]":<10s} \t{"G[kcal/mol]":<10s} \t{"B[cm-1]":<10s}')
        for i in confs:
            if not i.active: continue
            # log.debug('Creating the OUTPUT for each conformer calculated')
            log.info(i.create_log())

        confs = check_ensemble(confs, p, log)
        save_snapshot(f'ensemble_after_{p.number}.xyz', confs, log)

        log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')

    save_snapshot('final_ensemble.xyz', confs, log)

    return None


if __name__ == '__main__':
    
    args = parser_arguments()

    main(
        ensemble = args.ensemble, 
        protocol_file = args.protocol, 
        threshold_file = args.threshold,
        cpu = args.cpu,
        output = args.output
    ) 


    os.system('rm ORCA*')