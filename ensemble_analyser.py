#!/usr/bin/python3

import os
from tabulate import tabulate

import numpy as np
from scipy.constants import R


from ensemble_analyser.protocol import create_protocol
from ensemble_analyser.pruning import check_ensemble
from ensemble_analyser.ioFile import read_ensemble
from ensemble_analyser.logger import save_snapshot, create_log
from ensemble_analyser.protocol import load_threshold, load_protocol
from ensemble_analyser.parser_arguments import parser_arguments
from ensemble_analyser.launch import launch



def calculate_rel_energies(confs, T):
    c = [i for i in confs if i.active]
    ens = np.array([i.get_energy for i in confs if i.active])
    ens -= min(ens)
    bolz = np.exp((-ens*4186)/(R*T))
    pop = (bolz/np.sum(bolz))*100
    for idx, i in enumerate(list(ens)):
        c[idx]._last_energy['Erel'] = i 
        c[idx]._last_energy['Pop'] = pop[idx]



def main(ensemble: str, protocol_file: str , threshold_file: str , cpu:int, output, charge: int, multiplicity: int, temperature: float):

    log = create_log(output)

    log.debug('# Reading protocol and threshold files')
    thrs = load_threshold(threshold_file)
    protocol = load_protocol(protocol_file)
    protocol = create_protocol(protocol, thrs, log)    


    log.debug('# Reading ensemble file')
    confs = read_ensemble(ensemble, charge, multiplicity)

    for p in protocol:

        log.info(f'STARTING PROTOCOL {p.number}')
        log.info(f'\nActive conformers for this phase: {len([i for i in confs if i.active])}\n')
        for i in confs:
            if not i.active: continue
            launch(i, p, cpu, log)

        confs = sorted(confs)

        calculate_rel_energies(confs, temperature)

        log.info(tabulate([i.create_log() for i in confs if i.active], headers=['CONFS', 'E[Eh]' ,'G[Eh]', 'B[cm-1]', 'E. Rel [kcal/mol]', 'Pop [%]', 'Elap. time [sec]'], floatfmt=".6f"))

        confs = check_ensemble(confs, p, log)
        save_snapshot(f'ensemble_after_{p.number}.xyz', confs, log)

        log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')

    save_snapshot('final_ensemble.xyz', confs, log)
    log.info(f'{"="*15}\nCALCULATIONS ENDED\n{"="*15}\n\n')

    return None


if __name__ == '__main__':
    
    args = parser_arguments()

    main(
        ensemble = args.ensemble, 
        protocol_file = args.protocol, 
        threshold_file = args.threshold,
        cpu = args.cpu,
        output = args.output,

        charge=args.charge, 
        multiplicity=args.multiplicity, 
        temperature=args.temperature, 
    ) 


    os.system('rm ORCA*')