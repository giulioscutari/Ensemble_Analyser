#!/data/bin/python_env/bin/python3

import os, datetime
from tabulate import tabulate



from ensemble_analyser.protocol import create_protocol
from ensemble_analyser.pruning import check_ensemble, calculate_rel_energies
from ensemble_analyser.ioFile import read_ensemble
from ensemble_analyser.logger import save_snapshot, create_log
from ensemble_analyser.protocol import load_threshold, load_protocol
from ensemble_analyser.parser_arguments import parser_arguments
from ensemble_analyser.launch import launch



def main(ensemble: str, protocol_file: str , threshold_file: str , cpu:int, output, charge: int, multiplicity: int, temperature: float):

    log = create_log(output)

    log.debug('# Reading protocol and threshold files')
    thresholds = load_threshold(threshold_file)
    protocol = load_protocol(protocol_file)
    protocol = create_protocol(protocol, thresholds, log)    


    log.debug('# Reading ensemble file')
    conformers = read_ensemble(ensemble, charge, multiplicity, log)

    return 

    for p in protocol:

        log.info(f'STARTING PROTOCOL {p.number}')
        log.info(f'\nActive conformers for this phase: {len([i for i in conformers if i.active])}\n')
        for i in conformers:
            if not i.active: continue
            launch(i, p, cpu, log)

        conformers = sorted(conformers)

        calculate_rel_energies(conformers, temperature)

        log.info('')
        log.info('Ended Calculations')
        log.info('')
        log.info('Summary')
        log.info('')
        log.info(tabulate([i.create_log() for i in conformers if i.active], headers=['conformers', 'E[Eh]' ,'G[Eh]', 'B[cm-1]', 'E. Rel [kcal/mol]', 'Pop [%]', 'Elap. time [sec]'], floatfmt=".6f"))
        log.info('')
        log.info('Total elapsed time: ' + str(datetime.timedelta(seconds = sum([i._last_energy['time'] for i in conformers if i.active]))))

        log.info('Start Pruning')
        conformers = check_ensemble(conformers, p, log)
        save_snapshot(f'ensemble_after_{p.number}.xyz', conformers, log)

        log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')

    save_snapshot('final_ensemble.xyz', conformers, log)
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