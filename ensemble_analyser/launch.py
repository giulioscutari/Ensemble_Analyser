

from ensemble_analyser.conformer import Conformer
from ensemble_analyser.ioFile import read_ensemble
from ensemble_analyser.logger import create_log, ordinal, save_snapshot
from ensemble_analyser.parser_arguments import parser_arguments
from ensemble_analyser.parser_parameter import get_conf_parameters
from ensemble_analyser.IOsystem import SerialiseEncoder


import ase
import time, json
import datetime
from tabulate import tabulate
import sys, os

from ensemble_analyser.protocol import Protocol, Solvent, load_protocol, load_threshold
from ensemble_analyser.pruning import calculate_rel_energies, check_ensemble


def launch(conf, protocol, cpu, log, ensemble):

    log.info(f'Running {ordinal(int(protocol.number))} PROTOCOL -> CONF{conf.number}')
    try:
        st = time.perf_counter()

        calculator, label = protocol.get_calculator(cpu=cpu, charge=conf.charge, mult=conf.mult)
        atm = conf.get_ase_atoms(calculator)
        try:
            atm.get_potential_energy()
        except ase.calculators.calculator.PropertyNotImplementedError:
            pass

        end = time.perf_counter()

        os.rename(f'{label}.out', f'{conf.folder}/protocol_{protocol.number}.out')
        os.remove(f'{label}.gbw')


    except ase.calculators.calculator.CalculationFailed:
        log.error(f'Calulator error.')
        with open(f'{label}.out') as f:
            fl = f.read()
        log.error('\n'.join(fl.splitlines()[-6:-3]))
        raise RuntimeError('Some sort of error have been encountered during the calculation of the calcultor.')


    get_conf_parameters(conf, protocol.number, end-st, log)

    json.dump({i.number: i.__dict__ for i in ensemble}, open('checkpoint.json', 'w'), indent=4, cls=SerialiseEncoder)



def run_protocol(conformers, p, temperature, cpu, log):
    log.info(f'STARTING PROTOCOL {p.number}')
    log.info(f'\nActive conformers for this phase: {len([i for i in conformers if i.active])}\n')
    for i in conformers:
        if not i.active: continue
        if i.energies.get(str(p.number)): continue
        launch(i, p, cpu, log, conformers)

    conformers_tmp = sorted(conformers)

    calculate_rel_energies(conformers, temperature)

    log.info('')
    log.info('Ended Calculations')
    log.info('')
    log.info('Summary')
    log.info('')
    log.info(tabulate([i.create_log() for i in conformers_tmp if i.active], headers=['conformers', 'E[Eh]' ,'G[Eh]', 'B[cm-1]', 'E. Rel [kcal/mol]', 'Pop [%]', 'Elap. time [sec]'], floatfmt=".6f"))
    log.info('')
    log.info('Total elapsed time: ' + str(datetime.timedelta(seconds = sum([i._last_energy['time'] for i in conformers_tmp if i.active]))))

    log.info('Start Pruning')
    conformers = check_ensemble(conformers, p, log)
    save_snapshot(f'ensemble_after_{p.number}.xyz', conformers, log)

    log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')


def start_calculation(conformers, protocol, cpu:int, temperature: float, start_from: int, log):

    log.debug('# Reading protocol and threshold files')


    log.debug('# Reading ensemble file')

    for p in protocol[start_from:]:
        with open('last_protocol', 'w') as f:
            f.write(str(p.number))

        run_protocol(conformers, p, temperature, cpu, log)

    save_snapshot('final_ensemble.xyz', conformers, log)
    log.info(f'{"="*15}\nCALCULATIONS ENDED\n{"="*15}\n\n')

    return None



def restart():

    confs = json.load(open('checkpoint.json'))
    ensemble = [Conformer.load_raw(confs[i]) for i in confs]

    p = json.load(open('protocol_dump.json'))
    protocol = [Protocol.load_raw(p[i]) for i in p]

    with open('last_protocol') as f:
        start_from = int(f.readlines()[0])

    return ensemble, protocol, start_from



def create_protocol(p, thrs, log):
    protocol = []

    log.info('Loading Protocol\n')
    for idx, d in p.items():
        func    = d.get('func', None)
        basis   = d.get('basis', 'def2-svp')
        opt     = d.get('opt', False)
        freq    = d.get('freq', False)

        add_input = d.get('add_input', '')

        solv    = d.get('solv', None)
        if solv or solv.get('solvent', None): solv = Solvent(solv)

        thrG    = d.get('thrG', None)
        thrB    = d.get('thrB', None)
        thrGMAX = d.get('thrGMAX', None)

        if not func:
            log.critical(f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nFUNC key must be passed in order to calculate energy. DFT functional or HF for Hartree-Fock calculation or semi-empirical methods (XTB1/XTB2/PM3/AM1 or similar supported by the calculator) (Probelm at {ordinal(int(idx))} protocol definition)\n{'='*20}\n")
            raise IOError('There is an error in the input file with the definition of the functional. See the output file.')

        protocol.append(Protocol(number =idx, functional=func, basis=basis, solvent= solv, opt= opt, freq= freq, add_input = add_input, thrs_json= thrs, thrG = thrG, thrB = thrB, thrGMAX = thrGMAX))

    log.info('\n'.join((f"{i.number}: {str(i)} - {i.calculation_level}\n {i.thr}" for i in protocol)) + '\n')
    return protocol

def main():
    args = parser_arguments()

    try:
        settings = json.load(open('settings.json'))
    except FileNotFoundError:
        settings = {
            'output' : args.output,
            'cpu' : args.cpu,
            'temperature' : args.temperature,
        }
        json.dump(settings, open('settings.json', 'w'), indent=4)

    output = settings.get('output', args.output) if not args.restart else ''.join(settings.get('output', args.output).split('.')[:-1])+'_restart.out'
    cpu = settings.get('cpu', args.cpu)
    temperature = settings.get('temperature', args.temperature)


    log = create_log(output)

    if args.restart:
        conformers, protocol, start_from = restart()
    else:
        protocol = create_protocol(load_protocol(args.protocol), load_threshold(args.threshold), log)
        start_from = protocol[0].number
        json.dump({i.number: i.__dict__ for i in protocol}, open('protocol_dump.json', 'w'), indent=4, cls=SerialiseEncoder)
        conformers = read_ensemble(args.ensemble, args.charge, args.multiplicity, log)

    start_calculation(
        conformers = conformers,
        protocol = protocol,
        cpu = cpu,

        temperature= temperature,
        start_from= int(start_from),
        log = log,
    )


    os.system('rm ORCA*')






