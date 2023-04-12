

from ensemble_analyser.conformer import Conformer
from ensemble_analyser.ioFile import read_ensemble
from ensemble_analyser.logger import create_log, ordinal, save_snapshot
from ensemble_analyser.parser_arguments import parser_arguments
from ensemble_analyser.parser_parameter import get_conf_parameters
from ensemble_analyser.IOsystem import SerialiseEncoder


import ase
import time, json
import sys, os

from ensemble_analyser.protocol import Protocol, create_protocol, load_protocol, load_threshold, run_protocol


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


    get_conf_parameters(conf, protocol.number, end-st)

    json.dump({i.number: i.__dict__ for i in ensemble}, open('checkpoint.json', 'w'), indent=4, cls=SerialiseEncoder)

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
        json.dumps(settings, open('settings.json', 'w'), indent=4)

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




