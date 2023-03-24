from ensemble_analyser.logger import ordinal
from ensemble_analyser.parser_parameter import get_conf_parameters


import ase
import time
import sys, os


def launch(conf, protocol, cpu, log):

    log.info(f'Running {ordinal(int(protocol.number))} PROTOCOL -> CONF{conf.number}')
    try:
        st = time.perf_counter()

        calculator = protocol.get_calculator(cpu=cpu, charge=conf.charge, mult=conf.mult)
        atm = conf.get_ase_atoms(calculator)
        atm.get_potential_energy()

        end = time.perf_counter()

        os.rename(f'{label}.out', f'{conf.folder}/protocol_{protocol.number}.out')


    except ase.calculators.calculator.CalculationFailed:
        log.error(f'Calulator error.')
        with open('orca.out') as f:
            fl = f.read()
        log.error('\n'.join(fl.splitlines()[-6:-3]))
        raise RuntimeError('Some sort of error have been encountered during the calculation of the calcultor.')


    get_conf_parameters(conf, protocol.number, end-st)