from ensemble_analyser.logger import ordinal
from ensemble_analyser.parser_parameter import get_conf_parameters


import ase
import time
import sys


def launch(conf, protocol, cpu, log):

    log.info(f'Running {ordinal(int(protocol.number))} PROTOCOL -> CONF{conf.number}')
    try:
        st = time.perf_counter()
        if protocol.opt and protocol.freq:
            calculator = protocol.get_calculator(cpu=cpu, opt=True)
            atm = conf.get_ase_atoms(calculator)
            atm.get_potential_energy()
            conf.last_geometry = atm.get_positions()

        calculator = protocol.get_calculator(cpu=cpu)
        atm = conf.get_ase_atoms(calculator)
        atm.get_potential_energy()

        end = time.perf_counter()

    except ase.calculators.calculator.CalculationFailed:
        log.error(f'Calulator error.')
        with open('ORCA.out') as f:
            fl = f.read()
        log.error('\n'.join(fl.splitlines()[-6:-3]))
        raise RuntimeError('Some sort of error have been encountered during the calculation of the calcultor.')


    get_conf_parameters(conf, protocol.number, end-st)