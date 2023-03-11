#!/usr/bin/python3

import json, os
import sys
from ase.calculators.orca import ORCA

try:
    from .logger import ordinal
except ModuleNotFoundError:
    pass


DEBUG = os.getenv('DEBUG')


def load_protocol(file:str ):
    default = 'ensemble_analyser/parameters_file/default_protocol.json'
    return json.load(open(default if not file else file ))

def load_threshold(file:str ):
    default = 'ensemble_analyser/parameters_file/default_threshold.json'
    return json.load(open(default if not file else file ))




LEVEL_DEFINITION = {
    0 : 'SP'.lower()          , # mere energy calculation
    1 : 'OPT'.lower()         , # optimisation step
    2 : 'FREQ'.lower()        , # single point and frequency analysis
    3 : 'OPT+FREQ'.lower()      # optimisation and frequency analysis
}




class Solvent: 
    def __init__(self, solv:dict):
        self.solvent = solv['solvent']
        self.smd = solv['smd']

    def __str__(self):
        return f' - SMD({self.solvent})' if self.smd else  f' - CPCM({self.solvent})'
    def __repr__(self):
        return f' - SMD({self.solvent})' if self.smd else  f' - CPCM({self.solvent})'

    def orca_input_smd(self):
        if self.smd:
            return f'%cpcm smd true smdsolvent "{self.solvent}" end'
        return ''
    


class Protocol: 

    def __init__(self, number:int, functional:str, basis:str, solvent:Solvent, opt:bool, freq:bool, add_input:str , thrs_json, calculator='orca'):

        self.number = number
        self.functional = functional.upper()
        self.basis = basis.upper()
        self.solvent = solvent
        self.opt = opt
        self.freq = freq
        self.add_input = add_input
        self.get_thrs(thrs_json)
        self.calculator = calculator

    @property
    def calculation_level(self):
        return LEVEL_DEFINITION[self.number_level].upper()

    @property
    def level(self):
        return f'{self.functional}/{self.basis}' + (str(self.solvent) if self.solvent else '')
    
    @property
    def thr(self):
        return f'\tthrG    : {self.thrG} kcal/mol\n\tthrB    : {self.thrB} cm-1\n\tthrGMAX : {self.thrGMAX} kcal/mol\n'
    
    @property
    def number_level(self):
        c = 0
        if self.opt: c += 1
        if self.freq: c += 2
        return c

    def get_calculator(self, cpu, opt=False):
        calc = {
            'orca' : self.get_orca_calculator(cpu, opt)
        }

        return calc[self.calculator]
        

    def get_thrs(self, thr_json):
        c = LEVEL_DEFINITION[self.number_level]
        self.thrG = thr_json[c]['thrG']
        self.thrB = thr_json[c]['thrB']
        self.thrGMAX = thr_json[c]['thrGMAX']


    def __str__(self):
        return f'{self.functional}/{self.basis}' + (str(self.solvent) if self.solvent.solvent else '')
    def __repr__(self):
        return f'{self.functional}/{self.basis}' + (str(self.solvent) if self.solvent.solvent else '')
    

    def get_orca_calculator(self, cpu:int, opt:bool = False, charge:int = 0, mult:int = 1):
        
        # if FREQ and OPT, every calculation from ASE is also an HESSIAN calculation, so switching off FREQ when also optimizing
        freq = self.freq if not opt else False

        # possibilities for solvent definitions
        if self.solvent.solvent:
            if 'xtb' in self.functional.lower():
                solv = f"ALPB({self.solvent.solvent})"
            else:
                solv = f" cpcm({self.solvent.solvent})"
        else:
            solv = ''

        # ! B3LYP def2-SVP FREQ CPCM(solvent) ENGRAD
        # the optimisation is carried by ASE, ORCA is a gradient engine
        simple_input = f'{self.functional} {self.basis} {"freq" if freq else ""}{"opt" if opt else ""}{solv}'


        # %cpcm
        #     smd True
        #     smdSolvent solvent
        # end
        smd = ''
        if self.solvent and 'xtb' not in self.functional.lower(): smd = self.solvent.orca_input_smd()

        calculator = ORCA(
        label = "ORCA",
        orcasimpleinput = simple_input,
        orcablocks=f'%pal nprocs {cpu} end ' + smd,
        charge = charge, mult=mult, task='energy'
        )


        # log.debug(f'# Parameters of calculator: {calculator.parameters}\n')

        return calculator



def create_protocol(p, thrs, log):
    protocol = []

    log.info('Loading Protocol\n')
    for idx, d in p.items():
        func    = d.get('func', None)
        basis   = d.get('basis', 'def2-svp')
        opt     = d.get('opt', False)
        freq    = d.get('freq', False)

        add_input = d.get('add_input', None)

        solv    = d.get('solv', None)
        if solv or solv.get('solvent', None): solv = Solvent(solv)

        if not func:
            log.critical(f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nFUNC key must be passed in order to calculate energy. DFT functional or HF for Hartree-Fock calculation or semi-empirical methods (XTB1/XTB2/PM3/AM1 or similar supported by the calculator) (Probelm at {ordinal(int(idx))} protocol definition)\n{'='*20}\n")
            raise IOError('There is an error in the input file with the definition of the functional. See the output file.')

        protocol.append(Protocol(number =idx, functional=func, basis=basis, solvent= solv, opt= opt, freq= freq, add_input= add_input, thrs_json= thrs))

    log.info('\n'.join((f"{i.number}: {str(i)} - {i.calculation_level}\n {i.thr}" for i in protocol)) + '\n')
    return protocol
    #{'charge': 0, 'mult': 1, 'task': 'gradient', 'orcasimpleinput': 'B3LYP def2-svp ', 'orcablocks': ''}


if __name__ == '__main__':
    from logger import log

    thrs_json = load_threshold(None)

    p = Protocol(1, 'B3LYP', 'def2-svp', None, True, False, None, thrs_json)
    log.debug(f'''{p} THRS:
thrG    : {p.thrG} kcal/mol
thrB    : {p.thrB} cm-1
thrGMAX : {p.thrGMAX} kcal/mol
    ''')

    p.get_orca_calculator(1)

