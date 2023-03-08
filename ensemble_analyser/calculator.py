from ase.calculators.orca import ORCA
import os


DEBUG = os.getenv('DEBUG')




def get_calculator(protocol, opt=False, cpu:int =11, charge:int =0, mult:int =1):

    # Input
    func = protocol['func']
    basis = protocol['basis']
    s = protocol['solv']['solv']
    smd = f'%cpcm smd true smdsolvent "{s}" end' if protocol['solv']['smd'] else ''

    freq = protocol['freq'] if not opt else False

    calculator = ORCA(
        label = "ORCA",
        orcasimpleinput=f'{func} {basis} {"freq" if freq else ""} {"cpcm("+s+")" if s else ""}',
        orcablocks=(f'%pal nprocs {cpu} end ' if not DEBUG else '') + smd,
        charge = charge, mult=mult, task='gradient'
    )

    return calculator
