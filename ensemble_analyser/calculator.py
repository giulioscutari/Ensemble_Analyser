from ase.calculators.orca import ORCA
import os


DEBUG = os.getenv('DEBUG')




def get_calculator(func:str, basis:str, opt:bool=True, freq:bool=False, cpu:int =11, charge:int =0, mult:int =1,):
    calculator = ORCA(
        label = "ORCA",
        orcasimpleinput=f'{func} {basis} {"opt" if opt else ""} {"freq" if freq else ""} SMALLPRINT',
        orcablocks=f'%pal nprocs {cpu} end' if not DEBUG else '',
        charge = charge, mult=mult, task='gradient'
    )

    return calculator
