#!/data/bin/python_env/bin/python3

import numpy as np
from ase.atoms import Atoms
from ensemble_analyser.IOsystem import mkdir

class Conformer:

    def __init__(self, number: int, geom: np.array, atoms: np.array, charge : int= 0, mult : int = 1) -> None:
        self.number = number
        self._initial_geometry = geom
        self.charge = charge
        self.mult = mult
        
        self.last_geometry = geom
        self.atoms = atoms
        self.energies = {}
        self.active = True
        
        # IO
        self.folder = f'conf_{self.number}'
        mkdir(self.folder)
    
    def get_ase_atoms(self, calc=None):
        return Atoms(
            symbols = ''.join(list(self.atoms)),
            positions = self.last_geometry,
            calculator=calc
        )

    @property
    def rotatory(self):
        return self.energies[list(self.energies.keys())[-1]]['B']
    
    @property
    def moment(self):
        return self.energies[list(self.energies.keys())[-1]]['m']

    @property
    def get_energy(self):
        en = self.energies[list(self.energies.keys())[-1]]
        if en['G']: return en['G']
        return en['E']
    
    @property
    def _last_energy(self):
        return self.energies[list(self.energies.keys())[-1]]
        
    def write_xyz(self):
        if not self.active: return ''
        txt = f'{len(self.atoms)}\nCONFORMER {self.number} {"G : {:.6f} Eh".format(self._last_energy["G"]) if self._last_energy["G"] else "E : {:.6f} Eh".format(self._last_energy["E"])}\n'
        for a, pos in zip(self.atoms, self.last_geometry):
            x, y, z = pos
            txt += f' {a}\t{x:14f}\t{y:14f}\t{z:14f}\n'
        return txt.strip()

    def create_log(self):
        en = self._last_energy
        number, e, g, b, erel, time, pop = self.number, en.get('E', float(0)), en.get('G', float(0)), en.get('B', float(0)), en.get('Erel', float(0)), en.get('time'), en.get('Pop', float(0))
        if g: g /= 627.51
        return number, e/627.51, g, b, erel, pop, time



    def __str__(self) -> str:
        return self.write_xyz()
    
    def __repr__(self) -> str:
        return self.write_xyz()
    


    # Functions needed for sorting the conformers' ensemble

    def __lt__(self, other):
        if not self.active: return 0 < other.get_energy
        return self.get_energy < other.get_energy
        
    def __gt__(self, other):
        if not self.active: return 0 > other.get_energy
        return self.get_energy > other.get_energy
        
    def __eq__(self, other):
        if not self.active: return 0 == other.get_energy
        return self.get_energy == other.get_energy