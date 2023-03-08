import numpy as np
from ase.atoms import Atoms


class Conformer:

    def __init__(self, number: int, geom: np.array, atoms: np.array, charge : int= 0, mult : int = 1) -> None:
        self.number = number
        self._initial_geometry = geom
        self.last_geometry = geom
        self.atoms = atoms
        self.energies = {}
        self.active = True
    
    def get_ase_atoms(self, calc):
        return Atoms(
            symbols = ''.join(list(self.atoms)),
            positions = self.last_geometry,
            calculator=calc
        )
    
    @property
    def rotatory(self):
        return self.energies[list(self.energies.keys())[-1]]['B']

    @property
    def get_energy(self):
        en = self.energies[list(self.energies.keys())[-1]]
        if en['G']: return en['G']
        return en['E']
    
    def write_xyz(self):
        txt = f'{len(self.atoms)}\nCONFORMER {self.number}\n'
        for a, pos in zip(self.atoms, self.last_geometry):
            x, y, z = pos
            txt += f' {a}\t{x:14f}\t{y:14f}\t{z:14f}\n'
        return txt.strip()

    def __str__(self) -> str:
        return self.write_xyz()
    
    def __repr__(self) -> str:
        return self.write_xyz()
    


    # Functions needed for sorting the conformers' ensemble

    def __lt__(self, other):
        return self.get_energy < other.get_energy
        
    def __gt__(self, other):
        return self.get_energy > other.get_energy
        
    def __eq__(self, other):
        return self.get_energy == other.get_energy