#%%%
import ase
from ase.atoms import Atoms
from ase.calculators.orca import ORCA
from ase.optimize import LBFGS
from ase.vibrations import Vibrations
import os
import scipy as sp

calc = ORCA(label='ORCA',
              maxiter=2000,
              charge=0, mult=1,task='gradient',
              orcasimpleinput='PBE def2-SVP opt freq',
              orcablocks=''
              )


os.getcwd()

#%%
water = ase.io.read('struct/water.xyz')# calculator=calc)
water.calc = calc

opt = LBFGS(atoms=water, trajectory='struct/waret.trj', logfile='struct/water.log')
opt.run()


#%%
print(water.get_potential_energy())
print(opt.__dict__.keys())


# %%
import cclib

data = cclib.io.ccread('ORCA.out')
# %%
data.__dict__
print(data.writexyz())
# %%


os.system('rm ORCA*')
# %%
