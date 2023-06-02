from ase import build
from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS
import numpy as np

from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState

from ase.io.trajectory import Trajectory

import os


### 
#   Compute total energy for H2
### 

os.environ['ASE_ESPRESSO_COMMAND']="/Users/akgray/Documents/DFT/q-e/bin/pw.x -in PREFIX.pwi > PREFIX.pwo"
os.environ['ESPRESSO_PSEUDO']="/Users/akgray/Documents/DFT/q-e/pseudo"

# functionals = ['LDA', 'PBE', 'BEEF-vdW']

input_data = {
    'system': {
        #'input_dft': 'BEEF-vdW'
        #'ensemble_energies': '.true.'
        }
    }
    

# ASE_ESPRESSO_COMMAND="/path/to/pw.x -in PREFIX.pwi > PREFIX.pwo"
pseudopotentials = {'H': 'H.pbe-kjpaw_psl.1.0.0.UPF'}
model=build.molecule('H2', vacuum=5)
# calc = Espresso(pseudopotentials=pseudopotentials,
#                 tstress=True, tprnfor=True, kpts=(1, 1, 1)) #, xc='BEEF-vdW')

calc = Espresso(pseudopotentials=pseudopotentials,
                kpts=(1, 1, 1), input_data=input_data)

model.calc = calc

cell = model.get_cell()

volumes=[]
energies=[]

traj = Trajectory('H2.traj', 'w')
for x in np.linspace(0.95, 1.05, 10):
    model.set_cell(cell * x, scale_atoms=True)
    energies.append(model.get_potential_energy())
    volumes.append(model.get_volume())
    traj.write(model)


eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()

print(B / kJ * 1.0e24, 'GPa')
eos.plot('H2-eos.png')