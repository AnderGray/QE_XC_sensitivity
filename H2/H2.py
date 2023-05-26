from ase import build
from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS

import os


### 
#   Compute total energy for H2
### 

os.environ['ASE_ESPRESSO_COMMAND']="/Users/akgray/Documents/DFT/q-e/bin/pw.x -in PREFIX.pwi > PREFIX.pwo"
os.environ['ESPRESSO_PSEUDO']="/Users/akgray/Documents/DFT/q-e/pseudo"

# ASE_ESPRESSO_COMMAND="/path/to/pw.x -in PREFIX.pwi > PREFIX.pwo"
pseudopotentials = {'H': 'H.pbe-kjpaw_psl.1.0.0.UPF'}
model=build.molecule('H2', vacuum=5)
# calc = Espresso(pseudopotentials=pseudopotentials,
#                 tstress=True, tprnfor=True, kpts=(1, 1, 1)) #, xc='BEEF-vdW')

calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(1, 1, 1), xc='BEEF-vdW')

model.calc = calc
opt = BFGS(model, trajectory='H2.traj')


opt.run(fmax=0.05)

H2_tot_energy = model.get_total_energy()

###
#   Compute disassosiation energy by using the theoretical ground state for 2H
#
#   Compare to experimental results
#
#       Experimental results for H2 energy from: 
#           https://www.sciencedirect.com/science/article/pii/0022285261901114? ref=pdf_download&fr=RR-2&rr=7cd53f6bef527315
#
#
#       Experimental result for bond length from NIST: https://cccbdb.nist.gov/exp2x.asp?casno=1333740#1979HUB/HER
###

H_ground_state = -13.6 # ev  ASK DUC if this is the right thing to do

cm_2_ev = 1 / 8065.54

H2_experimental = 36113.0 
H2_experimental_error = 0.3 

H2_experimental_ev = 4.47718
H2_experimental_ev_error = 0.00012

disassosiation_energy = 2 * H_ground_state - H2_tot_energy

H2_bond_length_NIST = 0.7414 # (Å)

H2_DFT_bond_length = model.get_distance(0,1)

print()
print("*********************************")
print("H2 disassosiation energy")
print(f"DTF: {disassosiation_energy}  (ev)")
print(f"EXP: {H2_experimental_ev} ± {H2_experimental_ev_error}  (ev)")

print("*********************************")
print(f"DFT H2 total: {H2_tot_energy}  (ev)")
print("*********************************")
print(f"DFT bond length: {H2_DFT_bond_length}  (Å)")
print(f"EXP bond length: {H2_bond_length_NIST}  (Å)")

# model_single=build.molecule('H', vacuum=5)
# model_single.calc = calc

# model_single.get_total_energy()

# opt_single = BFGS(model_single, trajectory='H.traj')
# opt_single.run(fmax=0.05)