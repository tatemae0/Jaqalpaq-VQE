# Imports for QSCOUT
import jaqalpaq
from jaqalpaq.core import circuitbuilder
from jaqalpaq.core.circuit import normalize_native_gates
from jaqalpaq import emulator
from qscout.v1 import native_gates
# Imports for basic mathematical functionality
from math import pi
import numpy as np
# Imports for OpenFermion(-PySCF)
import openfermion as of
from openfermion.chem import MolecularData
from openfermionpyscf import run_pyscf
# Import for VQE optimizer
from scipy import optimize

# Build molecule

# Set the basis set, spin, and charge of CH4
basis = 'sto-6g'
multiplicity = 1
charge = 0

# Set calculation parameters
run_scf = 1
run_fci = 1
delete_input = True
delete_output = False
from openfermion.chem import geometry_from_pubchem
methane_geometry = geometry_from_pubchem('methane')


molecule = MolecularData(
    methane_geometry, basis, multiplicity, charge,    
    filename='./CH4_sto-6g_single')

molecule = run_pyscf(molecule,
                     run_scf=run_scf,
                     run_fci=run_fci,
                     verbose=False)

hamiltonian = molecule.get_molecular_hamiltonian()
with open("CH4_base_Hamiltonian.txt", "w") as hamiltonian_output:
    hamiltonian_output.writelines(str(hamiltonian))
hamiltonian_ferm = of.get_fermion_operator(hamiltonian)
with open("CH4_fermion_Hamiltonian.txt", "w") as hamiltonian_output:
    hamiltonian_output.writelines(str(hamiltonian_ferm))
hamiltonian_bk = of.bravyi_kitaev(hamiltonian_ferm)
with open("CH4_bk_Hamiltonian.txt", "w") as hamiltonian_output:
    hamiltonian_output.writelines(str(hamiltonian_bk))
for orbital in range(molecule.n_orbitals):
    print('Spatial orbital {} has energy of {} Hartree.'.format(
        orbital, molecule.orbital_energies[orbital]))