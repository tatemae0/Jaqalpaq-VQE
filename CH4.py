# Imports for QSCOUT
import jaqalpaq
from jaqalpaq.core import circuitbuilder
from jaqalpaq.core.circuit import normalize_native_gates
from jaqalpaq import emulator
from qscout.v1 import native_gates

# Imports for basic mathematical functionality
from math import pi
import numpy as np
from scipy.sparse import load_npz
from scipy.sparse import csr_matrix

# Imports for OpenFermion(-PySCF)
import openfermion as of
from openfermion.linalg import get_sparse_operator
from openfermion.chem import MolecularData
from openfermion.chem import geometry_from_pubchem
from openfermionpyscf import run_pyscf
from openfermion.ops import QubitOperator 
from openfermion.utils import commutator
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
methane_geometry = geometry_from_pubchem('methane')


molecule = MolecularData(
    methane_geometry, basis, multiplicity, charge,    
    filename='./CH4_sto-6g_single')

molecule = run_pyscf(molecule,
                     run_scf=run_scf,
                     run_fci=run_fci,
                     verbose=False)

hamiltonian = molecule.get_molecular_hamiltonian()
hamiltonian_ferm = of.get_fermion_operator(hamiltonian)
hamiltonian_bk = of.bravyi_kitaev(hamiltonian_ferm, n_qubits=18).terms
ham_bk_list = list(hamiltonian_bk.keys())
ham_bk_sparse_terms = []
commutators = []

# Test case below, using only the known sym1. Generalize to any symmetry. This is already
# like 8900+ evals (one sym). Will have to do pairwise combinations for multiple symmetries.
# To extend, replace all instances of sym1 with a for loop variable, and make a new loop
# outside.
sym1 = load_npz("CH4_sym1_save_npz.npz")
print(f"symmetry 1 is of type {type(sym1)}")
print(sym1.shape)

symmetries = [sym1]
for sym in symmetries:
    term_counter = 0
    for term in ham_bk_list:
        if term_counter <= 10 or (100 <= term_counter < 120):
            sparse_term = csr_matrix(get_sparse_operator(QubitOperator(term)))
            # ham_bk_sparse_terms.append(sparse_term)
            print(sparse_term.shape)
            # commutators.append(commutator(sparse_term, sym)) # useful for multiple iterations.
            # # type cast to scipy.sparse for ease of comparison
            # with open("sym_CH4_bk_commutator.txt", "w") as commutator_output:
            #     commutator_output.writelines(commutators)
        term_counter += 1
    for i in range(10): # Print out the smallest 10 commutations.
        print(commutators.min())
        pop_index = commutators.min.index()
        commutators.pop(pop_index)
    # TODO: resolve dimension conflicts of operators in commutation calculation. 
    # sparse_hamiltonian_bk = get_sparse_operator(hamiltonian_bk[0])


# ADDENDUM (delete later): Does the example "symmetry" actually commute with any term in H?
# This BK hamiltonian isn't hermitian... Some coefficients are imaginary (albeit e-8 order).

# with open("CH4_base_Hamiltonian.txt", "w") as hamiltonian_output:
#     hamiltonian_output.writelines(str(hamiltonian))
# with open("CH4_fermion_Hamiltonian.txt", "w") as hamiltonian_output:
#     hamiltonian_output.writelines(str(hamiltonian_ferm))
# with open("CH4_bk_Hamiltonian.txt", "w") as hamiltonian_output:
#     hamiltonian_output.writelines(str(hamiltonian_bk))
# for orbital in range(molecule.n_orbitals):
#     print('Spatial orbital {} has energy of {} Hartree.'.format(
#         orbital, molecule.orbital_energies[orbital]))