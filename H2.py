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
def ansatz(theta, sample_noise=False)


term_probs = []
for i in range(len(terms))
builder = circuitbuilder.CircuitBuilder(
native_gates=normalize_native_gates(native_gates.NATIVE_GATES))

# Create a qubit register
q = builder.register('q', 2)
# Define a hadamard macro
hadamard = circuitbuilder.SequentialBlockBuilder()
hadamard.gate('Sy', 'a')
hadamard.gate('Px', 'a')
builder.macro('hadamard', ['a'], hadamard)
# Prepare the Hartree Fock state
builder.gate('prepare_all')
builder.gate('Px', q[0])
# Apply the UCC Ansatz exp[-i*theta(X1 Y0)]
builder.gate('MS', q[1], q[0], 0, np.pi/2)
builder.gate('Rz', q[1], theta)
builder.gate('MS', q[1], q[0], 0, -np.pi/2)
# Change basis for measurement depending on term
for j, qubit in enumerate(terms[i]):
if qubit == 'X':
builder.gate('hadamard', ('array_item', q, j)),
if qubit == 'Y':
builder.gate('Sxd', ('array_item', q, j)),
builder.gate('measure_all')
circuit = builder.build()
# Format results of simulation as a list of lists
sim_result = emulator.run_jaqal_circuit(circuit)
sim_probs = sim_result.subcircuits[0].probability_by_int
if sample_noise:  # Sample circuits to determine probs
probs = np.zeros(4)  # Number of possible states
for k in range(n_samples):
sample = np.random.choice(4, p=sim_probs)
probs[sample] += 1  # Increment state counter
probs = probs/n_samples  # Determine probabilities from sampling
term_probs += [probs]  # Combine lists of probs of each term in Hamiltonian
else:  # Exact solution without sampling
term_probs += [sim_probs]
return term_probs
# Calculate energy of one term of the Hamiltonian for one possible state
def term_energy(term, state, coefficient, prob):


parity = 1
for i in range(len(term)):
# Change parity if state is occupied and is acted on by a pauli operator
if term[i] != None and state[i] == '1':
parity = -1*parity
return coefficient*prob*parity
# Calculate energy of the molecule for a given value of theta
def make_calculate_energy(sample_noise=False):
def calculate_energy(theta):


energy = 0
# Convert tuple (from optimization) to float
probs = ansatz(theta[0], sample_noise)
for i in range(len(terms)):  # For each term in the hamiltonian

for j in range(len(probs[0])):  # For each possible state
term = terms[i]
state = '{0:02b}'.format(j)[::-1]  # convert state to binary (# of qubits)
# binary must be inverted due to jaqalpaq convention
coefficient = cs[i].real
prob = probs[i][j]
energy += term_energy(term, state, coefficient, prob)
return energy
return calculate_energy
# Set the basis set, spin, and charge of the H2 molecule
basis = 'sto-3g'
multiplicity = 1
charge = 0
# Set calculation parameters
run_scf = 1
run_fci = 1
delete_input = True
delete_output = False
optimized_energies = [[], []]
exact_energies = []
# Loop over bond lengths from 0.3 to 1.3 angstroms
sample_noise = False
n_samples = 10000  # Sample circuit
n_pts = 11  # Number of points
bond_lengths = np.linspace(0.3, 1.3, n_pts)
for diatomic_bond_length in bond_lengths:
# Generate molecule at some bond length
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., diatomic_bond_length))]
molecule = MolecularData(
geometry, basis, multiplicity, charge,
description=str(round(diatomic_bond_length, 2)),
filename='./H2_sto-3g_single_dissociation')
# Run pyscf to generate new molecular data for sto-3g H2
molecule = run_pyscf(molecule,
run_scf=run_scf,
run_fci=run_fci,
verbose=False)
# Get the fermionic Hamiltonian for H2 and map it into qubits using the BK encoding
hamiltonian = molecule.get_molecular_hamiltonian()
hamiltonian_ferm = of.get_fermion_operator(hamiltonian)
hamiltonian_bk = of.bravyi_kitaev(hamiltonian_ferm)
# Define Pauli strings that appear in the reduced two-qubit Hamiltonian
# q0 , q1
terms = [[None, None], ['Z', None], [None, 'Z'],
    ['Z', 'Z'], ['X', 'X'], ['Y', 'Y']]
# Calculate effective coefficients for the reduced two-qubit Hamiltonian
# Derivation follows arXiv:1803.10238v2 appendix A-2
fs = hamiltonian_bk.terms  # Old coefficients from OpenFermion Hamiltonian
c0 = (fs[()] + fs[(1, 'Z'), ] + fs[(1, 'Z'), (3, 'Z'), ]).real
c1 = (fs[(0, 'Z'), ] + fs[(0, 'Z'), (1, 'Z'), ]).real
c2 = (fs[(2, 'Z'), ] + fs[(1, 'Z'), (2, 'Z'), (3, 'Z'), ]).real
c3 = (fs[(0, 'Z'), (2, 'Z'), ] + fs[(0, 'Z'), (1, 'Z'), (2, 'Z'), ] + fs[(0, 'Z'), (2, 'Z'),

for j in range(len(probs[0])):  # For each possible state
term=terms[i]
state='{0:02b}'.format(j)[::-1]  # convert state to binary (# of qubits)
# binary must be inverted due to jaqalpaq convention
coefficient=cs[i].real
prob=probs[i][j]
energy += term_energy(term, state, coefficient, prob)
return energy
return calculate_energy
# Set the basis set, spin, and charge of the H2 molecule
basis='sto-3g'
multiplicity=1
charge=0
# Set calculation parameters
run_scf=1
run_fci=1
delete_input=True
delete_output=False
optimized_energies=[[], []]
exact_energies=[]
# Loop over bond lengths from 0.3 to 1.3 angstroms
sample_noise=False
n_samples=10000  # Sample circuit
n_pts=11  # Number of points
bond_lengths=np.linspace(0.3, 1.3, n_pts)
for diatomic_bond_length in bond_lengths:
# Generate molecule at some bond length
geometry=[('H', (0., 0., 0.)), ('H', (0., 0., diatomic_bond_length))]
molecule=MolecularData(
geometry, basis, multiplicity, charge,
description=str(round(diatomic_bond_length, 2)),
filename='./H2_sto-3g_single_dissociation')
# Run pyscf to generate new molecular data for sto-3g H2
molecule=run_pyscf(molecule,
run_scf=run_scf,
run_fci=run_fci,
verbose=False)
# Get the fermionic Hamiltonian for H2 and map it into qubits using the BK encoding
hamiltonian=molecule.get_molecular_hamiltonian()
hamiltonian_ferm=of.get_fermion_operator(hamiltonian)
hamiltonian_bk=of.bravyi_kitaev(hamiltonian_ferm)
# Define Pauli strings that appear in the reduced two-qubit Hamiltonian
# q0 , q1
terms=[[None, None], ['Z', None], [None, 'Z'],
    ['Z', 'Z'], ['X', 'X'], ['Y', 'Y']]
# Calculate effective coefficients for the reduced two-qubit Hamiltonian
# Derivation follows arXiv:1803.10238v2 appendix A-2
fs=hamiltonian_bk.terms  # Old coefficients from OpenFermion Hamiltonian
c0=(fs[()] + fs[(1, 'Z'), ] + fs[(1, 'Z'), (3, 'Z'), ]).real
c1=(fs[(0, 'Z'), ] + fs[(0, 'Z'), (1, 'Z'), ]).real
c2=(fs[(2, 'Z'), ] + fs[(1, 'Z'), (2, 'Z'), (3, 'Z'), ]).real
c3=(fs[(0, 'Z'), (2, 'Z'), ] + fs[(0, 'Z'), (1, 'Z'), (2, 'Z'), ] + fs[(0, 'Z'), (2, 'Z')
