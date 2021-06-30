import openfermion as of

from openfermion.linalg import get_sparse_operator  
from openfermion.utils import commutator
from openfermion.ops import FermionOperator
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
from pygsti.tools import to_unitary

# Second quantized representation of H2 reflection (perpendicular to bond axis)
# t1 = FermionOperator('1 1^', -1)
# t2 = FermionOperator('2 2^', -1)
# t3 = FermionOperator('1 2^')
# t4 = FermionOperator('1^ 2') 

# Second quantized representation of CH4 fermionic mode swap (reflection sigma_d) 
t1 = FermionOperator('2^ 2', -1)
t2 = FermionOperator('12^ 12', -1)
t3 = FermionOperator('2^ 12')
t4 = FermionOperator('2 12^') 


full_fermionic_operator = t1 + t2 + t3 + t4
sparse_operator = get_sparse_operator(full_fermionic_operator).astype(int) # Cast as sparse operator
# dense_operator = get_sparse_operator(full_fermionic_operator).todense() # Optional: Dense rep.
R = np.identity(sparse_operator.shape[0]) + sparse_operator # Add identity element to yield a proper permutation matrix.

w, V = eigh(R) # Return eigenvalues, eigenvectors of permutation matrix. 

# Enforce unitary 
_ , V_u = to_unitary(V)
_, R_u = to_unitary(R)
S = V_u*R_u*np.conj(V_u) # S is a symmetry of the Hamiltonian.
# For multiple symmetries, consider simultaneous diagonalizations (i.e., same V matrix?)
print(S)
# Find how to rewrite this as a product of pauli Z operators.
# Then verify that S indeed commutes with pauli terms in Hamiltonian (using openfermion's 'commmutator')



# plt.imshow(dense_operator)
# plt.savefig("operator_visualization")