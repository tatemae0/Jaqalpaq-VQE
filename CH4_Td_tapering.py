import numpy as np
import openfermion as of
import matplotlib.pyplot as plt
from openfermion.utils import commutator
from openfermion.ops import FermionOperator
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, eigs
from scipy.sparse import identity, find, csr_matrix, save_npz
from pygsti.tools import to_unitary

def save_sparse_csr(filename, array):
    # note that .npz extension is added automatically
    np.savez(filename, data=array.data, indices=array.indices,
            indptr=array.indptr, shape=array.shape)

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

print(f"Type of R is {type(R)}")
w, V = eigh(R) # Return eigenvalues, eigenvectors of permutation matrix. 

# R = R.todense()

# try:
#     save_sparse_csr("V_matrix_save_sparse_Csr", V)
# except:
#     pass

# Enforce unitary 
_ , V_u = to_unitary(V)
_, R_u = to_unitary(R)
S = V_u*R_u*np.conj(V_u) # S is a symmetry of the Hamiltonian.
# 
# For multiple symmetries, consider simultaneous diagonalizations (i.e., same V matrix?)
print(S.shape)
print(type(S))
try:
    S = csr_matrix(S)
    print(type(S))
except:
    pass
with open("CH4_Td_sym1_matrix.txt", "w") as symmetry_output:
    symmetry_output.writelines(str(S))
try:
    save_npz("CH4_sym1_save_npz", S)
except:
    pass
try:
    save_sparse_csr("CH4_sym1_save_sparse_csr", S)
except:
    pass

try:
    with open("CH4_Td_sym1_values.txt", "w") as symmetry_output:
        symmetry_output.writelines(find(S))
except: 
    pass

# Inspect S matrix to determine the right Z operator. This time it's manual, but going forward 
# we need a programmatic way to make a vector of the diagonal of S (use np.diagonal?), 
# extract the -1 eigenvalues, and 
# translate that into a pauli string.

# Then verify that S indeed commutes with pauli terms in Hamiltonian (using openfermion's 'commmutator')



# plt.imshow(dense_operator)
# plt.savefig("operator_visualization")