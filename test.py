import numpy as np
import matplotlib.pyplot as plt
from openfermion.linalg import get_sparse_operator
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse import identity, find, load_npz

# A = csr_matrix([
# [7,0,0],
# [0,4,5],
# [0,0,0]
# ])

# print(str(A))

matrix = load_npz("CH4_sym1_save_npz.npz").todense()
# matrix = numpy.load("CH4_sym1_save_npz.npz")
print(matrix.shape)
# print(matrix.diagonal().value)
print(f"------- Count nonzero method------- \n{np.count_nonzero(matrix)}")
print(f"-----bincount ---------\n {np.unique(matrix, return_counts=True)}") 
# plt.imshow(matrix, cmap="plasma")
# plt.savefig("matrix_load_TEST")