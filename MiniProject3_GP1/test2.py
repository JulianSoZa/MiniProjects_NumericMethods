import numpy as np

from scipy.sparse import csr_matrix
dataM = [1, 1, 1, 1]
row = [1, 1, 0, 0]
col = [1, 1, 0, 0]
M = csr_matrix((dataM, (row, col)))

print(M.diagonal())   