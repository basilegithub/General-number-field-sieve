# This is the file to build the binary matrix

from utils import *
from polynomial_functions import *

# Create the sparse matrix:
# Consider the dense binary matrix A[i][j]
# The sparse matrix M is defined by M[i] contains j if and only if A[i][j] = 1
# For large numbers to factor, this reduces that memory requirement for the linear algebra part
# Each prime in the factor base corresponds to one line
# Each relation corresponds to one column
def build_sparse_matrix(pairs, rat, alg, qua, div_lead):
    matrix = []
    for r in pairs:
        index = 1
        
        if r[3] < 0: line = [0]
        else: line = []

        for p in rat:
            tmp = p
            tmp2 = 0
            while not r[3]%tmp:
                tmp *= p
                tmp2 ^= 1
            if tmp2: line.append(index)
            index += 1
            
        if r[1] < 0: line.append(index)
        index += 1

        for p in range(len(alg)):
            tmp = rat[p]
            tmp2 = 0
            while not r[1]%tmp:
                tmp *= rat[p]
                tmp2 ^= 1

            for i in range(len(alg[p])):
                if not eval_mod(r[5], alg[p][i], rat[p]):
                    if tmp2: line.append(index)
                index += 1
                
        for p in range(len(qua)):
            if compute_legendre_character(eval_mod(r[5], qua[p][1], qua[p][0]), qua[p][0]) == -1: line.append(index)
            index += 1
            
        for p in range(len(div_lead)):
            if r[7+p]:
                tmp = div_lead[p]
                tmp2 = 0

                while not r[1]%tmp:
                    tmp *= div_lead[p]
                    tmp2 ^= 1
                if tmp2: line.append(index)

            index += 1
            
        if r[6]&1: line.append(index)
        matrix.append(set(line))

    return matrix
    
# Reduce the sparse matrix, according to various rules:
# 1. If a line contains only one non-zero value, then it is deleted as well as the column containing the corresponding non-zero value
# 2. If a line contrains no non-zero value, then it is deleted
# 3. We only keep 10 more columns than lines, to ensure we still have solutions while reducing the matrix size
def reduce_sparse_matrix(matrix, pairs):
    flag = True
    while flag:
        flag = False

        singleton_queue = [index for index, row in enumerate(matrix) if len(row) == 1]
        active_cols = set(range(len(pairs)))

        flag = len(singleton_queue)

        while singleton_queue:
            index = singleton_queue.pop()

            if len(matrix[index]) != 1: continue
            coeff = next(iter(matrix[index]))

            for j, row in enumerate(matrix):
                if j != index and coeff in row:
                        row.remove(coeff)
                        if len(row) == 1: singleton_queue.append(j)

            matrix[index].clear()
            active_cols.discard(coeff)

        matrix = [row for row in matrix if row]

        pairs = [pairs[index] for index in sorted(active_cols)]

        mapping = {old: new for new, old in enumerate(sorted(active_cols))}
        matrix = [{mapping[c] for c in row if c in mapping} for row in matrix]

        active_cols = set(range(len(pairs)))
        target_n_cols = len(matrix)+10
        to_delete = len(pairs) - target_n_cols
        current_len_row = 2

        while to_delete >= current_len_row-1:

            index = -1
            for i in range(len(matrix)):
                if len(matrix[i]) == current_len_row:
                    flag = True
                    index = i
                    break

            if index != -1:
                for coeff in matrix[index]:
                    active_cols.discard(coeff)

                    for j, row in enumerate(matrix):
                        if j != index and coeff in row:
                            row.discard(coeff)

                to_delete -= current_len_row - 1

                matrix[index].clear()

            else:
                current_len_row += 1

        matrix = [row for row in matrix if row]

        pairs = [pairs[index] for index in sorted(active_cols)]

        mapping = {old: new for new, old in enumerate(sorted(active_cols))}
        matrix = [{mapping[c] for c in row if c in mapping} for row in matrix]

    return matrix, pairs
                        
# Build the dense binary matrix by computing explicitely every A[i][j]
def build_dense_matrix(pairs, rat, alg, qua, div_lead):
    matrix = []
    for r in pairs:
        line = []
        
        if r[3] < 0: line.append(1)
        else: line.append(0)

        for p in rat:
            tmp = p
            tmp2 = 0
            while not r[3]%tmp:
                tmp *= p
                tmp2 ^= 1
            line.append(tmp2)
            
        if r[1] < 0: line.append(1)
        else: line.append(0)

        for p in range(len(alg)):
            tmp = rat[p]
            tmp2 = 0
            while not r[1]%tmp:
                tmp *= rat[p]
                tmp2 ^= 1
            for i in range(len(alg[p])):
                if not eval_mod(r[5], alg[p][i], rat[p]):
                    line.append(tmp2)
                else:
                    line.append(0)
                
        for p in range(len(qua)):
            if compute_legendre_character(eval_mod(r[5], qua[p][1], qua[p][0]), qua[p][0]) == -1:
                line.append(1)
            else:
                line.append(0)
            
        for p in range(len(div_lead)):
            if r[7+p]:
                tmp = div_lead[p]
                tmp2 = 0
                while not r[1]%tmp:
                    tmp *= div_lead[p]
                    tmp2 ^= 1
                line.append(tmp2)
            else: line.append(0)
            
        line.append(r[6]&1)
        matrix.append(line)
        
    return matrix
    
def siqs_build_matrix_opt(M):
    """Convert the given matrix M of 0s and 1s into a list of numbers m
    that correspond to the columns of the matrix.
    The j-th number encodes the j-th column of matrix M in binary:
    The i-th bit of m[i] is equal to M[i][j].
    """
    m = len(M[0])
    cols_binary = [""] * m
    for mi in M:
        for j, mij in enumerate(mi):
            cols_binary[j] += "1" if mij else "0"
    return [int(cols_bin[::-1], 2) for cols_bin in cols_binary], len(M), m