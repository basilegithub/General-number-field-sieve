# This file contains the required functions to perform block Lanczos algorithm

from utils import add_vector, dense_multiply, sparse_multiply, transpose_dense, identity, concatenate, transpose_sparse
import random
import log

# Efficiently performs d[i] ^= 1 for all i
def switch_indices(d):
    return ~d
    
# d selects the coefficients of A
def multiply_d(A, d):
    return [A[i]&d for i in range(len(A))]
    
# Compute W_inv and the indices for d
# Basically performs gaussian elimination
def block(T, N):     
    M = concatenate(T, identity(N), N)
    S = []
    
    for j in range(N):
        for k in range(j, N):
            if (M[k] >> 2*N - j - 1)&1 != 0:
                M[k], M[j] = M[j], M[k]
                break
        
        if (M[j] >> 2*N - j -1)&1 != 0:
            for k in range(N):
                if k != j:
                    if (M[k] >> 2*N - j - 1)&1 != 0:
                        M[k] ^= M[j]
            S.append(j)
        else:
            for k in range(j, N):
                if (M[k] >> N - j - 1)&1 != 0:
                    M[k] ,M[j] = M[j], M[k]
                    break
                
            if (M[j] >> N - j - 1)&1 == 0:
                for e in M:
                    print(e)
                return False
            for k in range(N):
                if k != j:
                    if (M[k] >> N - j - 1)&1 == 1:
                        M[k] ^= M[j]
                M[j] = 0
                
    new = [M[z]&((1<<N)-1) for z in range(N)]

    d = [0]*N
    for index in S: d[index] = 1
    d_new = 0
    for i in range(len(d)-1):
        d_new ^= d[i]
        d_new <<= 1
        
    d_new ^= d[-1]
    
    return new, d_new
    
# The whole block lanczos algorithm
# See the README for sources
def block_lanczos(B, nb_relations, N, LOG_PATH):
    Y = [random.randint(0, (1<<N)-1) for _ in range(nb_relations)]

    X = [0]*nb_relations
    b = transpose_sparse(B, nb_relations)
    Vo = sparse_multiply(b, sparse_multiply(B, Y))
    i = 0
    
    P = [0 for _ in range(nb_relations)]
    V = Vo
    d = 1
    while d and i <= int(len(B)/(N-0.764))+10:
        Z = sparse_multiply(b, sparse_multiply(B, V))
        vAv = dense_multiply(transpose_dense(V, N), Z)
        vAAv = dense_multiply(transpose_dense(Z, N), Z)
        
        W_inv, d = block(vAv,N)
        
        X = add_vector(X, dense_multiply(V, dense_multiply(W_inv, dense_multiply(transpose_dense(V, N), Vo))))
        
        neg_d = switch_indices(d)
        
        c = dense_multiply(W_inv, add_vector(multiply_d(vAAv, d), multiply_d(vAv, neg_d)))
        
        tmp1 = multiply_d(Z, d)
        tmp2 = multiply_d(V, neg_d)
        tmp3 = dense_multiply(V, c)
        tmp4 = multiply_d(vAv, d)
        tmp5 = dense_multiply(P, tmp4)
        tmp6 = dense_multiply(V, W_inv)
        tmp7 = multiply_d(P, neg_d)
            
        V, P = add_vector(add_vector(add_vector(tmp1, tmp2), tmp3), tmp5), add_vector(tmp6, tmp7)
        
        i += 1
 
    log.write_log(LOG_PATH, "lanczos halted after "+str(i)+" iterations")
    x = add_vector(X, Y)
    Z = concatenate(x, V, N)
    matrix = transpose_dense(sparse_multiply(B, Z), N<<1)
    Z = transpose_dense(Z, N<<1)
    matrix, Z = solve(matrix, Z, len(B))

    solutions = []
    for i in range(len(matrix)):
        if matrix[i] == 0 and Z[i] != 0 and Z[i] not in solutions:
            solutions.append(Z[i])

    if len(solutions) == 0:
        solutions = block_lanczos(B,nb_relations,N<<1,LOG_PATH)

    return solutions
    
# Performs gaussian elimination
def solve(matrix, block, nb_relations):
    k = 0
    
    for l in range(nb_relations):
        for i in range(k, len(matrix)):
            if (matrix[i] >> nb_relations - i - 1)&1 != 0:
                matrix[k], matrix[i] = matrix[i], matrix[k]
                block[k], block[i] = block[i], block[k]
                k += 1
                break
                
        for z in range(i+1, len(matrix)):
            if (matrix[z] >> nb_relations - l - 1)&1 == 1:
                matrix[z] ^= matrix[k-1]
                block[z] ^= block[k-1]

    return matrix, block