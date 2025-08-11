# This is the file containing the functions that perform the block Wiedemann algorithm

import random
from utils import sparse_multiply, sparse_dot_prod, poly_div_bin, poly_prod_bin, poly_gcd_bin

# from a block of vectors, generate the associated block of numbers whose bits encode the vectors entries
def create_block(vectors,n,block_size):
    block = [0]*n
    for i in range(n):
        line = 0
        for j in range(block_size-1):
            line ^= vectors[j][i]
            line <<= 1
        line ^= vectors[-1][i]
        block[i] = line
    return block
    
# from a block of numbers, generate the associated block of vectors 
def block_to_vec(block,block_size,n):
    res = [0]*block_size
    for i in range(block_size):
        vec = []
        for j in range(n):
            vec.append(block[j] >> block_size-1-i & 1)
        res[i] = vec
    return res

# compute the minimal polynomial of the sequence using Berlekamp-Massey algorithm
def poly_anul(sequence,m):
    A = 1 << (m<<1)
    B = sequence
    C = 0
    D = 1
    while B.bit_length() > m:
        (Q,R) = poly_div_bin(A,B)
        E = C^poly_prod_bin(Q,D)
        C, D, A, B = D, E, B, R
    return D

# For each vector v in the block, and with the current estimate of the minimal polynomial P
# Finds the largest m such that P(A) = A^m * Q(A) for some polynomial Q
# Compute Q(A)v and then finds the highest power n such that A^n Q(A)v != 0 and A^(n+1) Q(A)v == 0
# Return A^n Q(A)v for each v (n depends on v)
def compute_kernel_vectors(block,A,poly,n,block_size,d):
    cpt = 0
    while not poly&1:
        cpt += 1
        poly >>= 1
    
    res = [0]*n
    for i in range(poly.bit_length()-1, 0, -1):
        tmp = (poly >> i)&1
        for j in range(n): res[j] ^= block[j]*tmp
        res = sparse_multiply(A,res)
    tmp = poly&1
    for j in range(n): res[j] ^= block[j]*tmp
    
    selected_columns = []
    null_space = [0]*n
    
    for _ in range(cpt):
        tmp = sparse_multiply(A, res)
        
        for i in range(block_size):
            if i not in selected_columns and all(not((tmp[j] >> i) & 1) for j in range(n)):
                selected_columns.append(i)
                for j in range(n):
                    null_space[j] ^= (res[j] >> i)&1
                    null_space[j] <<= 1
                    
        res = list(tmp)
        
    tmp = sparse_multiply(A, res)
    
    for i in range(block_size):
        if i not in selected_columns and all(not((tmp[j] >> i) & 1) for j in range(n)):
            selected_columns.append(i)
            for j in range(n):
                null_space[j] ^= (res[j] >> i)&1
                null_space[j] <<= 1
                    
    for j in range(n):
        null_space[j] >>= 1
    
    return null_space
    
# Implementation of the Wiedemann algorithm
# Not been able yet to implement the block version of the algorithm
# Here I only use the binary encoding to compute multiple scalar linear sequences at once
def wiedemann(A, n, block_size, mini_poly_estim):
    vectors = []
    for i in range(block_size):
        vec = [random.randint(0,1) for j in range(n)]
        vectors.append(vec)
    block = create_block(vectors,n,block_size)
    
    stored = [i for i in block] # stores the solution to AX = Y where Y is the block of vectors
    if mini_poly_estim.bit_length()-1 < len(A):
        block = sparse_multiply(A,block)
        # Sample random left projection
        tmp = [random.randint(0,1) for i in range(n)]
        lbd = []
        for i in range(n):
            if tmp[i]: lbd.append(i) # sparse encoding of the left projection tmp
            
        # Compute the sequence (lbd.A^i.block) for growing i
        # We will find the generator for each scalar sequence, and compute the lcm of generator and current estimate of minimal polynomial
        tmp2 = block
        sequence = [0]*block_size
        for i in range(2*n-1):
            tmp = sparse_dot_prod(lbd,tmp2)
            for j in range(block_size):
                if i < n<<1:
                    sequence[j] ^= (tmp >> block_size-j-1) & 1
                    if i+1 < n<<1: sequence[j] <<= 1
            tmp2 = sparse_multiply(A,tmp2)
        tmp = sparse_dot_prod(lbd,tmp2)
        for j in range(block_size):
            sequence[j] ^= (tmp >> block_size-j-1 & 1)
            
        for i in range(block_size):
            current_sequence = sequence[i]
            
            tmp_pi = poly_anul(current_sequence,n)
            deg = tmp_pi.bit_length()-1
            if deg > 0: # if we have found a non-trivial polynomial, update minimal polynomial estimate
                mini_poly_estim = poly_div_bin(poly_prod_bin(mini_poly_estim, tmp_pi), poly_gcd_bin(mini_poly_estim, tmp_pi))[0]
            
    stored = compute_kernel_vectors(stored,A,mini_poly_estim,n,block_size,mini_poly_estim.bit_length()-1)
    return block_to_vec(stored,block_size,n), mini_poly_estim
    
# Deletes duplicates in the list of null space vectors
def reduce_null_space_vectors(null_space):
    i = 0
    while i < len(null_space)-1:
        if not sum(null_space[i]): del null_space[i]
        else:
            j = i+1
            while j < len(null_space):
                if all(null_space[i][k] == null_space[j][k] for k in range(len(null_space[0]))) or sum(null_space[j]) == 0:
                    del null_space[j]
                else: j += 1
            i += 1
        
    return null_space