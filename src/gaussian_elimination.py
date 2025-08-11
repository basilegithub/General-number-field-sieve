# This file contains the functions that perform gaussian elimination on a dense matrix

from utils import *

def add_column_opt(M_opt, tgt, src):
    """For a matrix produced by siqs_build_matrix_opt, add the column
    src to the column target (mod 2).
    """
    M_opt[tgt] ^= M_opt[src]


def find_pivot_column_opt(M_opt, j):
    """For a matrix produced by siqs_build_matrix_opt, return the row of
    the first non-zero entry in column j, or None if no such row exists.
    """
    if M_opt[j] == 0:
        return None
    return lowest_set_bit(M_opt[j])


def siqs_solve_matrix_opt(M_opt, n, m):
    """
    Perform the linear algebra step of the SIQS. Perform fast
    Gaussian elimination to determine pairs of perfect squares mod n.
    Use the optimisations described in [1].

    [1] Koç, Çetin K., and Sarath N. Arachchige. 'A Fast Algorithm for
        Gaussian Elimination over GF (2) and its Implementation on the
        GAPP.' Journal of Parallel and Distributed Computing 13.1
        (1991): 118-122.
    """
    row_is_marked = [False] * n
    pivots = [-1] * m
    for j in range(m):
        i = find_pivot_column_opt(M_opt, j)
        if i is not None:
            pivots[j] = i
            row_is_marked[i] = True
            for k in range(m):
                if k != j and (M_opt[k] >> i) & 1:  # test M[i][k] == 1
                    add_column_opt(M_opt, k, j)
    perf_squares = []
    for i in range(n):
        if not row_is_marked[i]:
            perfect_sq_indices = [i]
            for j in range(m):
                if (M_opt[j] >> i) & 1:  # test M[i][j] == 1
                    perfect_sq_indices.append(pivots[j])
            perf_squares.append(perfect_sq_indices)
    return perf_squares