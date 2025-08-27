# This file contains the functions required for the sieving step

from polynomial_functions import *
import math

def sieve(length, f_x, rational_base, algebraic_base, m0, m1, b, logs, offset, leading,div):
    offset += div.bit_length()-1
    pairs = []
    tmp_poly = new_coeffs(f_x, b)

    sieve_array = [0]*(length<<1)
    for q, p in enumerate(rational_base):
        tmp_len = length%p
        log = logs[q]

        if m1%p:
            root = (tmp_len+b*m0*invmod(m1, p))%p
            for i in range(root, len(sieve_array), p): sieve_array[i] += log

        for r in algebraic_base[q]:
            root = (tmp_len+b*r)%p
            for i in range(root, len(sieve_array), p): sieve_array[i] += log

    if b&1:
        eval2 = -length*m1-b*m0
        a = -length
        tmp = [0]*len(tmp_poly)
        for j in range(len(tmp_poly)): tmp[j] = evaluate(tmp_poly, -length+j)
        for q in range(1, len(tmp_poly)):
            for k in range(len(tmp_poly)-1, q-1, -1): tmp[k] -= tmp[k-1]

        eval1 = tmp[0]
        for k in range(len(sieve_array)):
            if a and math.gcd(a, b) == 1 and eval2:
                eval = abs(eval1*eval2)
                if eval != 0 and sieve_array[k] > eval.bit_length()-offset:
                    pairs.append([[-b, a*leading], eval1, [[1 ,1]], eval2, [1], [-b,a], 1])
            a += 1
            eval2 += m1
            eval1 += tmp[1]
            for q in range(1, len(tmp_poly)-1): tmp[q] += tmp[q+1]
    else:
        init = 0
        if not length&1:
            length -= 1
            init = 1
        eval2 = -length*m1-b*m0
        a = -length
        tmp = [0]*len(tmp_poly)
        for j in range(len(tmp_poly)): tmp[j] = evaluate(tmp_poly, -length+2*j)
        for q in range(1, len(tmp_poly)):
            for k in range(len(tmp_poly)-1, q-1, -1): tmp[k] -= tmp[k-1]

        eval1 = tmp[0]
        for k in range(init, len(sieve_array), 2):
            if math.gcd(a, b) == 1 and eval2:
                eval = abs(eval1*eval2)
                if eval != 0 and sieve_array[k] > eval.bit_length()-offset:
                    pairs.append([[-b, a*leading], eval1, [[1, 1]], eval2, [1], [-b,a], 1])
            a += 2
            eval2 += m1<<1
            eval1 += tmp[1]
            for q in range(1, len(tmp_poly)-1): tmp[q] += tmp[q+1]

    return pairs