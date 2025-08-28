# This is the file conatining functions to generate the solutions

import math, sys
from square_roots import *
from polynomial_functions import *
from utils import my_norm
from datetime import datetime
import log

# Used when the gaussian pivot is used, convert sparse vector of indices into dense binary vector
def convert_to_binary(z, n):
    res = [0]*n
    for index in z: res[index] ^= 1
    return res
    
def convert_to_binary_lanczos(z, n):
    res = [0]*n
    for i in range(n):
        if (z >> n - i - 1)&1: res[i] = 1
    return res
    
def create_solution(pairs, null_space, n, len_primes, primes, f_x, m0, m1, inert, f_prime_sq, leading, f_prime_eval, u, LOG_PATH):
    f_norm = 0
    tmp = 1
    for x in f_x:
        f_norm += x*x*tmp
        tmp *= leading
    f_norm = int(math.sqrt(f_norm))+1
    fd = int(pow(len(f_x)-1,1.5))+1
    S = 0
    
    x = f_prime_eval
    x = x*create_rational(null_space, n, len_primes, primes, pairs)%n
    
    rational_square = [i for i in f_prime_sq]
    for k in range(len(null_space)):
        if null_space[k]:
            rational_square = div_poly(poly_prod(rational_square, pairs[k][0]), f_x)
            S += pairs[k][6]
            
    x = x*pow(leading, S>>1, n)%n
            
    coeff_bound = [fd*pow(f_norm, len(f_x)-1-i)*pow(2*(leading*u)*f_norm, S>>1) for i in range(len(f_x)-1)]
    
    y = square_root(f_x, rational_square, inert, m0, m1, leading, max(coeff_bound), LOG_PATH)
    y = y*pow(m1, S>>1, n)%n
    
    return x, y
    
def create_rational(null_space, n, len_primes, primes, pairs):
    x = 1
    vec = [0]*(len_primes)

    for z in range(len(null_space)):
        if null_space[z]:
            for large_prime in pairs[z][4]: x = x*large_prime%n
            cmp = pairs[z][3]
            for p in range(len_primes):
                tmp = primes[p]
                tmp2 = 0
                while not cmp%tmp:
                    tmp *= primes[p]
                    tmp2 += 1
                vec[p] += tmp2

    for i in range(len_primes): x = x*pow(primes[i], vec[i]>>1, n)%n

    return x
    
def create_solution_couveignes(pairs, null_space, n, len_primes, primes, f_x, f_prime, m0, m1, f_prime_sq, leading,
                               f_prime_eval, d, inert_set, zero, delta, u):
    f_norm = 0
    for x in f_x:
        f_norm += x*x
    f_norm = int(math.sqrt(f_norm))+1
    S = 0

    x = f_prime_eval*create_rational(null_space, n, len_primes, primes, pairs)%n
    norm_vec = [0]*(len(primes)+1)
    prod = [1]*d

    for k in range(len(null_space)):
        if null_space[k]:
            for i in range(d): prod[i] *= my_norm(evaluate(pairs[k][0], zero[i]))
            if pairs[k][1] < 0: norm_vec[0] += 1
            for p in range(len(primes)):
                tmp = primes[p]
                while not pairs[k][1]%tmp:
                    tmp *= primes[p]
                    norm_vec[p+1] += 1
            S += pairs[k][6]

    x = x*pow(leading, S>>1, n)%n
    coeffs = [0]*d
    for i in range(d):
        for k in range(d): coeffs[i] += delta[k][i]*(isqrt(prod[k]) + 1)
        
    #goal = 2*max([abs(c) for c in coeffs])
        
    fd = int(pow(d,1.5))+1
    coeff_bound = [fd*pow(f_norm, d-i)*pow(2*u*f_norm, S>>1) for i in range(len(f_x)-1)]
    
    goal = max(coeff_bound)

    sqrt_set, bounds_set = [], []

    target = pow(2, (int(math.sqrt(math.log2(goal)))))
    P = 1
    p = 0
    while P < goal and p < len(inert_set):
        P *= inert_set[p]
        bounds_set.append(inert_set[p])
        p += 1

    while P < goal:
        tmp = random.randint(target//1000,target*1000)
        if m1%tmp and is_prime(tmp) and irreducibility(f_x, tmp) and tmp not in inert_set:
            inert_set.append(tmp)
            P *= tmp
            bounds_set.append(tmp)
            
    iteration = 0
    for p in bounds_set:
        rational_square = [i%p for i in f_prime_sq]
        large = 1
        norm_sq = power(rational_square, f_x, p, (pow(p, d)-1)//(p-1))[0]%p
        for k in range(len(null_space)):
            if null_space[k]:
                rational_square = div_poly_mod(poly_prod_mod(rational_square, pairs[k][0], p), f_x, p)
                for large_prime in pairs[k][2]: large = large*large_prime[0]%p
                norm_sq = norm_sq*pairs[k][1]%p
                
        # The norm of an algebraic number in alpha is equal to f(a,b)/c_d
        # The norm of an algebraic number in omega is equal to f(a,b)*c_d^(d-1) = c_d^d * the norm of that same algebraic number but in alpha
        # The norm is computed as N(q(omega)) = pow(q, (p^d-1)//(p-1), p)
                
        while len(rational_square) < d: rational_square = [0]+rational_square
        
        root = compute_root_mod(rational_square, f_x, p)
        
        N = power(f_prime, f_x, p, (pow(p, d)-1)//(p-1))[0]
        
        if (norm_vec[0]>>1)&1: N = -N%p
        for i in range(len(primes)):
            N = N*pow(primes[i], norm_vec[i+1]>>1, p)%p
        
        N = N*large*pow(leading, (d-1)*(S>>1), p)%p
        
        if power(root, f_x, p, (pow(p, d)-1)//(p-1))[0] != N: root = [-i%p for i in root]

        sqrt_set.append(root)
        iteration += 1
        sys.stdout.write('\r'+str(iteration)+"/"+str(len(bounds_set))+" square roots computed")
        
    r = 0
    rest = 0
    y = 0
    for i in range(len(sqrt_set)):
        while len(sqrt_set[i]) < d: sqrt_set[i] = [0]+sqrt_set[i]
        sqrt_set[i] = eval_F(leading*m0, m1, sqrt_set[i], d-1)%bounds_set[i]
        tmp = invmod(P//bounds_set[i], bounds_set[i])
        tmp2 = tmp*sqrt_set[i]
        r += tmp2//bounds_set[i]
        rest += (tmp2%bounds_set[i])/bounds_set[i]
        y = (y+tmp2*(P//bounds_set[i])%n)%n
        
    r = r+round(rest)
    y = pow(m1, S>>1, n)*(y-r*P%n)%n
    return x, y

def compute_factors(pairs_used, vec, n, primes, g, g_prime, g_prime_sq, g_prime_eval, m0, m1, leading_coeff, d,
                    inert_set, zeros, delta, M, FLAG_SQUARE_ROOT_COUVEIGNES, time_1, LOG_PATH):
    if FLAG_SQUARE_ROOT_COUVEIGNES:
        x, y = create_solution_couveignes(pairs_used, vec, n, len(primes), primes, g, g_prime, m0, m1, g_prime_sq,
                                         leading_coeff, g_prime_eval, d, inert_set, zeros, delta, M<<1)
    
    else:
        x, y = create_solution(pairs_used, vec, n, len(primes), primes, g, m0, m1, inert_set[-1], g_prime_sq,
                              leading_coeff, g_prime_eval, M<<1, LOG_PATH)

    if x != y and math.gcd(x-y, n) != 1 and math.gcd(x+y, n) != 1:
        print_final_message(x, y, n, time_1, LOG_PATH)
        return True, 1

    return False, 0
    
def print_final_message(x, y, n, time_1, LOG_PATH):
    time_2 = datetime.now()
    log.skip_line(LOG_PATH)
    log.write_log(LOG_PATH, "found factor : "+str(math.gcd(x-y,n)))
    log.write_log(LOG_PATH, "found factor : "+str(math.gcd(x+y,n)))
    log.skip_line(LOG_PATH)
    log.write_log(LOG_PATH, "null space found in "+format_duration(time_2-time_1)+".\n")