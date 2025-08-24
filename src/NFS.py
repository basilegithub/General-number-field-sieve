# This is the file that execute the Quadratic Sieve algorithm

## import libraries

ZQDFEGSRHTDR = 35279032783773322923326661023

import math, time, os
from datetime import datetime
import log
import parse_config
from utils import *
from polynomial_functions import *
from generate_primes import *
import mono_cpu_polynomial_selection
import multi_cpu_polynomial_selection
import mono_cpu_sieve
import multi_cpu_sieve
from relations import *
from build_matrix import *
from gaussian_elimination import *
from wiedemann import *
from block_lanczos import *
import compute_solutions

## Set path to find config file

CONFIG_PATH = os.path.abspath("../config/config.ini")
    
## NFS functions
    
def initialize(n):
    
    d = int(pow(3*math.log(n)/math.log(math.log(n)),1/3))
    first_B = int(math.exp(pow((8/9)*math.log(n),1/3)*pow(math.log(math.log(n)),2/3)))
    B = 10*first_B//70
    
    return d, B
    
def initialize_2(f_x, n, m1, d, primes, leading_coeff):
    pairs_used,R_p,logs,tmp = [],[],[],10*d
    divide_leading,pow_div = [],[]
    for p in primes:
        if not leading_coeff%p:
            divide_leading.append(p)
            u = p
            while not leading_coeff%u: u *= p
            pow_div.append(u//p)
        logs.append(round(math.log2(p)))
        if p > tmp: R_p.append(find_roots_poly(f_x,p))
        else: R_p.append(fast_roots(f_x,p))
        if len(R_p[-1]) == d:
            if d&1: pairs_used.append([[leading_coeff*p],leading_coeff*pow(p,d),[[1,1]],p*m1,[1],[0,p],1])
            else: pairs_used.append([[leading_coeff*p],leading_coeff*p,[[1,1]],p*m1,[1],[0,p],1])
    for p in divide_leading:
        for i in range(len(pairs_used)): pairs_used[i].append(True)
        
    B_prime = sum([len(i) for i in R_p]) + 1
    
    return pairs_used, R_p, logs, divide_leading, pow_div, B_prime
    
def initialize_3(n, f_x, f_prime, const1, leading_coeff):
    k = 3*n.bit_length()
    Q = []
    q = const1+1
    if not q%2: q += 1
    while len(Q) < k:
        if is_prime(q) and leading_coeff%q:
            if not n%q: return q, n//q
            else:
                tmp = find_roots_poly(f_x,q)
                for r in tmp:
                    if eval_mod(f_prime,r,q): Q.append([q,r])
        q += 2
        
    return Q, k
    
def NFS(n):
    now = datetime.now()
    
    LOG_PATH = os.path.abspath("../logs/log_"+str(now.year)+str(now.month)+str(now.day)+"_"+str(now.hour)+str(now.minute)+str(now.second)+".txt")
    
    parameters = parse_config.parse_config(CONFIG_PATH)
    
    flag_use_batch_smooth_test = parameters[0].lower() in ["true"]
    flag_gaussian_pivot = parameters[1].lower() in ["true"]
    flag_lanczos = parameters[2].lower() in ["true"]
    flag_square_root_couveignes = parameters[3].lower() in ["true"]
    const = int(parameters[4])
    BLOCK_SIZE = int(parameters[5])
    nb_poly_coarse_eval = int(parameters[6])
    nb_poly_precise_eval = int(parameters[7])
    prime_bound = int(parameters[8])
    nb_roots = int(parameters[9])
    multiplier = int(parameters[10])
    NB_CPU_POLY_SELECTION = int(parameters[11])
    NB_CPU_SIEVE = int(parameters[12])
    
    n = int(n)
    
    d, B = initialize(n)
    
    primes = create_smooth_primes_base(B)
    for p in primes:
        if not n%p: return p,n//p

    prod_primes = math.prod(primes)
        
    const1, const2 = const*primes[-1], const*primes[-1]*primes[-1]
    
    if NB_CPU_POLY_SELECTION == 1:
        f_x,m0,m1,tmp,_ = mono_cpu_polynomial_selection.poly_search(n,primes,nb_roots,prime_bound,multiplier,
                                                                 int(pow(n,1/(d+1))),d,nb_poly_coarse_eval,
                                                                 nb_poly_precise_eval,LOG_PATH)
        
    else:
        f_x,m0,m1,tmp,_ = multi_cpu_polynomial_selection.poly_search(n,primes,nb_roots,prime_bound,multiplier,
                                                                 int(pow(n,1/(d+1))),d,nb_poly_coarse_eval,
                                                                 nb_poly_precise_eval,NB_CPU_POLY_SELECTION,LOG_PATH)

    log.write_log(LOG_PATH, "poly search completed, parameters : m0 = "+str(m0)+" ; m1 = "+str(m1)+" ; d = "+str(d)+"\n")
    
    f_x, m0, M = mono_cpu_polynomial_selection.evaluate_polynomial_quality(f_x,B,m0,m1,primes,LOG_PATH)
    
    leading_coeff = f_x[0]
    zeros_f = get_complex_roots(f_x)
    zeros = [leading_coeff*i for i in zeros_f]
    f_prime = get_derivative(f_x)
    
    g = [1,f_x[1]]
    for i in range(2,len(f_x)): g.append(f_x[i]*pow(leading_coeff,i-1))
    g_prime = get_derivative(g)

    g_prime_sq,g_prime_eval = div_poly(poly_prod(g_prime,g_prime),g),pow(leading_coeff, d-2, n)*eval_F(m0,m1,f_prime,d-1)%n
    
    delta = []
    for r in range(d):
        tmp = []
        for i in range(d):
            delt = 0
            for j in range(d-i):
                if zeros[r].real == 0 or zeros[r].imag == 0 or j == 0: delt += (abs(g[-i-j-2])*pow(leading_coeff,j)*pow(my_norm(zeros_f[r]),j))
                else: delt += (abs(g[-i-j-2])*pow(leading_coeff,j)*pow(my_norm(zeros_f[r]),j))+1
            tmp.append(delt)
        delta.append(tmp)

    
    inert_set = []
    for p in primes:
        if m1%p and leading_coeff%p and irreducibility(g,p): inert_set.append(p)
    
    pairs_used, R_p, logs, divide_leading, pow_div, B_prime = initialize_2(f_x, n, m1, d, primes, leading_coeff)
        
    Q, k = initialize_3(n, f_x, f_prime, B, leading_coeff)
        
    log.write_log(LOG_PATH, "components of factor base : ("+str(len(primes))+","+str(B_prime)+","+str(k)+","+str(len(divide_leading))+")")
    log.write_log(LOG_PATH, "starting with "+str(len(pairs_used))+" free relations")
    
    if NB_CPU_SIEVE <= 1:
        if NB_CPU_SIEVE < 1: print("NB_CPU parameter incorrectly set. Must be > 0. Sieving with 1 CPU.")
        pairs_used, V = mono_cpu_sieve.find_relations(f_x,leading_coeff,g,primes,R_p,Q,B_prime,divide_leading,
                                                      prod_primes,pow_div,pairs_used,const1,const2,logs,m0,m1,
                                                      M,d,n,flag_use_batch_smooth_test,LOG_PATH)

    else:
        pairs_used, V = multi_cpu_sieve.find_relations(f_x,leading_coeff,g,primes,R_p,Q,B_prime,divide_leading,
                                                       prod_primes,pow_div,pairs_used,const1,const2,logs,m0,m1,
                                                       M,d,n,flag_use_batch_smooth_test,LOG_PATH,NB_CPU_SIEVE)

    print("")
    log.write_log(LOG_PATH, "sieving complete, building matrix...")
    if flag_gaussian_pivot:
        matrix = build_dense_matrix(pairs_used, primes, R_p, Q, divide_leading)
        length = len(matrix[0])
        matrix, N, U = siqs_build_matrix_opt(matrix)
        log.write_log(LOG_PATH, "matrix built "+str(len(matrix))+"x"+str(len(pairs_used))+" finding kernel...")
    else:
        matrix = build_sparse_matrix(pairs_used,primes,R_p,Q,divide_leading)
        matrix = transpose_sparse(matrix, V)
        log.write_log(LOG_PATH, "matrix built "+str(len(matrix))+"x"+str(len(pairs_used))+" reducing...")
    if not flag_gaussian_pivot:
        matrix, pairs_used = reduce_sparse_matrix(matrix, pairs_used)
        log.write_log(LOG_PATH, "matrix built "+str(len(matrix))+"x"+str(len(pairs_used))+" finding kernel...")

    time_1 = datetime.now()

    if flag_gaussian_pivot:
        null_space = siqs_solve_matrix_opt(matrix, N, U)
        log.write_log(LOG_PATH, str(len(null_space))+" kernel vectors found")
        for vec in null_space:
            vec = compute_solutions.convert_to_binary(vec, length)
            flag, res = compute_solutions.compute_factors(pairs_used,vec,n,primes,f_x,g,g_prime,g_prime_sq,g_prime_eval,
                                                        m0,m1,leading_coeff,d,inert_set,zeros,delta,M,flag_square_root_couveignes,
                                                        time_1,LOG_PATH)
            if flag: return res

    else:
        if flag_lanczos:
            while True:
                null_space = block_lanczos(matrix,len(pairs_used),BLOCK_SIZE, LOG_PATH)
                
                log.write_log(LOG_PATH, str(len(null_space))+" kernel vectors found")
                for vec in null_space:
                    vec = compute_solutions.convert_to_binary_lanczos(vec, len(pairs_used))
                    flag, res = compute_solutions.compute_factors(pairs_used,vec,n,primes,f_x,g,g_prime,g_prime_sq,g_prime_eval,
                                                             m0,m1,leading_coeff,d,inert_set,zeros,delta,M,flag_square_root_couveignes,
                                                             time_1,LOG_PATH)
                    if flag: return res

        else:
            mini_poly_estim = 1
            while True:
                null_space, mini_poly_estim = wiedemann(matrix,len(pairs_used),BLOCK_SIZE,mini_poly_estim)
                null_space = reduce_null_space_vectors(null_space)
                log.write_log(LOG_PATH, str(len(null_space))+" kernel vectors found")
                for vec in null_space:
                    flag, res = compute_solutions.compute_factors(pairs_used,vec,n,primes,f_x,g,g_prime,g_prime_sq,g_prime_eval,
                                                             m0,m1,leading_coeff,d,inert_set,zeros,delta,M,flag_square_root_couveignes,
                                                             time_1,LOG_PATH)
                    if flag: return res