# This file contains all the functions required to do the polynomial selection step

import math
from polynomial_functions import *
from utils import *
import log
from datetime import datetime
import sys
import multiprocessing
from utils_polynomial_selection import *

def poly_search(n,primes,nb_roots,prime_bound,c,M,d,nb_poly_coarse_eval,nb_poly_precise_eval,NB_CPU_POLY_SELECTION,LOG_PATH):
    log.write_log(LOG_PATH, "Starting polynomial search")

    if d == 1:
        return degree_one_poly_selection(n, primes)
    else:
        return Kleinjung_poly_search(n,primes,nb_roots,prime_bound,c,M,d,nb_poly_coarse_eval,nb_poly_precise_eval,NB_CPU_POLY_SELECTION,LOG_PATH)

def degree_one_poly_selection(n, primes):

    best_poly, best_Lnorm = None, None

    m = isqrt(n)
    len_sieve = 10000
    wanted_skew = 1000
    nb_polys = 300

    cpt = 0

    for a in range(n//(2*m)-nb_polys//2, n//(2*m)+nb_polys//2):
        cpt += 1
        sys.stdout.write('\r'+str(cpt)+" polynomials tested")
        
        sieve_array = [1]*len_sieve
        for p in primes:
            if a%p:
                start = (n-a*m)*invmod(a, p)%p
                for i in range(start, len_sieve, p):
                    sieve_array[i] *= p

        maximum = 0
        for i in range(1, len_sieve):
            if abs((m+i)//sieve_array[i] - wanted_skew) < abs((m+i)//sieve_array[maximum] - wanted_skew):
                maximum = i

        prod = sieve_array[maximum]
        M = m + maximum

        b = (n - a*M)//prod

        f_x = [a, b]

        _, s = get_sieve_region(f_x, primes[-1])
        L_score = get_Lnorm(poly_prod(f_x,f_x), s, primes[-1])
        alpha = alpha_score(f_x, primes[:300])

        if best_poly is None:
            best_poly, best_Lnorm = [[a, -b], M, prod, L_score, alpha], L_score

        elif best_Lnorm > L_score:
            best_poly, best_Lnorm = [[a, -b], M, prod, L_score, alpha], L_score

    return best_poly

# Kleinjung first polynomial search algorithm
def Kleinjung_poly_search(n, primes, nb_roots, prime_bound, c, M, d, NB_POLY_COARSE_EVAL, NB_POLY_PRECISE_EVAL, NB_CPU_POLY_SELECTION, LOG_PATH):
    t1 = datetime.now()

    P = []
    polys = []
    for p in primes:
        if p > prime_bound: break
        if p%d == 1: P.append(p)
    a_d = c
    if d >= 4: admax = round(pow(pow(M, 2*d-2)/n, 1/(d-3)))
    else: admax = M

    cpu, sievers = min(NB_CPU_POLY_SELECTION, multiprocessing.cpu_count()), []

    ad_queue = multiprocessing.Queue()
    poly_queue = [multiprocessing.Queue() for _ in range(cpu-1)]

    for _ in range(2*(cpu-1)):
        ad_queue.put(a_d)
        a_d += c

    for i in range(cpu-1):
        pro = multiprocessing.Process(target=worker, args=(ad_queue, poly_queue[i], primes, n, M, d, nb_roots, P))
        sievers.append(pro)
        pro.start()

    cpt = 0
    avg = 0
    while a_d < admax and cpt < NB_POLY_COARSE_EVAL:

        if ad_queue.qsize() < (cpu-1):
            for _ in range(2*(cpu-1)):
                ad_queue.put(a_d)
                a_d += c

        for i in range(cpu-1):

            if cpt >= NB_POLY_COARSE_EVAL:
                        break

            try:
                poly_list = poly_queue[i].get(timeout=0)

                for to_test_poly in poly_list:
                    cpt += 1
                    sys.stdout.write('\r'+str(cpt)+" polynomials tested")

                    L_score = to_test_poly[-1]

                    avg += L_score
                    if not len(polys):
                        polys.append(to_test_poly)
                    else:
                        if L_score < polys[-1][3]:
                            a = 0
                            b = len(polys)-1
                            tmpu = (a+b)>>1
                            while a <= b:
                                if polys[tmpu][3] < L_score: a = tmpu+1
                                else: b = tmpu-1
                                tmpu = (a+b)>>1
                            polys.insert(a, to_test_poly)
                        elif len(polys) < NB_POLY_PRECISE_EVAL:
                            polys.append(to_test_poly)
                        if len(polys) > NB_POLY_PRECISE_EVAL: del polys[-1]

                    if cpt >= NB_POLY_COARSE_EVAL:
                        break

            except:
                pass

    while len(poly_queue):
        del poly_queue[0]

    for p in sievers:
        p.terminate()
        p.join()

    ad_queue.close()

    t2 = datetime.now()
    print("")
    log.write_log(LOG_PATH, "Polynomial search done in "+format_duration(t2-t1)+".\n")
    log.write_log(LOG_PATH, str(cpt)+" polynomials created, "+str(len(polys))+" kept for ranking")
    log.write_log(LOG_PATH, "Average L2 score = "+str(avg/cpt))
    log.write_log(LOG_PATH, "Ranking polynomials")

    return select_best_poly_candidate(polys, primes)

def select_best_poly_candidate(polys, primes):
    best_poly = None
    best_E = None
    table = get_dickman_table(10)

    for i in range(len(polys)):
        x_limit,s = get_sieve_region(polys[i][0], primes[-1])
        alpha = alpha_score(polys[i][0], primes[:300])
        polys[i].append(alpha)
        E_score = get_Epscore(polys[i][0], [polys[i][2], -polys[i][1]], alpha, primes[-1], x_limit, round(primes[-1]/math.sqrt(s)), table)
        if best_E == None or E_score > best_E:
            best_poly = polys[i]
            best_E = E_score
    return best_poly

def practical_d(n):
    min,mind = None,None
    for d in range(2,10):
        f = d*math.log(d)+math.sqrt(pow(d*math.log(d), 2)+4*math.log(n)*math.log(math.log(n)/(d+1))/d+1)
        if min == None or f < min: min, mind = f, d
    return mind

## Workers definition

def worker(ad_queue, output_queue, primes, n, M, d, nb_roots, P):
    while True:
        
        try:
            a_d = ad_queue.get(timeout=0)

            tmp = a_d
            for p in primes:
                while not tmp%p: tmp//= p

            if tmp == 1: # If a_d is not primes[-1] smooth

                mw = math.ceil(pow(n/a_d, 1/d))
                ad1max = round(M*M/mw)
                if d > 2: ad2max = pow(pow(M, 2*d-6)/pow(mw, d-4), 1/(d-2))
                else: ad2max = M

                Q = []
                roots = []
                for p in P:
                    if not a_d%p: continue

                    f = [a_d]+[0]*d
                    f[-1] = (-n)%p # Construct polynomial a_d*x^d - n (mod r)
                    root = find_roots_poly(f, p)
                    if len(root) > 0:
                        Q.append(p)
                        roots.append(root)

                if len(roots) >= nb_roots:

                    combinations = prime_combinations_with_indices(Q, nb_roots, ad1max)

                    for set in combinations:
                        Q_used = []
                        prod = 1
                        for i in range(nb_roots):
                            Q_used.append(Q[set[i]])
                            prod *= Q[set[i]]

                        root_used = [roots[set[i]] for i in range(nb_roots)]
                        for i in range(nb_roots): # Do some CRT
                            x = prod//Q_used[i]
                            tmp2 = x*invmod(x, Q_used[i])
                            for j in range(d): root_used[i][j] = root_used[i][j]*tmp2%prod

                        m0 = mw+(-mw)%prod
                        e = compute_e(m0, root_used, nb_roots, prod, a_d, n, d)
                        f, f0 = compute_f(n, a_d, m0, d, prod, root_used, nb_roots, e)

                        epsilon = ad2max/m0
                        array1 = create_first_array(nb_roots, f0, f, d)
                        len_vec = nb_roots>>1
                        array2 = create_second_array(nb_roots, len_vec, d, f)

                        to_return_array = []
                        
                        min = 0
                        for j in range(len(array2)):
                            while min < len(array1) and array2[j][0]-epsilon > array1[min][0]: min += 1
                            if min == len(array1): break
                            z = min
                            while z < len(array1) and abs(array2[j][0]-array1[z][0]) < epsilon:
                                tmp = [poly(m_mu(m0, root_used, array1[z][1]+array2[j][1], nb_roots), prod, a_d, n, d),
                                       m_mu(m0, root_used, array1[z][1]+array2[j][1], nb_roots),
                                       prod]
                                tmp[0], tmp[1] = local_opt(tmp[0], [prod,-tmp[1]],primes[-1])
                                _, s = get_sieve_region(tmp[0], primes[-1])
                                L_score = get_Lnorm(poly_prod(tmp[0], tmp[0]), s, primes[-1])
                                tmp.append(L_score)

                                to_return_array.append(tmp)

                                z += 1

                        output_queue.put(to_return_array)

        except:
            pass