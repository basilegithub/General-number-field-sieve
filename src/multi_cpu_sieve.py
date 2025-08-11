# This file contains the functions that realize the sieve step when multiprocessing is on (NB_CPU > 1)

from datetime import datetime
import time
import sys, os
import math
import log
import sieve
import find_smooth
from relations import *
from utils import format_duration
import multiprocessing
                
def siever(b_queue, tmp_pairs, f_x, primes, R_p, prod_primes, m0, m1, d, b, M, logs, offset, leading_coeff, divide_leading,
           pow_div, const_1, const_2, flag_use_batch_smooth_test):
    
    while True:
        try:
            b = b_queue.get(timeout=0)

            div = 1
            for q in range(len(divide_leading)):
                p = divide_leading[q]
                if not b%p:
                    tmp = p
                    while not b%tmp and tmp <= pow_div[q]:
                        div *= p
                        tmp *= p

            pairs = sieve.sieve(M,f_x,primes,R_p,m0,m1,d,b,logs,offset,leading_coeff,div)

            if flag_use_batch_smooth_test:
                smooth = find_smooth.batch_smooth_test(pairs, prod_primes, const_1, const_2, div)

            else:
                smooth = [0]*len(pairs)
                for i, pair in enumerate(pairs):
                    smooth[i] = find_smooth.trial(pair,primes,const_1,const_2,div)

            tmp_pairs.put((pairs, smooth))

        except:
            pass

def find_relations(f_x, leading_coeff, g, primes, R_p, Q, B_prime, divide_leading, prod_primes, pow_div, pairs_used,
                   const1, const2, logs, m0, m1, M, d, n, flag_use_batch_smooth_test, LOG_PATH, NB_CPU):
    fp, pf = {}, {}
    connected_components_fp, connected_components_pf = {}, {}
    node_component_fp, node_component_pf = {}, {}
    index_component_fp, index_component_pf = 0, 0
    full_found = 0
    partial_found_fp, partial_found_pf = 0, 0
    size_fp, size_pf = 0, 0
    graph_fp, graph_pf = {}, {}
    offset = 15+math.log2(const1)
    cycle_len = [0]*10
    
    V = 3+len(primes)+B_prime+len(Q)+len(divide_leading)
    
    log.write_log(LOG_PATH, "sieving...")
    log.write_log(LOG_PATH, "need to find at least "+str(V+10)+" relations\n")
    
    time_1 = datetime.now()

    cpu, sievers = min(NB_CPU, multiprocessing.cpu_count()), []
    
    b = 1
    b_queue = multiprocessing.Queue()
    tmp_pairs_queues = [multiprocessing.Queue() for _ in range(cpu-1)]

    for _ in range(2*(cpu-1)):
        b_queue.put(b)
        b += 1

    for i in range(cpu-1):
        pro = multiprocessing.Process(target=siever,args=(b_queue, tmp_pairs_queues[i], f_x, primes, R_p, prod_primes, m0, m1,
                                                          d, b, M, logs, offset, leading_coeff, divide_leading, 
                                                          pow_div, const1, const2, flag_use_batch_smooth_test))
        sievers.append(pro)
        pro.start()
        
    sys.stdout.write('\r'+"0/("+str(V)+"+10) relations found")
    while len(pairs_used) < V+10:

        if b_queue.qsize() < (cpu-1):
            for _ in range(2*(cpu-1)):
                b_queue.put(b)
                b += 1

        for k in range(cpu-1):
            try:
                pairs, smooth = tmp_pairs_queues[k].get(timeout=0)
                for i, pair in enumerate(pairs):
                    z = smooth[i]
                    if z[0]:
                        tmp = [u for u in pair]
                        for p in divide_leading: tmp.append(not pair[5][0]%p)
                        if z[2] == 1 and z[4] == 1:
                            pairs_used.append(tmp)
                            full_found += 1

                        elif z[4] > 1 and z[2] == 1:
                            pairs_used, fp, graph_fp, size_fp, index_component_fp, cycle_len, full_found, partial_found_fp = handle_large_fp(tmp,z,pairs_used,fp,graph_fp,size_fp,connected_components_fp,node_component_fp,index_component_fp,g,divide_leading,cycle_len,full_found,partial_found_fp)

                        elif z[4] == 1 and z[2] > 1:
                            pairs_used, pf, graph_pf, size_pf, index_component_pf, cycle_len, full_found, partial_found_pf = handle_large_pf(tmp, z, pairs_used, pf, graph_pf, size_pf, connected_components_pf, node_component_pf, index_component_pf, g, divide_leading, cycle_len, full_found, partial_found_pf)

                sys.stdout.write('\r'+"b = "+str(b)+" "+str(len(pairs_used))+"/("+str(V)+"+10) ; full relations = "+str(full_found)+" | partial found fp = "+str(partial_found_fp)+" ("+str(size_fp)+") | partial found pf = "+str(partial_found_pf)+" ("+str(size_pf)+")")

                if len(pairs_used) >= V+10: break

            except:
                pass
        
    print("\n")
    
    time_2 = datetime.now()

    while len(tmp_pairs_queues):
        del tmp_pairs_queues[0]

    for p in sievers:
        p.terminate()
        p.join()

    b_queue.close()
    
    log.write_log(LOG_PATH, "sieving done in "+format_duration(time_2-time_1)+".\n")
    log.write_log(LOG_PATH, "Distribution of cycle length:")
    log.write_log(LOG_PATH, "1-cycle: "+str(full_found))
    for i in range(9):
        log.write_log(LOG_PATH, str(i+2)+"-cycle: "+str(cycle_len[i]))
    log.write_log(LOG_PATH, "11+-cycle: "+str(cycle_len[-1])+"\n")
    log.write_log(LOG_PATH, str(partial_found_fp)+" partial relations fp found")
    log.write_log(LOG_PATH, str(partial_found_pf)+" partial relations pf found")
    #log.write_log(LOG_PATH, str(len(fp))+" smooth with large prime(s) found")
                
    return pairs_used, V