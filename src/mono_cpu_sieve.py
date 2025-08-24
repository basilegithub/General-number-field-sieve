# This file contains the functions that realize the sieve step in the case where only one cpu is used

from datetime import datetime
import sys
import math
import log
from polynomial_functions import *
import sieve
import find_smooth
from relations import *
from utils import format_duration

def find_relations(f_x, leading_coeff, g, primes, R_p, Q, B_prime, divide_leading, prod_primes, pow_div, pairs_used,
                   const1, const2, logs, m0, m1, M, d, n, flag_use_batch_smooth_test, LOG_PATH):
    fp, pf = {}, {}
    parent_fp, parent_pf = {}, {}
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
    
    b = 1
        
    sys.stdout.write('\r'+"0/("+str(V)+"+10) relations found")
    while len(pairs_used) < V+10:
        
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
            smooth = find_smooth.batch_smooth_test(pairs, prod_primes, const1, const2, div)

        for i, pair in enumerate(pairs):
            if flag_use_batch_smooth_test:
                z = smooth[i]
            else:
                z = find_smooth.trial(pair,primes,const1,const2,div)
            if z[0]:
                tmp = [u for u in pair]
                for p in divide_leading: tmp.append(not pair[5][0]%p)
                if z[2] == 1 and z[4] == 1:
                    pairs_used.append(tmp)
                    full_found += 1

                elif z[4] > 1 and z[2] == 1:
                    pairs_used, fp, graph_fp, size_fp, parent_fp, cycle_len, full_found, partial_found_fp = handle_large_fp(tmp,z,pairs_used,fp,graph_fp,size_fp,parent_fp,g,divide_leading,cycle_len,full_found,partial_found_fp)

                elif z[4] == 1 and z[2] > 1:
                    pairs_used, pf, graph_pf, size_pf, parent_pf, cycle_len, full_found, partial_found_pf = handle_large_pf(tmp, z, pairs_used, pf, graph_pf, size_pf, parent_pf, g, divide_leading, cycle_len, full_found, partial_found_pf)

        sys.stdout.write('\r'+"b = "+str(b)+" "+str(len(pairs_used))+"/("+str(V)+"+10) ; full relations = "+str(full_found)+" | partial found fp = "+str(partial_found_fp)+" ("+str(size_fp)+") | partial found pf = "+str(partial_found_pf)+" ("+str(size_pf)+")")
        b += 1
        
    print("\n")
    
    time_2 = datetime.now()
    
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