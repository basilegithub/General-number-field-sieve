# This is the file to find smooth numbers

import math, random
from utils import fermat_primality

# batch smooth test using a product tree
# more efficient than testing for smoothness every candidate : should be prefered
# Note that in this function, we test for smoothness abs(candidates), see QS.sieve_and_batch_smooth function
# Source: https://cr.yp.to/factorization/smoothparts-20040510.pdf (Bernstein)
    
def batch_smooth_test(candidates, prod_primes, cst_1, cst_2, div):
    L = len(candidates)
    # compute the smallest e such that 2**(2**e) >= max(candidates), ie e = int(log2(log2(max(candidates)))
    e = 0
    tmp = 2
    tmp_max = max([abs(candidate[1])//div for candidate in candidates] + [abs(candidate[3]) for candidate in candidates])
    while tmp < tmp_max:
        tmp = tmp*tmp
        e += 1
        
    # build the product tree
    tmp = [abs(candidate[1])//div for candidate in candidates] + [abs(candidate[3]) for candidate in candidates]
    prod_tree = [tmp]
    
    while len(prod_tree[-1]) > 1:
        line = [0]*(len(prod_tree[-1])>>1)
        for i in range(1,len(prod_tree[-1]),2):
            tmp = prod_tree[-1][i]*prod_tree[-1][i-1]
            if tmp <= prod_primes: line[i>>1] = tmp
            else: line[i>>1] = prod_primes+1
        if len(prod_tree[-1]) & 1: line.append(prod_tree[-1][-1])
        prod_tree.append(line)
    
    # update the product tree by computing the adequate remainders

    prod_tree[-1][0] = prod_primes%prod_tree[-1][0]
    for i in range(len(prod_tree)-1):
        for j in range(len(prod_tree[-i-1])-1):
            prod_tree[-i-2][(j<<1)] = prod_tree[-i-1][j]%prod_tree[-i-2][j<<1]
            prod_tree[-i-2][(j<<1)+1] = prod_tree[-i-1][j]%prod_tree[-i-2][(j<<1)+1]
        tmp = len(prod_tree[-i-1])-1
        prod_tree[-i-2][tmp<<1] = prod_tree[-i-1][tmp]%prod_tree[-i-2][tmp<<1]
        if (tmp<<1)+1 < len(prod_tree[-i-2]):
            prod_tree[-i-2][(tmp<<1)+1] = prod_tree[-i-1][tmp]%prod_tree[-i-2][(tmp<<1)+1]

    # test for smoothness
    smooth = [0]*L
    for i in range(L):
        tmp = 0
        j_1, j_2 = prod_tree[0][i], prod_tree[0][L+i]
        while (j_1 or j_2) and tmp < e:
            j_1 = (j_1*j_1)%(abs(candidates[i][1]//div))
            j_2 = (j_2*j_2)%abs(candidates[i][3])
            tmp += 1

        if not j_1 and not j_2: smooth[i] = [True, 1, 1, 1, 1]

        else:
            tmp_1 = (abs(candidates[i][1])//div)//math.gcd(abs(candidates[i][1]//div),j_1)
            tmp_2 = abs(candidates[i][3])//math.gcd(abs(candidates[i][3]),j_2)
            
            if tmp_1 > 1 and tmp_2 > 1: smooth[i] = [False, 1, 1, 1, 1]

            elif tmp_1 > 1:
                if tmp_1 < cst_1:
                    smooth[i] = [True, 1, tmp_1, 1, 1]
                elif tmp_1 < cst_2 and fermat_primality(tmp_1):
                    tmp_1 = pollard_rho(tmp_1)
                    smooth[i] = [True, min(tmp_1), max(tmp_1), 1, 1]
                else:
                    smooth[i] = [False, 1, 1, 1, 1]

            elif tmp_2 > 1:
                if tmp_2 < cst_1:
                    smooth[i] = [True, 1, 1, 1, tmp_2]
                elif tmp_2 < cst_2 and fermat_primality(tmp_2):
                    tmp_2 = pollard_rho(tmp_2)
                    smooth[i] = [True, 1, 1, min(tmp_2), max(tmp_2)]
                else:
                    smooth[i] = [False, 1, 1, 1, 1]



            else: smooth[i] = [False, 1, 1, 1, 1]
            
    return smooth

# Naive smoothness test: trial division by every prime in the base  
# Note that in this function, we test for smoothness (candidate) and not abs(candidate), see QS.sieve_and_smooth
def trial(pair, primes, const1, const2, div):
    large1 = 1
    large2 = 1
    
    result = abs(pair[3])
    for p in primes:
        while not result%p: result //= p
        if result == 1: break
    if result > 1:
        if result > const1: return False, 1, 1, 1, 1
        else: large1 = result
        
    result = abs(pair[1])//div
    for p in primes:
        while not result%p: result //= p
        if result == 1: return True,1,1,1,large1
    if result > 1:
        if large1 > 1: return False,1,1,1,1
        elif result < const1: large2 = result
        else: return False,1,1,1,1
        
    return True, 1, large2, 1, large1
    
    
# Pollard rho algorithm to factor small primes
# This is used to find the two large primes when needed
def pollard_rho(n):
    a = int(random.randint(1,n-3))
    s = int(random.randint(0,n-1))
    x = s
    y = s
    d = 1
    while d == 1:
        e = 1
        X = x
        Y = y
        for k in range(0,100):
            x = (x*x+a)%n
            y = (y*y+a)%n
            y = (y*y+a)%n
            e = e*(x-y)%n
        d = math.gcd(e,n)
    if d == n:
        x = X
        y = Y
        d = 1
        while d == 1:
            x = (x*x+a)%n
            y = (y*y+a)%n
            y = (y*y+a)%n
            d = math.gcd(x-y,n)
    if d == n:
        return pollard_rho(n)
    return d, n//d