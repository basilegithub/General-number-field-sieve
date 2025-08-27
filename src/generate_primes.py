# This is the file for generating prime numbers

from utils import isqrt

# Erasthotenes sieve algorithm

def create_smooth_primes_base(B):
    m = isqrt(B)+1
    
    liste = [True, False]*(B>>1)
    base=[2]
    
    for k in (k for k in range (3, m, 2) if liste[k-1]):
        base.append(k)
        for i in range (k*k, B, k<<1):
            liste[i-1]=False

    for k in (k for k in range (m, B) if liste[k-1]):
        base.append(k)

    return base