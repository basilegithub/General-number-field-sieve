# This file contains the functions required to compute the square roots in the NFS algorithm

from polynomial_functions import *
import time
    
def square_root(f, N, p, m0, m1, leading, bound):
    d = len(f)-1
    root = compute_root_mod_Newton(N, f, p)
    print("initial root computed, lifting...")
    root = Newton(root, N, f, p, bound)
    print(str(time.localtime()[3])+":"+str(time.localtime()[4])+":"+str(time.localtime()[5])+" lift done")
    return eval_F(leading*m0, m1, root, d-1)
                
def Newton(root, gamma, f, p, bound):
    modulo = p
    old = [1]
    while modulo < bound<<1:
        modulo *= modulo
        tmp = div_poly_mod(poly_prod(root, root), f, modulo)
        tmp = div_poly_mod(poly_prod(tmp, gamma), f, modulo)
        tmp = [-i%modulo for i in tmp]
        tmp[-1] = (tmp[-1]+3)%modulo
        tmp = div_poly_mod(poly_prod(root, tmp), f, modulo)
        tmp2 = invmod(2, modulo)
        root = [tmp2*i%modulo for i in tmp]
        
        tmp = div_poly_mod(poly_prod(root, gamma), f, modulo)
        for i in range(len(tmp)):
            if tmp[i] > modulo>>1: tmp[i] -= modulo

    return tmp
        
                
def compute_root_mod_Newton(number, f, p):
    n = div_poly_mod(number, f, p)
    for k in range(len(n)):
        if n[k]:
            n = n[k:]
            break
            
    s = pow(p, len(f)-1)-1
    r = 0
    while not s&1:
        s >>= 1
        r += 1
        
    z = [random.randint(0, p-1) for i in range(len(f)-1)]
    while quadratic_residue(z, f, p) != p-1:
        z = [random.randint(0, p-1) for i in range(len(f)-1)]
    
    lbd = power(n, f, p, s)
    omega = power(n, f, p, (s+1)>>1)
    zeta = power(z, f, p, s)
    while True:
        if lbd == [0]*len(lbd): return [0]
        if lbd == [1]: return power(omega, f, p,pow(p, len(f)-1)-2)
        k = 0
        for m in range(1, r):
            k = m
            if div_poly_mod(power(lbd, f, p, 1<<m), f, p)==[1]: break
            
        tmp = 1<<(r-k-1)
        tmp = power(zeta, f, p, tmp)
        omega = div_poly_mod(poly_prod(omega, tmp), f, p)
        
        tmp = 1<<(r-k)
        tmp = power(zeta, f, p, tmp)
        lbd = div_poly_mod(poly_prod(lbd, tmp), f, p)
        
def compute_root_mod(number, f, p):
    n = div_poly_mod(number, f, p)
    for k in range(len(n)):
        if n[k]:
            n = n[k:]
            break
            
    s = pow(p, len(f)-1)-1
    r = 0
    while not s&1:
        s >>= 1
        r += 1
        
    z = [random.randint(0, p-1) for i in range(len(f)-1)]
    while quadratic_residue(z, f, p) != p-1:
        z = [random.randint(0, p-1) for i in range(len(f)-1)]
    
    lbd = power(n, f, p, s)
    omega = power(n, f, p, (s+1)>>1)
    zeta = power(z, f, p, s)
    while True:
        if lbd == [0]*len(lbd): return [0]
        if lbd == [1]: return omega
        k = 0
        for m in range(1, r):
            k = m
            if div_poly_mod(power(lbd, f, p, 1<<m), f, p)==[1]: break
            
        tmp = 1<<(r-k-1)
        tmp = power(zeta, f, p, tmp)
        omega = div_poly_mod(poly_prod(omega, tmp), f, p)
        
        tmp = 1<<(r-k)
        tmp = power(zeta, f, p, tmp)
        lbd = div_poly_mod(poly_prod(lbd, tmp), f, p)