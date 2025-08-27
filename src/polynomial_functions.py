# This file contains all the functions required to handle the NFS polynomials
# Similar to utils, but there are so much functions I prefer a dedicated file

from utils import *
import math

## Creation of polynomials

def coefficients(n, m, d):
    tmp1, tmp2 = n, 1
    coeff = [0]*(d+1)
    i = 0
    for k in range(d+1):
        tmp3 = tmp1%(tmp2*m)
        coeff_k = (tmp3//tmp2)
        coeff[d-k] = coeff_k
        tmp1 -= tmp3
        tmp2 *= m
    return coeff

## Operations on one polynomial
    
def get_derivative(f):
    res = [0]*(len(f)-1)
    for i in range(len(f)-1):
        res[i] = (len(f)-1-i)*f[i]
    return res
    
def new_coeffs(f, x):
    b = 1
    tmp = [i for i in f]
    for i in range(len(f)-1):
        tmp[i] *= b
        b *= x
    tmp[-1] *= b
    return tmp
    
def power(poly, f, p, exp):
    if exp == 1: return poly

    tmp = power(poly, f, p, exp>>1)
    tmp = poly_prod(tmp, tmp)

    if exp&1:
        tmp = poly_prod(tmp, poly)
        return div_poly_mod(tmp, f, p)
    
    else: return div_poly_mod(tmp, f, p)
    
def shift(poly,k):
    res = [i for i in poly]
    for i in range(len(poly)-1):
        for j in range(i+1, len(poly)):
            res[j] += binom(j-i, len(poly)-1-i)*pow(k, j-i)*poly[i]
    return res
    
## Characteristics of one polynomial
    
def irreducibility(f, p):
    g = [1,0]
    for i in range((len(f)-1)//2+1):
        g = power(g, f, p, p)
        tmp2 = [u for u in g]
        if len(tmp2) == 1: tmp2 = [-1, tmp2[0]]
        else: tmp2[-2] -= 1
        tmp = gcd_mod(f, tmp2, p)
        if len(tmp) > 1: return False
    return True
    
def fast_roots(f, p):
    tmp_f = [i%p for i in f]
    for k in range(len(tmp_f)):
        if tmp_f[k]:
            tmp_f = tmp_f[k:]
            break

    roots = []
    tmp = [0]*len(f)
    for i in range(len(f)): tmp[i] = eval_mod(tmp_f, i, p)
    for q in range(1, len(tmp_f)):
        for k in range(len(tmp_f)-1, q-1, -1): tmp[k] -= tmp[k-1]

    res = tmp[0]
    if not res: roots.append(0)
    for q in range(1, p):
        res = (res+tmp[1])%p
        if not res:
            roots.append(q)
            if len(roots) == len(f)-1: return roots
        for k in range(1, len(tmp)-1): tmp[k] += tmp[k+1]
    return roots
    
    
def gcd_mod(f, poly, p):
    while poly != [0]*len(poly):
        (f,poly) = (poly, div_poly_mod(f, poly, p))

    return f
    
def find_roots_poly(f, p):
    tmp_f = [i%p for i in f]
    for k in range(len(f)):
        if tmp_f[k]:
            tmp_f = tmp_f[k:]
            break

    r = []
    tmp = [1,0]
    g = [1]
    tmp_p = p
    while tmp_p>1:
        if tmp_p&1:
            g = div_poly_mod(poly_prod(g, tmp), tmp_f, p)
        tmp = div_poly_mod(poly_prod(tmp, tmp), tmp_f, p)
        tmp_p >>= 1

    g = div_poly_mod(poly_prod(g, tmp), tmp_f, p)
    if len(g) == 1: g = [-1, g[0]]
    else: g[-2] -= 1
    g = gcd_mod(f, g, p)
    if g[-1] == 0:
        r.append(0)
        del g[-1]
    return r+roots(g, p)
    
def roots(g, p):
    if len(g) == 1: return []
    if len(g) == 2: return [-g[1]*invmod(g[0], p)%p]
    if len(g) == 3:
        tmp = (g[1]*g[1]-4*g[0]*g[2])%p
        if tmp == 0: return [-g[1]*invmod(2*g[0], p)%p]
        if compute_legendre_character(tmp, p) == -1: return []
        tmp = compute_sqrt_mod_p(tmp, p)*invmod(2*g[0], p)%p
        return [(-g[1]*invmod(g[0]<<1, p)+tmp)%p, (-g[1]*invmod(g[0]<<1, p)-tmp)%p]
    
    h = [1]
    while len(h) == 1 or h == g:
        a = random.randint(0, p-1)
        h = power([1, a], g, p, (p-1)>>1)
        for k in range(len(h)):
            if h[k]:
                h = h[k:]
                break
        h[-1] -= 1
        h = gcd_mod(h, g, p)
    r = roots(h, p)
    h = quotient_poly_mod(g, h, p)
    return r+roots(h, p)
    
def quadratic_residue(poly, f, p):
    return power(poly, f, p, (pow(p, len(f)-1)-1)>>1)[0]
    
def get_complex_roots(f):
    rho = 0
    for i in range(1, len(f)): rho += abs(f[i])
    rho /= abs(f[0])
    rho = max(1,rho)
    d = len(f)-1
    roots = []
    for i in range(d): roots.append(rho*pow(math.cos(2*math.pi/d)+math.sin(2*math.pi/d)*1j, i))
    next_roots = [None]*d

    for _ in range(1000):
        for i in range(d):
            prod = f[0]
            for k in range(d):
                if k != i: prod *= (roots[i]-roots[k])
            next_roots[i] = roots[i]-evaluate(f, roots[i])/prod
        roots = [i for i in next_roots]
    return roots
    
## Operations between two polynomials
    
def poly_prod(a, b):
    res = [0]*(max(len(a), len(b))+min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] += a[i]*b[j]

    return res
    
def poly_prod_mod(a, b, p):
    res = [0]*(max(len(a), len(b)) + min(len(a), len(b))-1)

    for i in range(len(a)):
        for j in range(len(b)):
            res[i+j] = (res[i+j]+a[i]*b[j])%p

    return res
    
def div_poly(a,b):
    remainder = [i for i in a]
    difference = len(a)-len(b)+1
    leading_coeff = b[0]

    for j in range(difference):
        quotient = -remainder[j]//leading_coeff
        for k in range(len(b)):
            remainder[j+k] += quotient*b[k]
            
    for k in range(len(remainder)):
        if remainder[k]: return remainder[k:]
        
    return [0]
    
def div_poly_mod(a, tmp_b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in tmp_b]
    
    #print(remainder, b)
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = invmod(-b[0], p)
    for j in range(difference):
        if remainder[j]:
            quotient = remainder[j]*coeff%p
            remainder[j] = 0
            for k in range(1,len(b)): remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(remainder)):
        if remainder[k]: return remainder[k:]
        
    return [0]
    
def quotient_poly_mod(a, b, p):
    remainder = [i%p for i in a]
    b = [i%p for i in b]
    
    while not b[0]: del b[0]
    
    difference = len(a)-len(b)+1
    coeff = invmod(-b[0], p)
    res = [0]*difference

    for j in range(difference):
        quotient = remainder[j]*coeff%p
        res[j] = -quotient
        for k in range(len(b)):
            remainder[j+k] = (remainder[j+k]+quotient*b[k]%p)%p
            
    for k in range(len(res)):
        if res[k]: return res[k:]
        
    return [0]
    
def gcd_mod(f, poly, p):
    while poly != [0]*len(poly):
        (f, poly) = (poly, div_poly_mod(f, poly, p))
    return f
    
## Evaluation of polynomials

def eval_mod(f, x, n):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x
        res %= n

    res += f[-1]
    res %= n

    return res
    
def evaluate(f, x):
    res = 0

    for i in range(len(f)-1):
        res += f[i]
        res *= x

    res += f[-1]

    return res
    
def eval_F(x, y, f, d):
    tmp = 0
    tmp2 = 1

    for k in range(d):
        tmp += f[k]*tmp2
        tmp *= x
        tmp2 *= y

    tmp += f[d]*tmp2

    return tmp
    
def eval_G(x, y, m):
    return x-y*m