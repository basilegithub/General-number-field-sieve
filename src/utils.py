# This is the file for the small useful functions

import random, math

## number theory functions

def lowest_set_bit(a):
    b = (a & -a)
    low_bit = -1
    while (b):
        b >>= 1
        low_bit += 1
    return low_bit

# Compute the largest integer m such that m² < n < (m+1)²
# Combines a binary search to find in O(log(log(n)) a good starting point for Newton algorithm, and the Newton algorithm
def isqrt(n):
    if n < 17: b, a = (n+1)>>1, n
    else:
        s, S = 1,1
        while (S<<4) < n:
            b, a, c = 16, 4, 2
            while b*S < n: b, a, c = b*b, b, a
            s *= c
            S = s*s
        if (S<<2) < n: b, a = (s<<2), n
        else: b, a = s<<1, n
    while b < a: a, b = b, (b*b+n)//(b<<1)
    return a

# Computes a^(-1) (mod m)
def invmod(a, m):
    (r, u, R, U) = (a, 1, m, 0)
    while R:
        q = r//R
        (r, u, R, U) = (R, U, r - q *R, u - q*U)
    return u%m
    
def compute_legendre_character(a, n):
    a = a%n
    t = 1
    while a:
        while not a&1:
            a = a>>1
            if n%8 == 3 or n%8 == 5: t = -t
        a, n = n, a
        if a%4 == n%4 and n%4 == 3: t = -t
        a = a%n
    if n == 1: return t
    return 0
        
# Compute x such that x² = n (mod p)
def compute_sqrt_mod_p(n, p):
    n %= p
    if n == 1 : return 1
    P = p-1
    z = int(random.randint(2, P))
    while compute_legendre_character(z, p) != -1:
        z = int(random.randint(2, P))
    r = 0
    while not P&1:
        P >>= 1
        r += 1
    s = P
    generator = pow(z, s, p)
    lbd = pow(n, s, p)
    omega = pow(n, (s+1)>>1, p)

    while True:
        if not lbd: return 0
        if lbd == 1: return omega
        for m in range(1, r):
            if pow(lbd, 1<<m, p)==1: break

        tmp = pow(2, r-m-1, p-1)
        lbd = lbd*pow(generator, tmp<<1, p)%p
        omega = omega*pow(generator, tmp, p)%p
        
# Use chinese remainder theorem to generate some sieve polynomial coefficients
def CRT(moduli, a ,second_part):
    sum = 0
    for i in range(len(moduli)):
        sum = (sum+moduli[i]*second_part[i])%a
    return sum      
    
def is_prime(n):
    if n == 2 or n == 3:
        return True
    if not n&1 or n==1:
        return False
    r, s = 0, n - 1
    while not s&1:
        r += 1
        s >>= 1
    for _ in range(0, 50):
        a = random.randrange(2, n - 1)
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(0, r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True
    
# Fast but not always true Fermat primality test
# There are very few known Fermat primes, so it is most likely almost always correct
# Returns true is the number is not prime
def fermat_primality(n): return pow(2, n-1, n)-1

def dickman(x, table):
    k,res = math.ceil(x),0
    delta = k-x
    if k-1 > len(table): table = get_dickman_table(k)
    tmp = 1
    for i in range(len(table[k-2])):
        res += table[k-2][i]*tmp
        tmp *= delta

    return res, table

def get_dickman_table(k):
    coeffs = [[1-math.log(2)]+[1/(i*(1<<i)) for i in range(1, 30)]]
    for i in range(3, k+1):
        new = [0]*30
        for u in range(1, 30):
            c = 0
            for j in range(u): c += coeffs[-1][j]/(u*pow(i, u-j))
            new[u] = c

        c = 0
        for j in range(1, len(new)): c += new[j-1]/(j+1)
        new[0] = c/(i-1)
        coeffs.append(new)

    return coeffs
    
def my_norm(z):
    if z.real == 0: return math.ceil(abs(z.imag))
    if z.imag == 0: return math.ceil(abs(z.real))

    return isqrt(math.ceil(abs(z.imag))**2 + math.ceil(abs(z.real))**2)+1

def fac(n):
    res = 1
    for i in range(2, n+1): res *= i
    return res
    
def binom(k, n):
    res = 1

    for i in range(n-k+1, n+1): res *= i
    for i in range(2, k+1): res //= i

    return res
        
## polynomial operations functions

def poly_div_bin(poly_a, poly_b):
    remainder = poly_a
    quotient = 0
    deg_a = poly_a.bit_length() - 1
    deg_b = poly_b.bit_length() - 1

    for j in reversed(range(deg_a - deg_b + 1)):
        if remainder & (1 << (deg_b + j)):
            quotient |= (1 << j)
            remainder ^= (poly_b << j)

    return quotient, remainder

def poly_prod_bin(poly_a, poly_b):
    res = 0

    for i in range(poly_a.bit_length()-1, 0, -1):
        if poly_a >> i & 1: res ^= poly_b
        res <<= 1

    if poly_a & 1: res ^= poly_b

    return res
    
def poly_gcd_bin(poly_a, poly_b):
    A = poly_a
    B = poly_b
    
    while B != 0:
        (Q, R) = poly_div_bin(A, B)
        A, B = B, R
        
    return A
    

## linear algebra functions
    
def add_vector(a, b):
    new = [a[i]^b[i] for i in range(len(a))]
    return new
    
def sparse_dot_prod(lbd, x):
    res = 0
    for i in lbd: res ^= x[i]
    
    return res
    
def sparse_multiply(A, b):
    new = [0]*len(b)
    for i in range(len(A)):
        tmp = 0
        for j in A[i]: tmp ^= b[j]
        new[i] = tmp
    return new
    
def dense_multiply(A, B):
    C = [0]*len(A)
    for i in range(len(A)):
        tmp = 0
        for j in range(len(B)):
            if (A[i] >> len(B)-j-1)&1:
                    tmp ^= B[j]
        C[i] = tmp
    return C
    
def transpose_dense(A, N):
    dim1_len = N
    B = [0]*dim1_len
    
    for i in range(len(A)-1):
        tmp = A[i]
        for j in range(dim1_len-1, -1, -1):
            if tmp&1: B[j] ^= 1
            B[j] <<= 1
            tmp >>= 1
            
    tmp = A[-1]
    for j in range(dim1_len-1, -1, -1):
        if tmp&1: B[j] ^= 1
        tmp >>= 1
        
    return B
    
def transpose_sparse(A, n):
    A_transpose = [set() for _ in range(n)]
    
    for i in range(len(A)):
        for index in A[i]:
            A_transpose[index].add(i)
        
    return A_transpose
    
def identity(n):
    new = [(1<<i) for i in range(n-1, -1, -1)]
    return new
    
def concatenate(A, B, N):
    return [(A[i]<<N)^B[i] for i in range(len(A))]


## Statistics function

def central(k,x): return pow(x, k/2-1)/(math.exp(x/2)*pow(2, k/2)*math.gamma(k/2))

def non_central(k, l ,x):
    if x <= 0: return 0
    res = 0
    for i in range(100): res += central(k+2*i, x)*pow(l/2, i)/fac(i)
    
    return res*math.exp(-l/2)
        
        
## time formating functions

def format_duration(delta):
    hours = delta.seconds//3600
    minutes = delta.seconds//60 - hours*60
    seconds = delta.seconds%60
    return str(delta.days)+" days, "+str(hours)+" hours, "+str(minutes)+" minutes, "+str(seconds)+" seconds"