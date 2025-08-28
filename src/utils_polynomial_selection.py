from utils import *
from polynomial_functions import *
import log
from numpy.polynomial.legendre import leggauss

def minimise_Lnom(f, s, B, m0, m1):
    k, res = 1, get_Lnorm(poly_prod(f,f),s,B)
    while k > 0:
        flag = False
        F = shift(f, -k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < res:
            res, f, k, flag, m0 = tmp, [i for i in F], k<<1, True, m0+k*m1

        F = shift(f, k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < res:
            res, f, k, flag, m0 = tmp, [i for i in F], k<<1, True, m0-k*m1
        if not flag: k >>= 1

    return f, m0
    
def evaluate_polynomial_quality(f_x, B, m0, m1, primes, LOG_PATH):
    _, s = get_sieve_region(f_x, B)

    if len(f_x) > 2:

        if (s*abs(f_x[-3])-abs(f_x[-2]))//m0 > 0:
            f_x = root_sieve2(f_x, [m1,-m0], primes[:150], (s*abs(f_x[-3])-abs(f_x[-2]))//m0, (s*s*abs(f_x[-3])-abs(f_x[-1]))//m0)
            _, s = get_sieve_region(f_x, B)
            f_x, m0 = minimise_Lnom(f_x, s, B, m0, m1)

        elif (s*abs(f_x[-2])-abs(f_x[-1]))//m0 > 0:
            f_x = root_sieve(f_x, [m1,-m0], primes[:150], (s*abs(f_x[-2])-abs(f_x[-1]))//m0)
            _, s = get_sieve_region(f_x, B)
            f_x, m0 = minimise_Lnom(f_x, s, B, m0, m1)

    M, s = get_sieve_region(f_x, B)
    table = get_dickman_table(20)

    alpha = alpha_score(f_x, primes[:350])
    E_score = get_Escore(f_x, [m1, -m0], alpha[0], B, M, round(B/math.sqrt(s)), table)
    E_score_2 = get_Epscore(f_x, [m1, -m0], alpha, B, M, round(B/math.sqrt(s)), table)
    L_norm = get_Lnorm(poly_prod(f_x, f_x), s, B)

    log.write_log(LOG_PATH, "f(x) = "+str(f_x)+"    g(x) = "+str([m1, -m0]))
    log.write_log(LOG_PATH, "alpha = "+str(alpha[0]))
    log.write_log(LOG_PATH, "E-score = "+str(E_score))
    log.write_log(LOG_PATH, "E-score 2 = "+str(E_score_2))
    log.write_log(LOG_PATH, "L²-norm = "+str(L_norm))
    log.write_log(LOG_PATH, "skew = "+str(s)+"\n")
    
    return f_x, m0, M
    
def get_sieve_region(f, B):
    F = poly_prod(f, f)
    d = len(f)-1
    ratios = [math.log(1e-7+abs(f[i+1]/(1+abs(f[i])))) for i in range(d)]
    s = 2*int(math.exp(sum(ratios)/d)) # skew factor
    k = 1 # shift
    best_norm = None

    while k > 0:
        updated = False

        if s-k > 0:
            norm = get_Lnorm(F, s-k, B)
            if best_norm == None or norm < best_norm:
                s = s-k
                k <<= 1
                best_norm = norm
                updated = True

        if s+k < B:
            norm = get_Lnorm(F, s+k, B)
            if best_norm == None or norm < best_norm:
                s = s+k
                k <<= 1
                best_norm = norm
                updated = True

        if not updated: k >>= 1

    x = round(B*math.sqrt(s))

    return x, s
    
def get_Lnorm(F, s, B):
    sqrt = math.sqrt(s)
    n = len(F)
    
    base_X, base_Y = B*sqrt/2, 1+B/sqrt
    current_X, current_Y = pow(base_X, n), base_Y
    base_X, base_Y = base_X*base_X, base_Y*base_Y

    res = 0.0

    for i in range(0, n, 2):
        res += (2*current_X)*F[i]*(current_Y-1)/((i+1)*(n-i))
        current_X /= base_X
        current_Y *= base_Y

    return math.log(abs(res))/2

def get_Escore(f, g, alpha, B, x_limit, y_limit, table):
    n = 32
    log_B = math.log(B)
    
    x, w = leggauss(n)  # nodes & weights on [-1, 1]
    a, b = 0, math.pi
    xm = 0.5 * (b - a) * x + 0.5 * (b + a)
    wm = 0.5 * (b - a) * w

    res = 0
    for i in range(n):
        angle = xm[i]
        X = x_limit*math.cos(angle)
        Y = (y_limit-1)*math.sin(angle)+1
        res1, table = dickman((math.log(abs(eval_F(X, Y, f, len(f)-1)))-alpha)/log_B, table)
        res2, table = dickman((math.log(abs(eval_F(X, Y, g, len(g)-1))))/log_B, table)
        res += wm[i]*res1*res2

    return res

def get_Epscore(f, g, alpha, B, x_limit, y_limit, table):
    k,l,c = alpha[1],alpha[2],alpha[3]

    n = 16
    log_B = math.log(B)

    # Gauss–Legendre for mu in [c, c + 5c]
    mu_nodes, mu_weights = leggauss(n)
    a_mu, b_mu = c, 6*c
    mu_x = 0.5*(b_mu-a_mu)*mu_nodes + 0.5*(b_mu+a_mu)
    mu_w = 0.5*(b_mu-a_mu)*mu_weights

    # Gauss–Legendre for theta in [0, π]
    th_nodes, th_weights = leggauss(n)
    a_th, b_th = 0, math.pi
    th_x = 0.5*(b_th-a_th)*th_nodes + 0.5*(b_th+a_th)
    th_w = 0.5*(b_th-a_th)*th_weights

    res = 0
    for i in range(n):
        Y = non_central(k, l, mu_x[i])
        for j in range(n):
            angle = th_x[j]
            X_eval = x_limit*math.cos(angle)
            Y_eval = (y_limit-1)*math.sin(angle)+1
            res1, table = dickman((math.log(abs(eval_F(X_eval, Y_eval, f, len(f)-1)))-(mu_x[i]-c))/log_B, table)
            res2, table = dickman((math.log(abs(eval_F(X_eval, Y_eval, g, len(g)-1))))/log_B, table)
            res += mu_w[i]*th_w[j]*res1*res2*Y

    return res

def alpha_score(f, primes):
    E, F = 0, 0
    evals = [0]*primes[-1]
    tmp = [0]*len(f)
    for j in range(len(f)): tmp[j] = evaluate(f, j)
    for q in range(1, len(f)):
        for k in range(len(f)-1, q-1, -1): tmp[k] -= tmp[k-1]

    eval = tmp[0]
    for k in range(primes[-1]):
        evals[k] = eval
        eval += tmp[1]
        for q in range(1, len(f)-1): tmp[q] += tmp[q+1]

    f_prime = get_derivative(f)
    upto = len(f)+10
    baseline_term = 0
    for p in primes:
        log_p = math.log(p)
        baseline_term += log_p/(p-1)
        Ep,Fp = 0,0
        ramified_roots = []
        ramified = False
        for r in range(p):
            if not evals[r]%p:
                if not eval_mod(f_prime, r, p):
                    ramified = True
                    ramified_roots.append(r)
                    Ep += 1/p
                    Fp += 1/p
                else:
                    Ep += 1/(p-1)
                    Fp += (p+1)/((p-1)**2)

        if ramified:
            tmp2 = p
            for i in range(2, upto):
                new = []
                for r in ramified_roots:
                    if not eval_mod(f, r, tmp2*p):
                        for k in range(p): new.append(r+k*tmp2)
                        Ep += 1/tmp2
                        Fp += (2*i-1)/tmp2
                tmp2 *= p
                ramified_roots = new.copy()

            Ep += len(ramified_roots)/(pow(p, upto-2)*(p-1))
            Fp += len(ramified_roots)*(upto*upto+2*upto/(p-1)+(p+1)/((p-1)**2))/tmp2

        if not f[0]%p:
            frev = f[::-1]
            fdrev = get_derivative(frev)
            if eval_mod(fdrev,0,p):
                Ep += 1/(p-1)
                Fp += (p+1)/((p-1)**2)
            else:
                Ep += 1/p
                Fp += 1/p
                ramified_roots = [0]
                tmp2 = p
                for i in range(2, upto):
                    new = []
                    for r in ramified_roots:
                        if not eval_mod(frev, r, tmp2*p):
                            for k in range(p): new.append(r+k*tmp2)
                            Ep += 1/tmp2
                            Fp += (2*i-1)/tmp2
                    tmp2 *= p
                    ramified_roots = new.copy()

                Ep += len(ramified_roots)/(pow(p, upto-2)*(p-1))
                Fp += len(ramified_roots)*(upto*upto+2*upto/(p-1)+(p+1)/((p-1)**2))/tmp2
                
        tmpE, tmpF = p*Ep/(p+1), p*Fp/(p+1)
        E += tmpE*log_p
        F += (tmpF-tmpE*tmpE)*log_p*log_p
        
    k = 2*E-F/2
    lbd = F/2-E
    return E-baseline_term, k, lbd, baseline_term
    
def local_opt(f, g, B):
    m0, m1 = -g[1], g[0]

    d = len(f)-1
    if d <= 5:
        poly = []
        for i in range(3): poly.append(f[i]*binom(2-i, d-i))
        poly = poly_prod(poly, poly)
        poly2 = get_derivative(poly)
        zeros = get_complex_roots(poly2)
        min,mink = f[2]**2,0
        for r in zeros:
            if r.imag == 0:
                if evaluate(poly, math.ceil(r.real)) < min:
                    min, mink = evaluate(poly, math.ceil(r.real)), math.ceil(r.real)
                if evaluate(poly, math.floor(r.real)) < min:
                    min, mink = evaluate(poly, math.floor(r.real)), math.floor(r.real)
    else:
        poly = []
        for i in range(4): poly.append(f[i]*binom(3-i, d-i))
        zero = get_complex_roots(poly)
        min, mink = abs(f[3]), 0
        for r in zero:
            if r.imag == 0:
                if abs(evaluate(poly, math.ceil(r.real))) < min:
                    min, mink = abs(evaluate(poly, math.ceil(r.real))), math.ceil(r.real)
                if abs(evaluate(poly, math.floor(r.real))) < min:
                    min, mink = abs(evaluate(poly, math.floor(r.real))), math.floor(r.real)

    f = shift(f, mink)
    m0 -= mink*m1

    k, u, v, iteration = 1, 1, 1, 0
    _, s = get_sieve_region(f, B)
    min_n = get_Lnorm(poly_prod(f, f), s, B)
    while iteration < 500 and (k > 0 or u > 0 or v > 0):
        F = shift(f, -k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        flag = False
        if tmp < min_n: f, min_n, k, flag, m0 = F.copy(), tmp, k<<1, True, m0+k*m1

        F = shift(f, k)
        tmp = get_Lnorm(poly_prod(F, F), s, B)
        if tmp < min_n: f, min_n, k, flag, m0 = F.copy(), tmp, k<<1, True, m0-k*m1
        if flag and not u: u = 1
        if flag and not v: v = 1
        elif not flag: k >>= 1

        if u:
            flag = False
            F = f.copy()
            F[-3] += u*m1
            F[-2] -= u*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, u, flag = F.copy(), tmp, u<<1, True

            F = f.copy()
            F[-3] -= u*m1
            F[-2] += u*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, u, flag = F.copy(), tmp, u<<1, True
            if flag and not k: k = 1
            if flag and not v: v = 1
            elif not flag: u >>= 1

        if v:
            flag = False
            F = f.copy()
            F[-2] += v*m1
            F[-1] -= v*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, v, flag = F.copy(), tmp, v<<1, True

            F = f.copy()
            F[-2] -= v*m1
            F[-1] += v*m0
            tmp = get_Lnorm(poly_prod(F, F), s, B)
            if tmp < min_n: f, min_n, v, flag = F.copy(), tmp, v<<1, True
            if flag and not k: k = 1
            if flag and not u: u = 1
            elif not flag: v >>= 1

        iteration += 1
        _, s = get_sieve_region(f, B)

    return f, m0

def prime_combinations_with_indices(Q, l, B):
    n = len(Q)
    indices = [0] * l  # reuse buffer

    def backtrack(start, depth, product):
        if depth == l:
            yield tuple(indices)
            return
        for i in range(start, n - (l - depth) + 1):
            p = Q[i]
            if product * p >= B:
                break  # sorted Q means all further i will be too large
            indices[depth] = i
            yield from backtrack(i + 1, depth + 1, product * p)

    yield from backtrack(0, 0, 1)

def compute_e(m0, root_used, NB_ROOTS, prod, a_d, n, d):
    e = []
    tmp_m = m_mu(m0, root_used, [0]*NB_ROOTS, NB_ROOTS)
    base = poly(tmp_m, prod, a_d, n, d)[1]%prod
    for i in range(NB_ROOTS):
        line = [0]*d
        for j in range(d):
            if not i:
                tmp_m = m_mu(m0, root_used, [j]+[0]*(NB_ROOTS-1), NB_ROOTS)
                line[j] = poly(tmp_m, prod, a_d, n, d)[1]%prod
            elif j:
                tmp_m = m_mu(m0, root_used, [0]*i+[j]+[0]*(NB_ROOTS-i-1), NB_ROOTS)
                line[j] = (poly(tmp_m, prod, a_d, n, d)[1]-base)%prod
        e.append(line)

    return e

def compute_f(n, a_d, m0, d, prod, root_used, NB_ROOTS, e):
    f0 = (n-a_d*pow(m0, d))/(prod*prod*pow(m0, d-1))
    f = []

    for i in range(NB_ROOTS):
        line = [0]*d
        for j in range(d):
            line[j] = -(a_d*d*root_used[i][j]/pow(prod, 2)+e[i][j]/prod)
        f.append(line)

    return f, f0

def create_first_array(NB_ROOTS, f0, f, d):
    vec = [0]*(NB_ROOTS>>1)
    array1 = []

    while vec[-1] < d:
        U = (f0 + sum([f[j][vec[j]] for j in range(NB_ROOTS>>1)]))%1

        if not len(array1) or U > array1[-1][0]: array1.append([U, [i for i in vec]])
        else:
            tmp_a = 0
            tmp_b = len(array1)-1
            tmp = (tmp_a+tmp_b)>>1
            while tmp_a <= tmp_b:
                if array1[tmp][0] > U: tmp_b = tmp-1
                else: tmp_a = tmp+1
                tmp = (tmp_a+tmp_b)>>1
            array1.insert(tmp_a, [U, [i for i in vec]])
        vec[0] += 1
        for j in range(len(vec)-1):
            if vec[j] == d:
                vec[j] = 0
                vec[j+1] += 1
            else: break

    return array1

def create_second_array(NB_ROOTS, len_vec, d, f):
    vect = [0]*(NB_ROOTS-len_vec)
    array2 = []
    while vect[-1] < d:
        U = -sum([f[len_vec+j][vect[j]] for j in range(len(vect))])%1

        if not len(array2) or U > array2[-1][0]: array2.append([U, [i for i in vect]])
        else:
            tmp_a = 0
            tmp_b = len(array2)-1
            tmp = (tmp_a+tmp_b)>>1
            while tmp_a <= tmp_b:
                if array2[tmp][0] > U: tmp_b = tmp-1
                else: tmp_a = tmp+1
                tmp = (tmp_a+tmp_b)>>1
            array2.insert(tmp_a, [U, [i for i in vect]])
        vect[0] += 1
        for j in range(len(vect)-1):
            if vect[j] == d:
                vect[j] = 0
                vect[j+1] += 1
            else: break

    return array2

def m_mu(m0,roots,vec,l):
    return m0 + sum([roots[i][vec[i]] for i in range(l)])

def poly(m0, m1, a_d, n, d):
    c = [a_d]
    r = [n]
    for i in range(d-1, -1, -1):
        r.append((r[-1]-c[-1]*pow(m0, i+1))//m1)
        delta = -r[-1]*m1*invmod(m1,pow(m0, i))%(m1*pow(m0, i))
        c.append((r[-1]+delta)//pow(m0, i))
    for i in range(1, len(c)):
        if c[i] > m0//2:
            c[i] -= m0
            c[i-1] += m1
        elif c[i] < -m0//2:
            c[i] += m0
            c[i-1] -= m1
    return c

def root_sieve(f, g, primes, U):
    array = [0]*(U<<1)
    for p in primes:
        k = 1
        P = p
        while P < primes[-1]:
            for x in range(P):
                eval1, eval2 = eval_mod(f, x, P), eval_mod(g, x, P)
                for u in range(P):
                    if (eval1+u*eval2)%P == 0:
                        for z in range((u+U)%P, len(array), P): array[z] += math.log(p)/(pow(p, k-1)*(p+1))
            k += 1
            P *= p

    max,maxi = array[0],0
    for i in range(1, len(array)):
        if array[i] > max: max, maxi = array[i], i
    u = maxi-U
    print(u)

    f[-1] += u*g[-1]
    f[-2] += u*g[-2]

    return f

def root_sieve2(f, g, primes, U, V):
    array = [[0 for _ in range(V<<1)] for __ in range(U<<1)]
    for p in primes:
        k = 1
        P = p
        while P < primes[-1]:
            for x in range(P):
                eval1, eval2 = eval_mod(f, x, P), eval_mod(g, x, P)
                for u in range(P):
                    if eval2%p:
                        v = -(u*x+eval1*invmod(eval2, P))%P
                        for i in range((u+U)%P, len(array), P):
                            for j in range((v+V)%P, len(array[0]), P): array[i][j] += math.log(p)/(pow(p, k-1)*(p+1))
            k += 1
            P *= p

    max, maxu, maxv = array[0][0], 0, 0
    for i in range(U<<1):
        for j in range(V<<1):
            if array[i][j] > max: max, maxu, maxv = array[i][j], i, j

    u, v = maxu-U, maxv-V
    print(u, v)

    f[-1] += v*g[-1]
    f[-2] += v*g[-2]+u*g[-1]
    f[-3] += u*g[-2]

    return f