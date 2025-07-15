from math import comb, factorial
ROUNDING_FACTOR = 2**512
RF = RealField(prec=1024)


def round_to_rational(x):
    A = ZZ(round(x * ROUNDING_FACTOR))
    return QQ(A) / QQ(ROUNDING_FACTOR)


def average_variance(D):
    mu = 0.
    s = 0.

    for (v, p) in D.items():
        mu += v * p
        s += v * v * p

    s -= mu * mu
    return round_to_rational(mu), round_to_rational(s)


def binomial(x, y):
    try:
        binom = factorial(x) // factorial(y) // factorial(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    return binomial(2 * k, x + k) / 2.**(2 * k)


def build_centered_binomial_law(k):
    D = {}
    for i in range(-k, k + 1):
        D[i] = centered_binomial_pdf(k, i)
    return D


def build_kyber_compress_distribution(q, d):
    tmp_l = [int(int(round(int(int(round(vi * 2 ** d / q)) % 2**d) * (q/2**d))-vi) % q) for vi in range(q)]
    tmp_l = [vi if int(vi) < int(-vi % q) else -int(-vi%q) for vi in tmp_l]
    tmp_pdf = {}
    for vi in set(tmp_l):
        tmp_pdf[vi] = RR(tmp_l.count(vi) / q)
    tmp_mean = int(round(average_variance(tmp_pdf)[0]))
    new_pdf = {}
    for vi in tmp_pdf:
        new_pdf[vi-tmp_mean] = tmp_pdf[vi]
    return new_pdf


def build_uniform_distribution(start, end):
    tmp_D = {}
    for i in range(start, end+1):
        tmp_D[i] = 1/(end-start+1)
    return tmp_D


def build_ternary_distribution(total_num, one_num, _one_num):
    tmp_D = {-1: _one_num / total_num, 1: one_num / total_num}
    tmp_D[0] = 1 - tmp_D[1] - tmp_D[-1]
    return tmp_D


def calculate_convolution_pdf(n, D1, D2, bound0=1e-64, bound1=1e-64, bound2=1e-64, bound3=1e-10):
    new_D = {}
    for a1 in D1.keys():
        if D1[a1] < bound0:
            continue
        for b1 in D2.keys():
            if D2[b1] < bound0:
                continue
            tmp = a1*b1
            if tmp not in new_D:
                new_D[tmp] = D1[a1] * D2[b1]
            else:
                new_D[tmp] += D1[a1] * D2[b1]
    PR = PolynomialRing(RR, 'X')
    X = PR.gens()[0]
    rho_x = 0
    for j in new_D.keys():
        if new_D[j] < bound1:
            continue
        rho_x += new_D[j] * X**j
    rho_x_up = rho_x.numerator() ** n
    D_after_convolution = {}
    down_degree = rho_x.denominator().degree() * n
    for j in range(rho_x_up.degree()+1):
        tmp = rho_x_up[j]
        if tmp < bound2:
            continue
        D_after_convolution[j-down_degree] = tmp
    normalization = sum([RF(di) for di in list(D_after_convolution.values())])
    if abs(normalization-1) > bound3:
        raise ValueError('Not enough precision')
    else:
        for k1 in range(min(D_after_convolution.keys()), max(D_after_convolution.keys()) + 1):
            D_after_convolution[k1] /= normalization
            return D_after_convolution


def add_distribution(pdf1, pdf2, prob_bound=1e-64):
    pdf3 = {}
    for p1 in pdf1.keys():
        if pdf1[p1] < prob_bound:
            continue
        for p2 in pdf2.keys():
            if pdf2[p2] < prob_bound:
                continue
            tmp = p1 + p2
            if tmp not in pdf3:
                pdf3[tmp] = pdf1[p1] * pdf2[p2]
            else:
                pdf3[tmp] += pdf1[p1] * pdf2[p2]
    return pdf3


def mul_distribution(pdf1, pdf2, prob_bound=1e-64):
    pdf3 = {}
    for p1 in pdf1.keys():
        if pdf1[p1] < prob_bound:
            continue
        for p2 in pdf2.keys():
            if pdf2[p2] < prob_bound:
                continue
            tmp = p1 * p2
            if tmp not in pdf3:
                pdf3[tmp] = pdf1[p1] * pdf2[p2]
            else:
                pdf3[tmp] += pdf1[p1] * pdf2[p2]
    return pdf3


def compute_prob_on_correct_bits(_dim, err_prob, correct_bits):
    total_prob = 0
    for k in range(correct_bits+1):
        total_prob += RF(comb(_dim, k)) * RF(err_prob)**k * RF(1-err_prob)**(_dim-k)
    return total_prob


def my_erf(sigma, x, u):
    return 1/2 * (erf((x-u)/(sqrt(2)*sigma)) - erf((-x-u)/(sqrt(2)*sigma)))


def multiply_distribution_with_factor(_D, _fac):
    new_D = {}
    for vi in _D.keys():
        new_D[vi * _fac] = _D[vi]
    return new_D


def calculate_DRF(use_D, bound, dimension, correct_bits):
    if bound > max(use_D.keys()):
        print('No decryption error')
        pass
    single_dfr = 0
    for i in range(bound+1, max(use_D.keys())+1):
        if i in use_D:
            single_dfr += use_D[i]
    single_dfr *= 2
    total_dfr = 1-compute_prob_on_correct_bits(dimension, RF(single_dfr), correct_bits)
    return single_dfr, total_dfr