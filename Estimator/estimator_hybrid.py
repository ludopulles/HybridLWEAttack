from functools import partial
from sage.functions.log import log
from sage.functions.other import sqrt, binomial
from sage.misc.all import cached_function
from sage.rings.all import RR, ZZ
from sage.rings.infinity import PlusInfinity
from sage.symbolic.all import pi
import logging

# from estimator import *
from estimator import BKZ, SDis, _dual, lattice_reduction_cost, \
    lattice_reduction_opt_m, reduction_default_cost, stddevf
import pprint

oo = PlusInfinity()


@cached_function
def B_leq_h(n, h):
    """
    Returns the number of sparse ternary secrets of length n with norm <= h.
    """
    result = 0
    for i in range(h + 1):
        result += binomial(n, i) * 2**i
    return result


def hyb_estimate(f, n, alpha, q, bound, HW, mem_bound=2**80, **kwds):
    '''
    * Input 'bound' means,
    the lattice reduction will finds short vectors of (log) size < 'log(q) - bound'.

    * Our Strategy:
    take s = s1 + s2
    If HW(s1) <= l, HW(s2) <= k - l, then HW(s) <= k

    1. Store S = {(s1, b - As1) : s1} first
    sizeS = 2 ** (dim / 2)
    cost_post = tau * sizeS

    2. Compute As2 and check noisy collision.

    :param HW: fixed norm ("hamming weight") of the (sparse) ternary secret.
    :param mem_bound: Bound on the memory capacity for MITM (default: 2**80)
    '''

    # too small a step size leads to an early abort, too large a step
    # size means stepping over target
    step_size = int(n / 8)

    C = binomial(n, HW)
    B = (2 + 1/sqrt(2 * pi)) * alpha * q * q / (sqrt(2) * (2 ** bound))

    best_k = None
    k = int(n / 2 - step_size)
    # Optimize over k:
    while step_size > 0:
        best_h1 = None
        h1 = 5
        # Optimize over h1:
        while h1 <= min(k, HW):
            # <constants>
            # - cost_lat      : cost of lattice reduction (dual lattice attack)
            # - cost_post     : cost of postprocess (want to compute)
            # - cost_mem = N  : size of MITM table T

            h2 = min(HW - h1, k - h1, h1 + 4)
            current = f(n - k, alpha, q, bound, k, h1, h2, HW, **kwds)
            cost_lat = current["red"]
            tau = current["repeat"]

            # Dimension-Error Trade-Off creates a new LWE sample with `k` many
            # secret coordinates and `tau` many LWE samples.
            cost_mem = tau * B_leq_h(k, h1)
            if cost_mem > mem_bound:
                # The instantation of the attack uses too much memory!
                break

            cost_tableConstruct = B_leq_h(k, h1) * (tau * k + tau)
            cost_query = B_leq_h(k, h2) * 2**(4 * tau * B / q)
            # Numbers acquired from Table 2 using `m = tau` and `n = k`.

            probability = 0
            for i in range(0, h1 + h2 + 1):
                probability += binomial(n - k, HW - i) * binomial(k, i)
            probability = probability / C

            cost_post = cost_tableConstruct + cost_query
            total_cost = (cost_post + cost_lat) / probability

            current["rop"] = total_cost
            current["post"] = cost_post
            current["prob_inv"] = int(1/probability)
            current["k"] = k
            current["h1"] = h1
            current["h2"] = h2
            current["mem"] = cost_mem
            current = current.reorder(["rop"])

            logging.getLogger("guess").debug(current)

            if int(cost_query) > int(cost_lat):
                best_h1 = current
                best_h1["rop"] = oo
                #print('mem = %5.1f, cost_tableConstruct = %5.1f, cost_query = %5.1f, cost_lat = %5.1f' % (log(mem,2), log(cost_tableConstruct,2),log(cost_query,2), log(cost_lat,2)))
                break

            if best_h1 is None or current["rop"] < best_h1["rop"]:
                best_h1 = current
                # print('k = %4d, h1 = %4d, h2 = %4d  -->  rop = %5.1f' % (best_h1["k"], best_h1["h1"], best_h1["h2"], log(best_h1["rop"], 2)))
            h1 += 1

        if best_h1 is not None and (best_k is None or best_h1["rop"] < best_k["rop"]):
            best_k = best_h1
            k += step_size
            #print('k = %4d  -->  rop = %5.1f, h1 = %4d, h2 = %4d' % (best_k["k"], log(best_k["rop"], 2), best_k["h1"], best_k["h2"]))
        else:
            # we go back and half the step size
            step_size = step_size // 2
            k = best_k["k"] - step_size
            if k <= 0:
                k = step_size
    return best_k


def dual_scale_hyb(n, alpha, q, bound, k, h1, h2, HW,
                   m=oo, reduction_cost_model=reduction_default_cost, c=None,
                   use_lll=True):

    T = B_leq_h(k, h1)
    Q = B_leq_h(k, h2)
    B = (2 + 1/sqrt(2 * pi)) * alpha * q * q / (sqrt(2) * (2 ** bound))
    tau = log(T * Q, 2) / (1 - 4 * B / q)

    if c is None:
        c = RR(stddevf(alpha * q) * sqrt(2 * n - n) / sqrt(HW))

    best = _dual(n=n, alpha=alpha, q=q, m=m, reduction_cost_model=reduction_cost_model)
    delta_0 = best["delta_0"]

    if use_lll:
        scale = 2
    else:
        scale = 1

    while True:
        m_optimal = lattice_reduction_opt_m(n, q / c, delta_0)
        m_ = ZZ(min(m_optimal, m + n))
        v = scale * delta_0 ** m_ * (q / c) ** (n / m_)
        if log(v, 2) <= log(q, 2) - bound:
            break
        delta_0 = delta_0 - RR(0.00005)

    m_optimal = lattice_reduction_opt_m(n, q / c, delta_0)
    m_ = ZZ(min(m_optimal, m + n))
    ret = lattice_reduction_cost(reduction_cost_model, delta_0, m_, B=log(q, 2))

    if use_lll:
        lll = BKZ.LLL(m_, log(q, 2))
    else:
        lll = None
    ret = ret.repeat(times=tau, lll=lll)

    ret["m"] = m_
    ret["repeat"] = tau
    ret["d"] = m_
    ret["c"] = c

    ret = ret.reorder(["rop", "m"])
    logging.getLogger("repeat").debug(ret)
    return ret


def MITM_estimator(n, alpha, q, h=64, start_bound=10, max_bound=13,
                   step_size=1, mem_bound=2**80, verbose=True,
                   reduction_cost_model=reduction_default_cost):
    bound = float(start_bound)
    mitm_hyb = partial(hyb_estimate, dual_scale_hyb)
    best = None
    level = 0
    step = float(step_size)

    if verbose:
        print('Chosen Parameters : ')
        print('     n = %5d, log(q) = %5.1f, stddev = %5.2f, HW(s) = %4d' % (n, log(q, 2), RR(stddevf(alpha * q)), h))
        print()
        print('Start Estimation . . .')
        print()

    while level <= 2:
        res = mitm_hyb(n, alpha, q, bound, mem_bound=mem_bound, HW=h, reduction_cost_model=reduction_cost_model)

        if verbose:
            print("Optimizing with beta = %4d . . ." % res["beta"])

        if best is None or res["rop"] < best["rop"]:
            best = res
            while bound + step > max_bound and level <= 2:
                step = float(step/2)
                level += 1
            bound += step
        else:
            step = float(step/2)
            level += 1
            if level <= 2:
                bound -= step

    if verbose:
        print()
        print('== Bit-security : %5.1f with optimal parameters' % log(best["rop"], 2))
        print('     k = %5d, h1 = %2d, h2 = %2d, beta = %4d, mem = %5.1f' % (best["k"], best["h1"], best["h2"], best["beta"], log(best["mem"], 2)))
        print('             (For simplicity, we set k1 = k2 = k/2)')
        print()
        print("== Details")
        pprint.pprint(best)
    return best
