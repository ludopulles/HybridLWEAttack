# HybridLWEAttack

## Estimator

Our "estimator_hybrid.py" sage module estimates the bit-security of the hybrid MITM attack, given LWE parameters.
We uses functions in "estimator.py" made by Albrecht et al [[APS15]](https://eprint.iacr.org/2015/046)
and the high-level structure of our estimator "estimator_hybrid.py" is from "estimator.py".

### Note 
The current version is somewhat unstable for estimating small parameters such as n < 1000 and log(q) < 30.

### How to Run

We recommend to use our **MITM_estimator** function in "estimator_hybrid.py". This function estimates the cost of ring operations, given a bounded memory capacity 2^80.

EXAMPLE :

    > sage: from estimator_hybrid import *
    > sage: n = 8192; q = next_prime(2^125); alpha = 8/q; h = 64;
    > sage: MITM_estimator(n, alpha, q, h, reduction_cost_model=BKZ.sieve)
    
    Chosen Parameters :
             n =  8192, log(q) = 125.0, stddev =  3.19, HW(s) =   64
     
    Start Estimation . . .

    Optimizing with beta =  268 . . .
    Optimizing with beta =  276 . . .
    Optimizing with beta =  272 . . .
    Optimizing with beta =  244 . . .

    == Bit-security : 131.2 with optimal parameters
        k =  4606, h1 =  7, h2 = 11, beta =  268, mem =  72.9
                (For simplicity, we set k1 = k2 = k/2)

    == Details
             rop:  2^131.2
               m:   2^12.9
             red:  2^110.5
         delta_0: 1.005218
            beta:      268
          repeat: 1.031899
               d:   2^12.9
               c:   23.890
            post:  2^108.6
        prob_inv:   2^20.3
               k:   2^12.2
              h1:        7
              h2:       11
             mem:   2^72.9
             
## Attack Implementation

"mitm.py" is an attack code for Sage.
This is just for a proof-of-concept of the hybrid attack, and hence it is definitely not optimized.

### EXAMPLE :

    > sage: from Mitm import *
    > sage: n = 45; q = next_prime(2^14); h = 6; beta = 20; k = 30; tau = 30;
    > sage: A, c, u, s = gen_instance(n, q, h, m = tau * n)
    
    (A, c) is n-dim LWE samples (with secret s) / (A, u) is uniform samples
    
    > sage: A_k, c_k, u_k, B = dim_error_tradeoff(A, c, u, beta, h, k)
    
    (A_k, c_k) is k-dim LWE samples (with secret s[-k:]) / (A_k, u_k) is uniform samples. 
    
    > sage: Mitm_on_LWE(A_k, c_k, u_k, B, h)
    
    Table size = 34281

    ** Mitm on (A_k, c_k) ** 

    Number of noisy searches = 193
     - Input is LWE with secret
    (0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)

    ** Mitm on (A_k, u_k) **

    Number of noisy searches = 34281
     - Input is uniform
     
    > sage: s[-k:]
    
    (0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)




