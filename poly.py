# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:05:42 2020

@author: Ria
"""

DILITHIUM_MODE = 2
from params import *

import numpy as np
from reduce import * 
class poly:
    def __init__(self):
        self.coeffs = []# Poly Lib Implement vector poly vec[L]
        
        
""" /*************************************************
* Name:        poly_reduce
*
* Description: Inplace reduction of all coefficients of polynomial to
*              representative in [-6283009,6283007].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/"""

def poly_reduce(poly_a):
    for i in range (N):
        poly_a.coeffs[i] = reduce32(poly_a.coeffs[i])


""" 
/*************************************************
* Name:        poly_caddq
*
* Description: For all coefficients of in/out polynomial add Q if
*              coefficient is negative.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/"""

def poly_cadq(poly_a):
    for i in range(N):
        poly_a.poly.coeffs[i] = caddq(poly_a.poly.coeffs[i])
        
        
        
"""/*************************************************
* Name:        poly_freeze
*
* Description: Inplace reduction of all coefficients of polynomial to
*              standard representatives.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/        """


def poly_freeze(poly_a):
    for i in range (N):
        poly_a.coeffs[i] = freeze(poly_a.coeffs[i])
        
        
        
""" /*************************************************
* Name:        poly_add
*
* Description: Add polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
**************************************************/"""

def poly_add(poly_c, poly_a, poly_b):
    for i in range(N):
        poly_c.coeffs[i] = poly_a.coeffs[i] + poly_b.coeffs[i]
        
    
""" /*************************************************
* Name:        poly_sub
*
* Description: Subtract polynomials. No modular reduction is
*              performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtraced from first i""" 

def poly_sub(poly_c, poly_a, poly_b):
    for i in range(N):
        poly_c.coeffs[i] = poly_a.coeffs[i] - poly_b.coeffs[i]


""" /*************************************************
* Name:        poly_shiftl
*
* Description: Multiply polynomial by 2^D without modular reduction. Assumes
*              input coefficients to be less than 2^{31-D} in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/"""

def poly_shiftl(poly_a):
    for i in range (N):
        poly_a.coeffs[i] <<=D
        
        
"""/*************************************************
* Name:        poly_ntt
*
* Description: Inplace forward NTT. Coefficients can grow by
*              8*Q in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/"""

def poly_ntt(poly_a):
    ntt(poly_a.coeffs)
    
""" /*************************************************
* Name:        poly_invntt_tomont
*
* Description: Inplace inverse NTT and multiplication by 2^{32}.
*              Input coefficients need to be less than Q in absolute
*              value and output coefficients are again bounded by Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/"""

def poly_invntt_tomont(poly_a):
    invntt_tomont(poly_acoeffs)
    
"""/*************************************************
* Name:        poly_pointwise_montgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              by 2^{-32}.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/"""

def poly_pointwise_montgomery(poly_c, poly_a, poly_b):
    for i in range (N):
        poly_c.coeffs[i] = montgomery_reduce(poly_a.coeffs[i]*poly_b.coeffs[i])
        
        
"""/*************************************************
* Name:        poly_power2round
*
* Description: For all coefficients c of the input polynomial,
*              compute c0, c1 such that c mod Q = c1*2^D + c0
*              with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
*              standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/"""

def poly_power2round(poly_a1,poly_a0,poly_a):
    for i in range(N):
        poly_a1.coeffs[i] = power2round(poly_a0.coeffs[i], poly_a.coeffs[i])
        
""" /*************************************************
* Name:        poly_decompose
*
* Description: For all coefficients c of the input polynomial,
*              compute high and low bits c0, c1 such c mod Q = c1*ALPHA + c0
*              with -ALPHA/2 < c0 <= ALPHA/2 except c1 = (Q-1)/ALPHA where we
*              set c1 = 0 and -ALPHA/2 <= c0 = c mod Q - Q < 0.
*              Assumes coefficients to be standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/"""       
    
def poly_decompose(poly_a1,poly_a0,poly_a):
    for i in range(N):
        poly_a1.coeffs[i] = decompose(poly_a0.coeffs[i], poly_a.coeffs[i])
        
""" /*************************************************
* Name:        poly_make_hint
*
* Description: Compute hint polynomial. The coefficients of which indicate
*              whether the low bits of the corresponding coefficient of
*              the input polynomial overflow into the high bits.
*
* Arguments:   - poly *h: pointer to output hint polynomial
*              - const poly *a0: pointer to low part of input polynomial
*              - const poly *a1: pointer to high part of input polynomial
*
* Returns number of 1 bits.
**************************************************/"""
        
def poly_make_hint(poly_h,poly_a0, poly_a1):
    for i in range(N):
        poly_h.coeffs[i] = make_hint(poly_a0.coeffs[i], poly_a1.coeffs[i])
        s += poly_h.coeffs[i]
    return s

""" /*************************************************
* Name:        poly_use_hint
*
* Description: Use hint polynomial to correct the high bits of a polynomial.
*
* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
*              - const poly *a: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
**************************************************/"""

def poly_use_hint(poly_b,poly_a,poly_h):
    for i in range(N):
        poly_b.coeffs[i] = use_hint(poly_a.coeffs[i], poly_h.coeffs[i])
        
""" /*************************************************
* Name:        poly_chknorm
*
* Description: Check infinity norm of polynomial against given bound.
*              Assumes input coefficients were reduced by reduce32().
*
* Arguments:   - const poly *a: pointer to polynomial
*              - int32_t B: norm bound
*
* Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
**************************************************/"""
                                                   
def poly_chknorm(poly_a,B):
    
    if(B> (Q-1)/8):
        return 1
    
    for i in range (N):
        t = poly_a.coeffs[i] >> 31
        t = poly_a.coeffs[i] - (t & 2*poly_a.coeffs[i])
        
        if(t >=B):
            return 1
        
    
    return 0


""" /*************************************************
* Name:        rej_uniform
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/"""

def rej_uniform(a,leng, buf, buflen):
    
    ctr = pos = 0
    while(ctr < len and pos+3 <= buflen):
        pos+=1
        t = buf[pos]
        pos+=1
        t |= buf[pos] << 8
        pos+=1
        t |= buf[pos] << 16
        
        if (t < Q):
            ctr += 1
            a[ctr+1] = t
            
    return ctr


POLY_UNIFORM_NBLOCKS = ((768 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)




























































