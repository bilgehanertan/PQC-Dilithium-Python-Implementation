# -*- coding: utf-8 -*-
DILITHIUM_MODE = 2
from params import *
from poly import *
import numpy as np
class polyvecLstruct:
    def __init__(self):
        self.poly.vec[K]# Poly Lib Implement vector poly vec[L]
        
class polyvecKstruct:
    def __init__(self):
        self.poly.vec[K]# Poly Lib Implement vector poly vec[L]
                
"""
/*************************************************
* Name:        expand_mat
*
* Description: Implementation of ExpandA. Generates matrix A with uniformly
*              random coefficients a_{i,j} by performing rejection
*              sampling on the output stream of SHAKE128(rho|j|i)
*              or AES256CTR(rho,j|i).
*
* Arguments:   - polyvecl mat[K]: output matrix
*              - const uint8_t rho[]: byte array containing seed rho
**************************************************/
"""        
def polyvec_matrix_expand(polyvecL, rho_SEEDBYTES_ELEMENT):
                                                              
    for i in range(K):
        for j in range(L):
            print("a")
            poly_uniform(mat[i].vec[j], rho, (i << 8) + j); # to be implemented
        
    
def polyvec_matrix_pointwise_montgomery(polyvecK, polyvecL, polyvecL2):
    for i in range(K):
        polyvec_matrix_pointwise_montgomery(polyvecK.poly.vec[i], polyvecL.poly.vec[i], polyvecL2)
        
    """ /**************************************************************/
/************ Vectors of polynomials of length L **************/
/**************************************************************/"""
    
def polyvecl_uniform_eta(polyvecL,seed,nonce): # pass seed[SEEDBYTES]
    for i in range(L):
         poly_uniform_eta(polyvecL.poly.vec[i],seed[SEEDBYTES],L*nonce+i)

     

def polyvecl_uniform_gamma1(polyvec, seed, nonce):
    for i in range (L):
        poly_uniform_gamma1(polyVecL.poly.vec[i],seed[SEEDBYTES],L*nonce+i)
        

def polyvecl_reduce(polyvecL):
    for i in range (L):
        poly_reduce(polyvecL.poly.vec[i])
        
        
        """/*************************************************
* Name:        polyvecl_freeze
*
* Description: Reduce coefficients of polynomials in vector of length L
*              to standard representatives.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/ """


def polyvecl_freeze(polyvecL):
    for i in range (L):
        poly_freeze(polyvecL.poly.vec[i])
        
        

"""/*************************************************
* Name:        polyvecl_add
*
* Description: Add vectors of polynomials of length L.
*              No modular reduction is performed.
*
* Arguments:   - polyvecl *w: pointer to output vector
*              - const polyvecl *u: pointer to first summand
*              - const polyvecl *v: pointer to second summand
**************************************************/"""

def polyvecl_add(polyvecl_w, polyvecl_u, polyvecl_v):
    for i in range (L):
        poly_add(polyvecl_w.poly.vec[i],polyvecl_u.poly.vec[i],polyvecl_v.poly.vec[i]) # to be implemented



""" /*************************************************
* Name:        polyvecl_ntt
*
* Description: Forward NTT of all polynomials in vector of length L. Output
*              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/"""


def polyvecl_ntt(polyvecl_v):
    for i in range (L):
        poly_ntt(polyvecl_v.poly.vec[i])

def polyvecl_invntt_tomont(polyvecl_v):
    for i in range(L):
        poly_invntt_tomont(polyvecl_v.poly.vec[i])


def polyvecl_pointwise_poly_montgomery(polyvecl_r,poly_a,polyvecl_v):
    for i in range(L):
        poly_pointwise_montgomery(polyvecl_r.poly.vec[i],poly_a,polyvecl_v.poly.vec[i])
        
        

""" /*************************************************
* Name:        polyvecl_pointwise_acc_montgomery
*
* Description: Pointwise multiply vectors of polynomials of length L, multiply
*              resulting vector by 2^{-32} and add (accumulate) polynomials
*              in it. Input/output vectors are in NTT domain representation.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvecl *u: pointer to first input vector
*              - const polyvecl *v: pointer to second input vector
**************************************************/"""








































        