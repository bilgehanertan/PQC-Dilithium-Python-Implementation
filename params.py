# -*- coding: utf-8 -*-



from general import DILITHIUM_MODE

SEEDBYTES = 32
CRHBYTES = 48
N = 256
Q = 8380417
D = 13
ROOT_OF_UNITY = 1753

if (DILITHIUM_MODE==2):
  K = 4
  L = 4
  ETA = 2
  TAU = 39
  BETA = 78
  GAMMA1 = (1 << 17)
  GAMMA2 = ((Q-1)/88)
  OMEGA = 80

elif (DILITHIUM_MODE == 3):
    K = 6
    L = 5
    ETA = 4
    TAU = 49
    BETA = 196
    GAMMA1 = (1 << 19)
    GAMMA2 = ((Q-1)/32)
    OMEGA = 55


elif (DILITHIUM_MODE==5):
    K = 8
    L = 7
    ETA = 2
    TAU = 60
    BETA = 120
    GAMMA1 = (1 << 19)
    GAMMA2 = ((Q-1)/32)
    OMEGA = 75
    
POLYT1_PACKEDBYTES = 320
POLYT0_PACKEDBYTES = 416
POLYVECH_PACKEDBYTES = (OMEGA + K)

if (GAMMA1 == (1<<17)):
    POLYZ_PACKEDBYTES =  576
    
elif (GAMMA1 == (1 << 19)):
    POLYZ_PACKEDBYTES  = 640
    
    
if (GAMMA2== (Q-1)/88):
    POLYW1_PACKEDBYTES = 192

elif (GAMMA2== (Q-1)/32):
    POLYW1_PACKEDBYTES = 128
    

if (ETA==2):
    POLYETA_PACKEDBYTES  = 96

elif (ETA==4):
    POLYETA_PACKEDBYTES = 128
    
    
CRYPTO_PUBLICKEYBYTES = (SEEDBYTES + K*POLYT1_PACKEDBYTES)   

CRYPTO_SECRETKEYBYTES = (2*SEEDBYTES + CRHBYTES + L*POLYETA_PACKEDBYTES  + K*POLYETA_PACKEDBYTES + K*POLYT0_PACKEDBYTES)
                               
CRYPTO_BYTES = (SEEDBYTES + L*POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    