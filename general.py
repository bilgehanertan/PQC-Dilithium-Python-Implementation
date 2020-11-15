DILITHIUM_MODE = 2

import numpy as np
from params import *
# -*- coding: utf-8 -*-





def Power2Round(r,d):

  a1 = (d+(1<<(D-1))-1)>>D
  r = d-(a1<<D)
  return a1


def makeHint(lowBits,highBits):
    
  if(lowBits<=GAMMA2 or a0>Q-GAMMA2 or (lowBits==Q-GAMMA2 and highBits==0)):
      
    return 0
  return 1


def use_hint(a, hint):
    
  a1 = decompose(a0,a)
  if(hint==0):
      
    return a1

  if (GAMMA2 == (Q-1)/32):
      
    if(a0>0):
        
      return (a1+1) and 15
  
    else:
        
      return (a1-1) and 15
  
  elif(GAMMA2==(Q-1)/88):
      
    if(a0>0):
        
      if(a1==43):
          
        return 0
    
      else:
          
        return a1+1
    
    else:
        
      if(a1==0):
          
        return 43
    
      else:
          
        return a1-1




