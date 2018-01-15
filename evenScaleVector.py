# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 09:31:51 2018

@author: cyinbox@gmail.com

Function to even scale numerical vector Tn (length N) to Tm (length M)
Inputs: numerical series: Tn, length to be scaled: M. Note: M can be larger or smaller than N.
Output: even scaled numerical vector of length M

Citation:
Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. 
Journal of theoretical biology, 382, 99-110.
"""

import math

def evenScaleVector(Tn,M):
 Tm = [0]*M
 Tm[0]=Tn[0]
 N = len(Tn)
 
 for k in range(1,M):
  Q = (k+1)*N/M
  R = math.floor((k+1)*N/M)
  print(k,Q,R)
  
  if Q==R: # When Q is an integer
    Tm[k]=Tn[R-1]
    print(Tm[k])
  else:      #When Q is not an integer
    if R==0: #when R==0, R-1 is the last element, need to set it as 1
     R=1
     Tm[1] = Tn[0] + (Q-R)*(Tn[1]-Tn[0])
    else:
     Tm[k] = Tn[R-1] + (Q-R)*(Tn[R]-Tn[R-1])
    
    print(Tm[k],Tn[R],Tn[R-1])
 return Tm

#Test
'''
Tn=[2,4,6,8,10,22,28,-10]
#Tn=[1,9,10,20,2,3,5,9,13,17,12]
M=20
Tm=evenScaleVector(Tn,M)
print(Tm)
'''


