# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%% 1. IMPORT LIBRARIES
# import seaborn as sns
# from matplotlib.lines import Line2D
import time 
import matplotlib.pyplot as plt
import random
# from tqdm import tqdm
import numpy as np 
# import scipy.stats as stats
# import sys
import pandas as pd
import math
#%% 2. RISK VECTOR REALIZATION
'''
LIST OF FUNCTIONS 
FOR GENERATING RISK VECTOR
'''
def coin_flip():
   randnum = round(random.uniform(0,100),2)
   if 0.00 <= randnum <= 16.34:
     return 0.09
   elif randnum <= 33.45:
     return 0.48
   elif randnum <= 48.89:
     return 1.07
   elif randnum <= 64.91:
     return 0.81
   elif randnum <= 80.78:
     return 0.42
   else:
     return 0.18
 

def get_risk_vector():
  '''
  fill risk_vec
  using coin flips
  sort risk_vec
  '''
  for i in range(0,pop):
    risk_vec[i]=round(coin_flip()/100,4)
    
#%% 4.a NON DILUTION CASE PERFORMANCE MEASURES
'''
LIST OF 
PERFORMANCE MEASURES
'''

def Exp_FalseNeg_ND(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return (1-const_sens)*risk_vec[Node_1]
  else:
    wt = 0
    for ele in range(Node_1,Node_2):
      wt = wt + risk_vec[ele]
    return (1-const_sens**2)*wt  
    
def Exp_FalsePos_ND(Node_1,Node_2):
  if Node_2==Node_1+1: 
    return (1-spec)*(1-risk_vec[Node_1])
  else:
    sum = 0 
    for ele in range(Node_1,Node_2):
      sum = sum + (1-risk_vec[ele])
    wt = 1
    for ele in range(Node_1,Node_2):
      wt = wt*(1-risk_vec[ele]) 
    return ((1-spec)*const_sens*sum)-((Node_2-Node_1)*(1-spec)*(const_sens+spec-1)*(wt))

def Exp_Test_Num_ND(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return 1
  else:
    wt = 1
    for ele in range(Node_1,Node_2):
      wt = wt*(1-risk_vec[ele])
    return 1+(Node_2-Node_1)*(const_sens-(const_sens+spec-1)*(wt))
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_ND(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_ND(Node_1,Node_2))+(lambda_2*Exp_FalsePos_ND(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_ND(Node_1,Node_2))
  return weight

#%% 4.b DILUTION CASE PERFORMANCE MEASURES

'''
LIST OF 
PERFORMANCE MEASURES
'''

def sens(n,k):
  '''
  Input- 
    n(int): number of people in group
    k(int): number of infected people

  Output- 
    returns sensitivity(float) from dilution equation
  '''
  return 1-spec*alpha**(k/(n**gamma))

def alpha_pr(Node_1,Node_2):
  return (1-alpha**(1/((Node_2-Node_1)**gamma)))


def Exp_FalseNeg_D(Node_1,Node_2): 
    if Node_2==Node_1+1: 
       return (1-sens(1,1))*risk_vec[Node_1]
    else:
        sum_1 = 0 
        for ele in range(Node_1,Node_2):
          sum_1 = sum_1 + (risk_vec[ele])
        sum_2 = 0
        for ele in range(Node_1,Node_2):
          sum_2 = sum_2 + ((risk_vec[ele])/(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])) 
        wt = 1
        for ele in range(Node_1,Node_2):
          wt = wt*(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])
        return  (1-sens(1,1))*sum_1+sens(1,1)*spec*alpha**(1/((Node_2-Node_1)**gamma))*wt*sum_2

def Exp_FalseNeg_D_Mod(Node_1,Node_2): 
    if Node_2==Node_1+1: 
       return (1-sens(1,1))*risk_vec[Node_1]
    else:
        sum_1 = sum(risk_vec[Node_1:Node_2])
        sum_2 = sum([((risk_vec[ele])/(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])) for ele in range(Node_1,Node_2)]) 
        wt = math.prod([(1-alpha_pr(Node_1,Node_2)*risk_vec[ele]) for ele in range(Node_1,Node_2)])
        return  (1-sens(1,1))*sum_1+sens(1,1)*spec*alpha**(1/((Node_2-Node_1)**gamma))*wt*sum_2


    
def Exp_FalsePos_D(Node_1,Node_2):
    if Node_2==Node_1+1: 
        return (1-spec)*(1-risk_vec[Node_1])
    else:
        sum_1 = 0 
        for ele in range(Node_1,Node_2):
          sum_1 = sum_1 + (1-risk_vec[ele])
        sum_2 = 0
        for ele in range(Node_1,Node_2):
          sum_2 = sum_2 + ((1-risk_vec[ele])/(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])) 
        wt = 1
        for ele in range(Node_1,Node_2):
          wt = wt*(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])
        return ((1-spec)*sum_1)-((1-spec)*spec*wt*sum_2)
    
    
def Exp_FalsePos_D_Mod(Node_1,Node_2):
    if Node_2==Node_1+1: 
        return (1-spec)*(1-risk_vec[Node_1])
    else:
        sum_1 = Node_2-Node_1-sum(risk_vec[Node_1:Node_2]) 
        sum_2 = sum([((1-risk_vec[ele])/(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])) for ele in range(Node_1,Node_2)])
        wt = math.prod([(1-alpha_pr(Node_1,Node_2)*risk_vec[ele]) for ele in range(Node_1,Node_2)])
        return ((1-spec)*sum_1)-((1-spec)*spec*wt*sum_2)    

def Exp_Test_Num_D(Node_1,Node_2): 
    if Node_2==Node_1+1: 
        return 1
    else:
        wt = 1
        for ele in range(Node_1,Node_2):
          wt = wt*(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])
        return 1+(Node_2-Node_1)*(1-spec*(wt))

def Exp_Test_Num_D_Mod(Node_1,Node_2): 
    if Node_2==Node_1+1: 
      return 1
    else:
        wt = math.prod([(1-alpha_pr(Node_1,Node_2)*risk_vec[ele]) for ele in range(Node_1,Node_2)])
        return 1+(Node_2-Node_1)*(1-spec*(wt))
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_D(Node_1,Node_2):
  # weight = (lambda_1*Exp_FalseNeg_D_Mod(Node_1,Node_2))+(lambda_2*Exp_FalsePos_D_Mod(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_D_Mod(Node_1,Node_2))
  return (lambda_1*Exp_FalseNeg_D_Mod(Node_1,Node_2))+(lambda_2*Exp_FalsePos_D_Mod(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_D_Mod(Node_1,Node_2))

def weight_D_old(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_D(Node_1,Node_2))+(lambda_2*Exp_FalsePos_D(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_D(Node_1,Node_2))
  return weight
#%%
lambda_1=0.3
lambda_2=0.3
pop=100
risk_vec=[0]*pop
get_risk_vector()

#Parameters
spec = 0.99
alpha= 0.0514
gamma= 0.2416
const_sens = sens(1,1)

#Confidence Interval
CI = 0.05

#Simulations
simulations=50

#%%
a=list(range(0,100000))
b=np.empty(10)
for sim in range(0,10):
    start_time=time.time()
    summ=0
    for i in range(len(a)):
        summ+=a[i]
    b[sim]=(time.time() - start_time)
#%%
b=np.empty(10)
for sim in range(0,10):
    start_time=time.time()
    summ=sum(a)
    b[sim]=(time.time() - start_time)
#%%
a=[0.999929]*1000
b=np.empty(10)
for sim in range(0,10):
    start_time=time.time()
    prod=1
    for i in range(len(a)):
        prod=prod*a[i]
    b[sim]=(time.time() - start_time)
#%%
b=np.empty(10)
for sim in range(0,10):
    start_time=time.time()
    log_a=[np.log(a[i]) for i in range(len(a))]
    summ=np.exp(sum(log_a))
    b[sim]=(time.time() - start_time)    

#%%
b=np.empty(10)
for sim in range(0,10):
    start_time=time.time()
    mpprod=math.prod([(a[i]) for i in range(len(a))])
    b[sim]=(time.time() - start_time)    
#%%
# a=np.empty(100)
b=np.empty(100)
for sim in range(0,100):
    start_time=time.time()
    for i in range(len(a)):
        a[i]=weight_D(5,30)
    b[sim]=(time.time() - start_time) 
#%%
a=np.empty(100)
b=np.empty(100)
for sim in range(0,100):
    start_time=time.time()
    for i in range(len(a)):
        a[i]=weight_D_old(5,30)
    b[sim]=(time.time() - start_time) 
#%%
# get_risk_vector()
b=np.empty(50)
for sim in range(0,50):
    pop=100
    start_time=time.time()
    # g_d = Graph(pop+1) 
    # wt_D =np.array([[None if i>=j else weight_D(i,j) for i in range(pop+1)] for j in range(pop+1)])
    # wt_D =[[weight_D(i,j) for i in range(0,pop)] for j in range(i+1,pop+1)]
    wt_D =np.array([[None for i in range(pop+1)]for i in range(pop)])
    for i in range(0,pop):
        for j in range(i+1,pop+1):
              wt_D[i][j]=Exp_FalseNeg_D(i,j)
    # for i in range(0,pop):
    #     for j in range(i+1,pop+1):
    #         g_d.addEdge(i,j,wt_D[i][j])
    b[sim]=(time.time() - start_time) 
#%%# get_risk_vector()
# alpha_p=np.array([(1-alpha**(1/((i)**gamma))) for i in range(1,pop+1)])
b=np.empty(50)
for sim in range(0,50):
    pop=100
    start_time=time.time()
    # g_d = Graph(pop+1) 
    # wt_D =np.array([[None if i>=j else weight_D(i,j) for i in range(pop+1)] for j in range(pop+1)])
    # wt_D =[[weight_D(i,j) for i in range(0,pop)] for j in range(i+1,pop+1)]
    wt_D2 =np.array([[Exp_FalseNeg_D_Mod(j,i) if i>j else None for i in range(pop+1)]for j in range(pop)])
    # for i in range(0,pop):
    #     for j in range(i+1,pop+1):
    #         g_d.addEdge(i,j,wt_D[i][j])
    b[sim]=(time.time() - start_time) 