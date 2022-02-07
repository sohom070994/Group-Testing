# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%% 1. IMPORT LIBRARIES
# import seaborn as sns
# from matplotlib.lines import Line2D
# import time 
# import matplotlib.pyplot as plt
# import random
# from tqdm import tqdm
import numpy as np 
# import scipy.stats as stats
# import sys
import pandas as pd
# import math
#%% 2. DILUTION CASE PERFORMANCE MEASURES PER PATIENT FOR PREVALENCE RATE

'''
LIST OF 
PERFORMANCE MEASURES
'''

def sens(n,k):
  return 1-spec*alpha**(k/(n**gamma))

def alpha_pr(Node_1,Node_2):
  return (1-alpha**(1/((Node_2-Node_1)**gamma)))

def funcgs(gs,pr):
   return 1-(1-alpha**(1/(gs**gamma)))*pr 
    
    
def Exp_FalseNeg_DH(gs,pr): 
    if gs==1: 
        return (1-sens(1,1))*pr
    else:
        return (1-sens(1,1))*pr+sens(1,1)*spec*(alpha**(1/(gs**gamma)))*((funcgs(gs,pr)**(gs-1)))*pr

def Exp_FalsePos_DH(gs,pr):
    if gs==1: 
        return (1-spec)*(1-pr)
    else:
        return (1-spec)*(1-pr)*(1-spec*(funcgs(gs,pr)**(gs-1)))    

def Exp_Test_Num_DH(gs,pr): 
    if gs==1: 
        return 1
    else:
        return (1/gs)+(1-spec*(funcgs(gs,pr)**gs))
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_DH(gs,pr):
  return (lambda_1*Exp_FalseNeg_DH(gs,pr))+(lambda_2*Exp_FalsePos_DH(gs,pr))+((1-lambda_1-lambda_2)*Exp_Test_Num_DH(gs,pr))
#%% 3. ESTABLISH CONSTANT PARAMETERS
#population
pop = 100
prev= 0.49/100

#COVID Parameters
spec = 0.99
alpha= 0.0514
gamma= 0.2416
const_sens = sens(1,1)
#%% 4. SIMULATION


results = pd.DataFrame(columns= ['lambda1',
                                 'lambda2',
                                 'minOF',
                                 'minGS'])
for w1 in range(0,12,2):
   for w2 in range(0,12,2):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10
           
           OF_DH=[1000]*(pop)
           for sim in range(1,pop+1):
               OF_DH[sim-1]=weight_DH(sim, prev)
           minwt, mings = np.min(OF_DH), np.argmin(OF_DH)+1
          
           results = results.append({'lambda1':lambda_1,
                       'lambda2':lambda_2,
                       'minOF':minwt, 
                       'minGS':mings},ignore_index=True)