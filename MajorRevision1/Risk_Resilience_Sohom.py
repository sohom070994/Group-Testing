# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""

# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%% 
# 1. IMPORT LIBRARIES
# from matplotlib import pyplot
import random
from tqdm import tqdm
import numpy as np 
# import scipy.stats as stats
# import sys
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

#%% 
# 2. RISK VECTOR REALIZATION
'''
LIST OF FUNCTIONS 
FOR GENERATING RISK VECTOR
'''
def perturb_rv(rv,delta):
    a=[]
    for ele in range(0,len(rv)):
        a.append(random.uniform(rv[ele]*(1-delta/100),rv[ele]*(1+delta/100)))
        # a.append(random.uniform(rv[ele]*(1+delta/100),rv[ele]*(1+2*delta/100)))
    return a

# def get_risk_vector():
#     # COVID Risk Distribution
#     elements = [0.0198,0.0138,0.046,0.0395,0.0126]
#     probabilities = [0.025,0.3715,0.2685,0.0031,0.3319]
#     a=np.random.choice(elements, pop, p=probabilities)
#     return a


# HBV risk distribution
def get_risk_vector():
    elements = [0.0009,0.0048,0.0107,0.0081,0.0042,0.0018]
    probabilities = [0.16339801983321,0.171090634447457,0.15443199329378,0.160180315971968,0.158686396851769,0.192212639601815]
    a=np.random.choice(elements, pop, p=probabilities)
    return a
#%% 
# 3. Debugging Block
# risk_vec=get_risk_vector()
# b=perturb_rv(risk_vec,10)
# # print(b)
# # print(risk_vec)
# c=[(b[ele]-risk_vec[ele])/risk_vec[ele]*100 for ele in range(0,len(b))]
# print(c)
#%% 
# 4.a NON DILUTION CASE PERFORMANCE MEASURES
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

#%% 
# 4.b DILUTION CASE PERFORMANCE MEASURES

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

def Exp_Test_Num_D(Node_1,Node_2): 
    if Node_2==Node_1+1: 
      return 1
    else:
        wt = 1
        for ele in range(Node_1,Node_2):
          wt = wt*(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])
        return 1+(Node_2-Node_1)*(1-spec*(wt))
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_D(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_D(Node_1,Node_2))+(lambda_2*Exp_FalsePos_D(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_D(Node_1,Node_2))
  return weight
#%% 
# 4.C U-RDDS PERFORMANCE MEASURES
'''
LIST OF 
PERFORMANCE MEASURES
'''
def risk_sum(Node_1,Node_2):
    s = 0 
    for ele in range(Node_1,Node_2):
      s = s + (1-risk_vec[ele])
    return s  

def gen_sensBM(n,k):
  '''
  Input- 
    n(int): number of people in group
    k(int): number of infected people

  Output- 
    returns sensitivity(float) from dilution equation
  '''
  return 1*(k/n)**0.0507


def Exp_FalseNeg_UB(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return (1-const_sens)*risk_vec[Node_1]
  else:
    sum = risk_sum(Node_1,Node_2)
    return gen_sensBM(1, 1)*(1-gen_sensBM(Node_2-Node_1, sum))*sum
    
def Exp_FalsePos_UB(Node_1,Node_2):
  if Node_2==Node_1+1: 
    return (1-spec)*(1-risk_vec[Node_1])
  else:
    sum = risk_sum(Node_1,Node_2)
    return (1-spec)*(Node_2-Node_1-sum)*gen_sensBM(Node_2-Node_1, sum)

def Exp_Test_Num_UB(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return 1
  else:
    sum = risk_sum(Node_1,Node_2)
    return 1+(Node_2-Node_1)*gen_sensBM(Node_2-Node_1, sum)
    
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_UB(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_UB(Node_1,Node_2))+(lambda_2*Exp_FalsePos_UB(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_UB(Node_1,Node_2))
  return weight

#%% 
# 4.D POISSON-BINOMIAL PERFORMANCE MEASURES

def PBD(m):
    summ=0.0
    for n in range(0,N+1):
        prod=1.0
        for k in range(0,N):
            prod=prod*(p[k]*cmath.exp(j*2*pi*n/(N+1.0))+(1-p[k]))
        summ=summ+cmath.exp(-j*2*pi*n*m/(N+1.0))*prod            
    return summ/(N+1.0)

def Exp_FalseNeg_pbd(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return (1-const_sens)*risk_vec[Node_1]
  else:
    sum = risk_sum(Node_1,Node_2)
    return gen_sensBM(1, 1)*(1-gen_sensBM(Node_2-Node_1, sum))*sum
    
def Exp_FalsePos_pbd(Node_1,Node_2):
  if Node_2==Node_1+1: 
    return (1-spec)*(1-risk_vec[Node_1])
  else:
    sum = risk_sum(Node_1,Node_2)
    return (1-spec)*(Node_2-Node_1-sum)*gen_sensBM(Node_2-Node_1, sum)

def Exp_Test_Num_pbd(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return 1
  else:
    sum = risk_sum(Node_1,Node_2)
    return 1+(Node_2-Node_1)*gen_sensBM(Node_2-Node_1, sum)
    
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_pbd(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_pbd(Node_1,Node_2))+(lambda_2*Exp_FalsePos_pbd(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_pbd(Node_1,Node_2))
  return weight

#%% 
# 5. SORTING AND SHORTEST PATH
'''
TOPOLOGICAL SORT ALGORITHM
FOR SHORTEST PATH
'''
from collections import defaultdict

'''
DEFINING GRAPH DICTIONARY 
AND RELATED FUNCTIONS
'''
# Graph is represented using adjacency list. Every 
# node of adjacency list contains vertex number of 
# the vertex to which edge connects. It also contains 
# weight of the edge
class Graph:
  def __init__(self,vertices):

    # No. of vertices 
    self.V= vertices

    # dictionary containing adjacency List 
    self.graph = defaultdict(list)

  # function to add an edge to graph 
  def addEdge(self,u,v,w):
    self.graph[u].append((v,w))


  # A recursive function used by shortestPath function
  def topologicalSortUtil(self,v,visited,stack):

    # Mark the current node as visited.
    visited[v] = True


    # Recur for all the vertices adjacent to this vertex
    if v in self.graph.keys():
      for node,weight in self.graph[v]:
        if visited[node] == False:
          self.topologicalSortUtil(node,visited,stack)

    # Push current vertex to stack which stores topological sort
    stack.append(v)


  #The function to find shortest paths from given vertex. 
  #It uses recursive topologicalSortUtil() to get topological 
  #sorting of given graph.
  def shortestPath(self,s):

    # Mark all the vertices as not visited 
    visited = [False]*self.V
    stack=[]

    # Call the recursive helper function to store Topological 
    # Sort starting from source vertex 
    for i in range(self.V):
      if visited[i] == False:
        self.topologicalSortUtil(s,visited,stack)

    # Initialize distances to all vertices as infinite and 
    # distance to source as 0
    dist = [float("Inf")]*(self.V)
    dist[s] = 0
    part = [0]*(self.V)

    # print(dist)
    # print(stack)

    # Process vertices in topological order 
    while stack:

      # Get the next vertex from topological order
      i=stack.pop()

      # Update distances of all adjacent vertices
      for node,weight in self.graph[i]:
        if dist[node] > dist[i] + weight:
          dist[node] = dist[i] + weight
          part[node] = (node,i)

    # Print the distance from source to Node_100
    # for i in range(self.V):
    # print ("Dist - %f, Node - %d" %(dist[100],100)) #if dist[i] != float ("Inf") else "Inf" ,
    node_list=[]
    a=part[pop]
    while (a[1]!=0):
        node_list.append(a[1])
        a=part[a[1]]
    node_list.append(0)  
    
    # et=0
    # et+=Exp_Test_Num(node_list[0],pop)
    # efp=0
    # efp+=Exp_FalsePos(node_list[0],pop)
    # efn=0
    # efn+=Exp_FalseNeg(node_list[0],pop)
    # for ele in range(0,len(node_list)-1):
    #     et+=Exp_Test_Num(node_list[ele+1],node_list[ele])
    #     efp+=Exp_FalsePos(node_list[ele+1],node_list[ele])
    #     efn+=Exp_FalseNeg(node_list[ele+1],node_list[ele])
        
    # print(node_list)
    
    # if (flag==0): 
    #     #0 is for the diluted case
    #     return dist[pop],et,efp,efn
    # if (flag==1):
    #     #1 is for the non-diluted case
    return node_list

def calc_meas(node_list):

    of=0
    of+=weight_D(node_list[0],pop)
    et=0
    et+=Exp_Test_Num_D(node_list[0],pop)
    efp=0
    efp+=Exp_FalsePos_D(node_list[0],pop)
    efn=0
    efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])
        of+=weight_D(node_list[ele+1],node_list[ele])
        
    return of,et,efp,efn

def calc_measof(node_list):

    of=0
    of+=weight_D(node_list[0],pop)
    # et=0
    # et+=Exp_Test_Num_D(node_list[0],pop)
    # efp=0
    # efp+=Exp_FalsePos_D(node_list[0],pop)
    # efn=0
    # efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        # et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        # efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        # efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])
        of+=weight_D(node_list[ele+1],node_list[ele])
        
    return of    

def Bellman(cost,B):   
    d=[[None for i in range(pop+1)] for i in range(B+1)]
    d[0]=[float("inf") for i in range(pop+1)]
    d[0][0]=0.0
    pathNode=[[None for i in range(pop+1)] for i in range(B+1)]
    for k in range(B):
        for i in range(pop+1):
            tempMin=float("inf")
            tempIndex=-1.0
            for j in range(i):
                if d[k][j]+cost[j][i]<tempMin:
                    tempMin=d[k][j]+cost[j][i]
                    tempIndex=j
            if tempMin<d[k][i]:
                d[k+1][i]=tempMin
                pathNode[k+1][i]=tempIndex
            else:
                d[k+1][i]=d[k][i]
                pathNode[k+1][i]=pathNode[k][i]
    
    path=[]
    node=pop
    # path.append(node)
    k=B
    while node>0:
        node=pathNode[k][node]
        path.append(node)
        k=k-1
        
    node_list=path
    # path.reverse()
    
    # node_list=[]
    # a=pathNode[B][pop]
    # while (a!=0):
    #     node_list.append(a)
    #     a=pathNode[B][a]
    # node_list.append(0) 
    # node_list.reverse()

    of=0
    of+=weight_D(node_list[0],pop)
    et=0
    et+=Exp_Test_Num_D(node_list[0],pop)
    efp=0
    efp+=Exp_FalsePos_D(node_list[0],pop)
    efn=0
    efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        of+=weight_D(node_list[ele+1],node_list[ele])
        et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])  
    
    # print(node_list)  
    # print(d[B][pop])
    return of,et,efp,efn,node_list


#%% 
# 6. ESTABLISH CONSTANT PARAMETERS
#population
pop = 100

# #COVID Parameters
# spec = 1
# alpha= 0.012595
# gamma= 0.257138
# const_sens = sens(1,1)

# HBV Parameters
spec = 0.99
alpha= 0.0514
gamma= 0.2416
const_sens = sens(1,1)

#Confidence Interval
CI = 0.05

#Simulations
simulations=100
lambda_1=0
lambda_2=0
#%% 
# 7A. SIMULATIONS Normal RV COVID
           


results = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'delta',
                                 'OF_D_EST',
                                 'OF_D_ACT',
                                 'OF_CS5_EST',
                                 'OF_CS5_ACT',
                                 'OF_CS8_EST',
                                 'OF_CS8_ACT'])
           
for w1 in tqdm(list(range(0, 12, 2)),desc='W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2),position=0,leave=True):
   for w2 in tqdm(list(range(0,1,2))):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10
           
          #  print('\n')
          #  # print('------------------------------------------------------------------')
          #  print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
          #  # print('------------------------------------------------------------------')
          #  print('\n')
           
           of_d_mean_est=np.empty(simulations)
           of_d_act=np.empty(simulations) 
           of_cs8_act=np.empty(simulations)
           of_cs5_act=np.empty(simulations)
           of_cs8_mean_est=np.empty(simulations)
           of_cs5_mean_est=np.empty(simulations)
           for dl in range(0,120,20): 
               for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
                   risk_vec=get_risk_vector()
                   uns_risk_vec=risk_vec
                   perturb_risk_vec=perturb_rv(risk_vec,dl)
                   # CURRENT PRATICE MEAN EST
                   current_node_list=list(range(0,100,8))
                   current_node_list.reverse()
                      
                   of_cs8_mean_est[sim]=calc_measof(current_node_list)
                   
                   current_node_list=list(range(0,100,5))
                   current_node_list.reverse()
                      
                   of_cs5_mean_est[sim]=calc_measof(current_node_list)
                   
                   #CURRENT PRACTICE ACTUAL
                   risk_vec=perturb_risk_vec
                   
                   current_node_list=list(range(0,100,8))
                   current_node_list.reverse()
                      
                   of_cs8_act[sim]=calc_measof(current_node_list)
                   
                   current_node_list=list(range(0,100,5))
                   current_node_list.reverse()
                      
                   of_cs5_act[sim]=calc_measof(current_node_list)
                   
                   
                   #Optimal Scenario
                   risk_vec=uns_risk_vec
                   sor_risk_vec=np.sort(risk_vec)
                   risk_vec=sor_risk_vec
                   
                   wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
                   for i in range(0,pop):
                      for j in range(i+1,pop+1):
                            wt_D[i][j]=weight_D(i,j)            
            
                   s = 0
            
                   g_d = Graph(pop+1) 
          
                   for i in range(0,pop):
                       for j in range(i+1,pop+1):
                           g_d.addEdge(i,j,wt_D[i][j])
                   
                   
                   # of_d[sim]=Bellman(wt_D,20)[0]
                   
                  #  node_list=g_d.shortestPath(s) 
                     # part_d[sim]=len(node_list)
                   # of_d_mean_est[sim]=calc_measof(node_list)
                   of_d_mean_est[sim],*_,node_list=Bellman(wt_D,20)
                   
                   risk_vec=uns_risk_vec
                   per_sort_rv=[perturb_risk_vec[a] for a in np.argsort(risk_vec)]
                   risk_vec=per_sort_rv
                   of_d_act[sim]=calc_measof(node_list)

               # histogram data
               if((lambda_1==0)&(lambda_2==0)&(dl==100)):
                    hist_ETps_est=of_d_mean_est
                    hist_ETps_act=of_d_act
                    hist_ETcs5_est=of_cs5_mean_est
                    hist_ETcs5_act=of_cs5_act
                    hist_ETcs8_est=of_cs8_mean_est
                    hist_ETcs8_act=of_cs8_act

               results = results.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'delta':dl,
                                     'OF_D_EST':round(np.mean(of_d_mean_est),4),
                                     'OF_D_ACT':round(np.mean(of_d_act),4),
                                     'OF_CS5_EST':round(np.mean(of_cs5_mean_est),4),
                                     'OF_CS5_ACT':round(np.mean(of_cs5_act),4),
                                     'OF_CS8_EST':round(np.mean(of_cs8_mean_est),4),
                                     'OF_CS8_ACT':round(np.mean(of_cs8_act),4)},ignore_index=True)
# print("Time:")
# print(time.time() - start_time) 
#%%
results.to_excel('COVID_RiskResilience1002.xlsx')   

#%%
df = pd.read_excel('COVID_ET_HistResultsDelta100.xlsx')
hist_ETps_est=df['hist_ETps_est']
hist_ETps_act=df['hist_ETps_act']
hist_ETcs5_est=df['hist_ETcs5_est']
hist_ETcs5_act=df['hist_ETcs5_act']
hist_ETcs8_est=df['hist_ETcs8_est']
hist_ETcs8_act=df['hist_ETcs8_act']

#%%    
# 7B. PLOT HISTOGRAM
sns.distplot(hist_ETps_est, hist=True, kde=False, 
             bins=10, color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='PS Estimated')
mean_ps_est=np.mean(hist_ETps_est)
plt.axvline(mean_ps_est+.02, c= 'blue',lw=1, ls= '--')

sns.distplot(hist_ETps_act, hist=True, kde=False, 
             bins=14, color = 'indianred', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='PS Actual')
mean_ps_act=np.mean(hist_ETps_act)
plt.axvline(mean_ps_act, c= 'tomato',lw=1, ls= '--')

sns.distplot(hist_ETcs5_est, hist=True, kde=False, 
             bins=7, color = 'darkgreen', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='CS$_{5}$ Estimated')
mean_cs5_est=np.mean(hist_ETcs5_est)
plt.axvline(mean_cs5_est, c= 'green', lw=1, ls= '--')

sns.distplot(hist_ETcs5_act, hist=True, kde=False, 
             bins=13, color = 'darkorange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='CS$_{5}$ Actual')
mean_cs5_act=np.mean(hist_ETcs5_act)
plt.axvline(mean_cs5_act, c= 'orange', lw=1, ls= '--')


# sns.distplot(ET_URDDS, hist=True, kde=False, 
#              bins=10, color = 'darkred', 
#              hist_kws={'edgecolor':'red'},
#              kde_kws={'linewidth': 4}, label='URDDS',axlabel='ET')

# ax_hist.axvline(mean_ps_est, color='blue', linestyle='--')   
# ax_hist.axvline(mean_cs, color='r', linestyle='--')

plt.ylabel("Frequency")
plt.xlabel('Expected number of tests, $\it{E[T]}$')
# plt.axvline(mean_ub, c= 'red', ls='--', label='Mean (URDDS)')

# plt.xlim(0.29, 0.48)
# plt.ylim(0, 30)
# plt.title('Percentage Difference: '+str(round((np.mean(ET_URDDS)-np.mean(opti))/np.mean(opti)*100,3))+'%' ) 

plt.legend(bbox_to_anchor=(0.3, 1),prop={'size': 8.5}, borderaxespad=.5)
# plt.show()
# sns.plt.show()
plt.savefig('COVIDRRHist.png', dpi=500,bbox_inches='tight')    
#%% 7A. SIMULATIONS HBV

results = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'delta',
                                 'OF_D_EST',
                                 'OF_D_ACT',
                                 'OF_CS16_EST',
                                 'OF_CS16_ACT'])

           
for w1 in range(0, 12, 2):
   for w2 in range(0,12,2):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10
           
           print('\n')
           # print('------------------------------------------------------------------')
           print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
           # print('------------------------------------------------------------------')
           print('\n')
           
           of_d_mean_est=np.empty(simulations)
           of_d_act=np.empty(simulations) 
           of_cs16_act=np.empty(simulations)
           of_cs16_mean_est=np.empty(simulations)
           for dl in tqdm(range(100,120,20), desc='delta loop', position=0, leave=True): 
               for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
                   risk_vec=get_risk_vector()
                   uns_risk_vec=risk_vec
                   perturb_risk_vec=perturb_rv(risk_vec,dl)
                   # CURRENT PRATICE MEAN EST
                   current_node_list=list(range(0,100,16))
                   current_node_list.reverse()
                      
                   of_cs16_mean_est[sim]=calc_measof(current_node_list)
                   
                   #CURRENT PRACTICE ACTUAL
                   risk_vec=perturb_risk_vec
                   
                   current_node_list=list(range(0,100,16))
                   current_node_list.reverse()
                      
                   of_cs16_act[sim]=calc_measof(current_node_list)
                   
                   
                   
                   #Optimal Scenario
                   risk_vec=uns_risk_vec
                   sor_risk_vec=np.sort(risk_vec)
                   risk_vec=sor_risk_vec
                   
                   wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
                   for i in range(0,pop):
                      for j in range(i+1,pop+1):
                            wt_D[i][j]=weight_D(i,j)            
            
                   s = 0
            
                   g_d = Graph(pop+1) 
          
                   for i in range(0,pop):
                       for j in range(i+1,pop+1):
                           g_d.addEdge(i,j,wt_D[i][j])
                   
                   
                   # of_d[sim]=Bellman(wt_D,20)[0]
                   
                  #  node_list=g_d.shortestPath(s) 
                    #  part_d[sim]=len(node_list)
                   # of_d_mean_est[sim]=calc_measof(node_list)
                  #  of_d_mean_est[sim]=Bellman(wt_D,8)[0]
                   of_d_mean_est[sim],*_,node_list=Bellman(wt_D,8)
                   
                   risk_vec=uns_risk_vec
                   per_sort_rv=[perturb_risk_vec[a] for a in np.argsort(risk_vec)]
                   risk_vec=per_sort_rv
                   of_d_act[sim]=calc_measof(node_list)

               # histogram data
               if((lambda_1==0)&(lambda_2==0)&(dl==100)):
                    hist_ETps_est=of_d_mean_est
                    hist_ETps_act=of_d_act
                    hist_ETcs16_est=of_cs16_mean_est
                    hist_ETcs16_act=of_cs16_act 
                    break 

               results = results.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'delta':dl,
                                     'OF_D_EST':round(np.mean(of_d_mean_est),4),
                                     'OF_D_ACT':round(np.mean(of_d_act),4),
                                     'OF_CS16_EST':round(np.mean(of_cs16_mean_est),4),
                                     'OF_CS16_ACT':round(np.mean(of_cs16_act),4)},ignore_index=True)   
#%% 7B. PLOT HISTOGRAM
sns.distplot(hist_ETps_est, hist=True, kde=False, 
             bins=10, color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='PS Estimated')
mean_ps_est=np.mean(hist_ETps_est)
plt.axvline(mean_ps_est, c= 'blue',lw=1, ls= '--')

sns.distplot(hist_ETps_act, hist=True, kde=False, 
             bins=14, color = 'indianred', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='PS Actual')
mean_ps_act=np.mean(hist_ETps_act)
plt.axvline(mean_ps_act, c= 'tomato',lw=1, ls= '--')

sns.distplot(hist_ETcs16_est, hist=True, kde=False, 
             bins=10, color = 'darkgreen', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='CS$_{16}$ Estimated')
mean_cs16_est=np.mean(hist_ETcs16_est)
plt.axvline(mean_cs16_est, c= 'green', lw=1, ls= '--')

sns.distplot(hist_ETcs16_act, hist=True, kde=False, 
             bins=13, color = 'darkorange', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='CS$_{16}$ Actual')
mean_cs16_act=np.mean(hist_ETcs16_act)
plt.axvline(mean_cs16_act, c= 'orange', lw=1, ls= '--')


# sns.distplot(ET_URDDS, hist=True, kde=False, 
#              bins=10, color = 'darkred', 
#              hist_kws={'edgecolor':'red'},
#              kde_kws={'linewidth': 4}, label='URDDS',axlabel='ET')

# ax_hist.axvline(mean_ps_est, color='blue', linestyle='--')   
# ax_hist.axvline(mean_cs, color='r', linestyle='--')

plt.ylabel("Frequency")
plt.xlabel('Expected number of tests, $\it{E[T]}$')
# plt.axvline(mean_ub, c= 'red', ls='--', label='Mean (URDDS)')

plt.xlim(10.5, 15)
# plt.ylim(0, 30)
# plt.title('Percentage Difference: '+str(round((np.mean(ET_URDDS)-np.mean(opti))/np.mean(opti)*100,3))+'%' ) 

plt.legend(bbox_to_anchor=(0.27, 1),prop={'size': 7.5}, borderaxespad=.5)
# plt.show()
# sns.plt.show()
plt.savefig('HBVRRHist.png', dpi=500,bbox_inches='tight')   
#%%
df=pd.DataFrame()
df[0]=of_d_mean_est
df[1]=of_d_act
df[2]=of_cs5_mean_est
df[3]=of_cs5_act
df[4]=of_cs8_mean_est
df[5]=of_cs8_act
df.to_excel("Delta1530.xlsx")
#%%

bins = np.linspace(20, 35, 100)

pyplot.hist(df[0], bins, alpha=0.5, label='mean_est')
pyplot.hist(df[1], bins, alpha=0.5, label='act')
pyplot.legend(loc='upper right')
pyplot.show()
#%% 7B. SIMULATIONS Perturbed RV
           


results_per = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'OF_D',
                                 'CI_OFD']) 
           
for w1 in tqdm(list(range(0, 12, 2)),desc='lambda loop',position=0,leave=True):
   for w2 in list(range(0,1,1)):
       if(w1+w2==0):
           lambda_1=w1/10
           lambda_2=w2/10
           
           print('\n')
           # print('------------------------------------------------------------------')
           print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
           # print('------------------------------------------------------------------')
           print('\n')
           
           of_d_per=np.empty(simulations)
            
            
           for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
               risk_vec=get_risk_vector()
               risk_vec=perturb_rv(risk_vec,30)
               #Optimal Scenario
               risk_vec=np.sort(risk_vec)
               
               wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
               for i in range(0,pop):
                  for j in range(i+1,pop+1):
                        wt_D[i][j]=weight_D(i,j)            
        
               s = 0
        
               g_d = Graph(pop+1) 
      
               for i in range(0,pop):
                   for j in range(i+1,pop+1):
                       g_d.addEdge(i,j,wt_D[i][j])
               
               
               # of_d[sim]=Bellman(wt_D,20)[0]
               
               node_list=g_d.shortestPath(s) 
                 # part_d[sim]=len(node_list)
               of_d_per[sim]=calc_measof(node_list)


           results_per = results_per.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'OF_D':round(np.mean(of_d_per),4),
                                     'CI_OFD':1.96*np.std(of_d_per)/np.size(of_d_per)},ignore_index=True)
print("Time:")
print(time.time() - start_time) 
#%% 7. SIMULATIONS
results = pd.DataFrame(columns= ['I',
                                 'OF_CS']) 
for i in range(1,101,1):
    of_cs=np.empty(simulations)
    for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
               risk_vec=get_risk_vector()
               current_node_list=list(range(0,100,i))
               current_node_list.reverse()
               of_cs[sim]=calc_measof(current_node_list)
    results = results.append({'I':i,
                          'OF_CS':round(np.mean(of_cs),4)},ignore_index=True)


#%% 7. SIMULATIONS
# sys.stdout = open("FINALRESULTS.txt", "w")
start_time=time.time()
results = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'OF_D',
                                 'CI_OFD',
                                 'OF_CS8',
                                 'CI_OFCS8',
                                 'OF_CS5',
                                 'CI_OFCS5']) 
           
for w1 in tqdm(list(range(0, 12, 2)),desc='lambda loop',position=0,leave=True):
   for w2 in list(range(0,1,1)):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10

           print('\n')
           # print('------------------------------------------------------------------')
           print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
           # print('------------------------------------------------------------------')
           print('\n')

           
           of_d=np.empty(simulations)
           of_cs8=np.empty(simulations)
           of_cs5=np.empty(simulations)

           
           for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
               risk_vec=get_risk_vector()
               
               #Current Scenario
               # if ((lambda_1+lambda_2) == 1):
               #     current_node_list=list(range(0,100,1))
               #     current_node_list.reverse()
               # else:
               current_node_list=list(range(0,100,8))
               current_node_list.reverse()
                  
               of_cs8[sim]=calc_measof(current_node_list)
               
               current_node_list=list(range(0,100,5))
               current_node_list.reverse()
                  
               of_cs5[sim]=calc_measof(current_node_list)
               
               #Optimal Scenario
               risk_vec=np.sort(risk_vec)
               
               wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
               for i in range(0,pop):
                  for j in range(i+1,pop+1):
                        wt_D[i][j]=weight_D(i,j)            
        
               s = 0
        
               g_d = Graph(pop+1) 
      
               for i in range(0,pop):
                   for j in range(i+1,pop+1):
                       g_d.addEdge(i,j,wt_D[i][j])
               
               
               # of_d[sim]=Bellman(wt_D,20)[0]
               
               node_list=g_d.shortestPath(s) 
                 # part_d[sim]=len(node_list)
               of_d[sim]=calc_measof(node_list)


           results = results.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'OF_D':round(np.mean(of_d),4),
                                     'CI_OFD':1.96*np.std(of_d)/np.size(of_d),
                                     'OF_CS8':round(np.mean(of_cs8),4),
                                     'CI_OFCS8':1.96*np.std(of_cs8)/np.size(of_cs8),
                                     'OF_CS5':round(np.mean(of_cs5),4),
                                     'CI_OFCS5':1.96*np.std(of_cs5)/np.size(of_cs5)},ignore_index=True)
print("Time:")
print(time.time() - start_time) 







    
    
    #%% 5. SORTING AND SHORTEST PATH
'''
TOPOLOGICAL SORT ALGORITHM
FOR SHORTEST PATH
'''
from collections import defaultdict

'''
DEFINING GRAPH DICTIONARY 
AND RELATED FUNCTIONS
'''
# Graph is represented using adjacency list. Every 
# node of adjacency list contains vertex number of 
# the vertex to which edge connects. It also contains 
# weight of the edge
class Graph:
  def __init__(self,vertices):

    # No. of vertices 
    self.V= vertices

    # dictionary containing adjacency List 
    self.graph = defaultdict(list)

  # function to add an edge to graph 
  def addEdge(self,u,v,w):
    self.graph[u].append((v,w))


  # A recursive function used by shortestPath function
  def topologicalSortUtil(self,v,visited,stack):

    # Mark the current node as visited.
    visited[v] = True


    # Recur for all the vertices adjacent to this vertex
    if v in self.graph.keys():
      for node,weight in self.graph[v]:
        if visited[node] == False:
          self.topologicalSortUtil(node,visited,stack)

    # Push current vertex to stack which stores topological sort
    stack.append(v)


  #The function to find shortest paths from given vertex. 
  #It uses recursive topologicalSortUtil() to get topological 
  #sorting of given graph.
  def shortestPath(self,s):

    # Mark all the vertices as not visited 
    visited = [False]*self.V
    stack=[]

    # Call the recursive helper function to store Topological 
    # Sort starting from source vertex 
    for i in range(self.V):
      if visited[i] == False:
        self.topologicalSortUtil(s,visited,stack)

    # Initialize distances to all vertices as infinite and 
    # distance to source as 0
    dist = [float("Inf")]*(self.V)
    dist[s] = 0
    part = [0]*(self.V)

    # print(dist)
    # print(stack)

    # Process vertices in topological order 
    while stack:

      # Get the next vertex from topological order
      i=stack.pop()

      # Update distances of all adjacent vertices
      for node,weight in self.graph[i]:
        if dist[node] > dist[i] + weight:
          dist[node] = dist[i] + weight
          part[node] = (node,i)

    # Print the distance from source to Node_100
    # for i in range(self.V):
    # print ("Dist - %f, Node - %d" %(dist[100],100)) #if dist[i] != float ("Inf") else "Inf" ,
    node_list=[]
    a=part[pop]
    while (a[1]!=0):
        node_list.append(a[1])
        a=part[a[1]]
    node_list.append(0)  
    
    # et=0
    # et+=Exp_Test_Num(node_list[0],pop)
    # efp=0
    # efp+=Exp_FalsePos(node_list[0],pop)
    # efn=0
    # efn+=Exp_FalseNeg(node_list[0],pop)
    # for ele in range(0,len(node_list)-1):
    #     et+=Exp_Test_Num(node_list[ele+1],node_list[ele])
    #     efp+=Exp_FalsePos(node_list[ele+1],node_list[ele])
    #     efn+=Exp_FalseNeg(node_list[ele+1],node_list[ele])
        
    # print(node_list)
    
    # if (flag==0): 
    #     #0 is for the diluted case
    #     return dist[pop],et,efp,efn
    # if (flag==1):
    #     #1 is for the non-diluted case
    return node_list

def calc_meas(node_list):

    of=0
    of+=weight_D(node_list[0],pop)
    et=0
    et+=Exp_Test_Num_D(node_list[0],pop)
    efp=0
    efp+=Exp_FalsePos_D(node_list[0],pop)
    efn=0
    efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])
        of+=weight_D(node_list[ele+1],node_list[ele])
        
    return of,et,efp,efn

def calc_measof(node_list):

    of=0
    of+=weight_D(node_list[0],pop)
    # et=0
    # et+=Exp_Test_Num_D(node_list[0],pop)
    # efp=0
    # efp+=Exp_FalsePos_D(node_list[0],pop)
    # efn=0
    # efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        # et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        # efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        # efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])
        of+=weight_D(node_list[ele+1],node_list[ele])
        
    return of


def Bellmanof(cost,B):   
    d=[[None for i in range(pop+1)] for i in range(B+1)]
    d[0]=[float("inf") for i in range(pop+1)]
    d[0][0]=0.0
    pathNode=[[None for i in range(pop+1)] for i in range(B+1)]
    for k in range(B):
        for i in range(pop+1):
            tempMin=float("inf")
            tempIndex=-1.0
            for j in range(i):
                if d[k][j]+cost[j][i]<tempMin:
                    tempMin=d[k][j]+cost[j][i]
                    tempIndex=j
            if tempMin<d[k][i]:
                d[k+1][i]=tempMin
                pathNode[k+1][i]=tempIndex
            else:
                d[k+1][i]=d[k][i]
                pathNode[k+1][i]=pathNode[k][i]
    
    path=[]
    node=pop
    # path.append(node)
    k=B
    while node>0:
        node=pathNode[k][node]
        path.append(node)
        k=k-1
        
    node_list=path
    # path.reverse()
    
    # node_list=[]
    # a=pathNode[B][pop]
    # while (a!=0):
    #     node_list.append(a)
    #     a=pathNode[B][a]
    # node_list.append(0) 
    # node_list.reverse()

    of=0
    of+=weight_D(node_list[0],pop)
    # et=0
    # et+=Exp_Test_Num_D(node_list[0],pop)
    # efp=0
    # efp+=Exp_FalsePos_D(node_list[0],pop)
    # efn=0
    # efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        of+=weight_D(node_list[ele+1],node_list[ele])
        # et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        # efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        # efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])  
    
    # print(node_list)  
    # print(d[B][pop])
    return of,len(node_list)


def Bellmanmodof(m,cost,B):   
    d=[[None for i in range(pop+1)] for i in range(B+1)]
    d[0]=[float("inf") for i in range(pop+1)]
    d[0][0]=0.0
    pathNode=[[None for i in range(pop+1)] for i in range(B+1)]
    for k in range(B):
        for i in range(pop+1):
            tempMin=float("inf")
            tempIndex=-1.0
            for j in range(i):
                if (i-j)>m:
                    continue
                elif d[k][j]+cost[j][i]<tempMin:
                    tempMin=d[k][j]+cost[j][i]
                    tempIndex=j
            if tempMin<d[k][i]:
                d[k+1][i]=tempMin
                pathNode[k+1][i]=tempIndex
            else:
                d[k+1][i]=d[k][i]
                pathNode[k+1][i]=pathNode[k][i]
    
    path=[]
    node=pop
    # path.append(node)
    k=B
    while node>0:
        node=pathNode[k][node]
        path.append(node)
        k=k-1
        
    node_list=path
    # path.reverse()
    
    # node_list=[]
    # a=pathNode[B][pop]
    # while (a!=0):
    #     node_list.append(a)
    #     a=pathNode[B][a]
    # node_list.append(0) 
    # node_list.reverse()

    of=0
    of+=weight_D(node_list[0],pop)
    # et=0
    # et+=Exp_Test_Num_D(node_list[0],pop)
    # efp=0
    # efp+=Exp_FalsePos_D(node_list[0],pop)
    # efn=0
    # efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        of+=weight_D(node_list[ele+1],node_list[ele])
        # et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        # efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        # efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])  
    
    # print(node_list)  
    # print(d[B][pop])
    return of,len(node_list)

def Bellman(cost,B):   
    d=[[None for i in range(pop+1)] for i in range(B+1)]
    d[0]=[float("inf") for i in range(pop+1)]
    d[0][0]=0.0
    pathNode=[[None for i in range(pop+1)] for i in range(B+1)]
    for k in range(B):
        for i in range(pop+1):
            tempMin=float("inf")
            tempIndex=-1.0
            for j in range(i):
                if d[k][j]+cost[j][i]<tempMin:
                    tempMin=d[k][j]+cost[j][i]
                    tempIndex=j
            if tempMin<d[k][i]:
                d[k+1][i]=tempMin
                pathNode[k+1][i]=tempIndex
            else:
                d[k+1][i]=d[k][i]
                pathNode[k+1][i]=pathNode[k][i]
    
    path=[]
    node=pop
    # path.append(node)
    k=B
    while node>0:
        node=pathNode[k][node]
        path.append(node)
        k=k-1
        
    node_list=path
    # path.reverse()
    
    # node_list=[]
    # a=pathNode[B][pop]
    # while (a!=0):
    #     node_list.append(a)
    #     a=pathNode[B][a]
    # node_list.append(0) 
    # node_list.reverse()

    of=0
    of+=weight_D(node_list[0],pop)
    et=0
    et+=Exp_Test_Num_D(node_list[0],pop)
    efp=0
    efp+=Exp_FalsePos_D(node_list[0],pop)
    efn=0
    efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        of+=weight_D(node_list[ele+1],node_list[ele])
        et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])  
    
    # print(node_list)  
    # print(d[B][pop])
    return of,et,efp,efn,len(node_list)


#%% 6. ESTABLISH CONSTANT PARAMETERS
#population
pop = 100

#COVID Parameters
spec = 1
alpha= 0.012595
gamma= 0.257138
const_sens = sens(1,1)

#Confidence Interval
CI = 0.05

#Simulations
simulations=1000
lambda_1=0
lambda_2=0
#%% 7. SIMULATIONS
results = pd.DataFrame(columns= ['I',
                                 'OF_CS']) 
for i in range(1,101,1):
    of_cs=np.empty(simulations)
    for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
               risk_vec=get_risk_vector()
               current_node_list=list(range(0,100,i))
               current_node_list.reverse()
               of_cs[sim]=calc_measof(current_node_list)
    results = results.append({'I':i,
                          'OF_CS':round(np.mean(of_cs),4)},ignore_index=True)


#%% 7. SIMULATIONS
# sys.stdout = open("FINALRESULTS.txt", "w")
start_time=time.time()
results = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'OF_D',
                                 'CI_OFD',
                                 'OF_CS8',
                                 'CI_OFCS8',
                                 'OF_CS5',
                                 'CI_OFCS5']) 
           
for w1 in tqdm(list(range(0, 12, 2)),desc='lambda loop',position=0,leave=True):
   for w2 in list(range(0,1,1)):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10

           print('\n')
           # print('------------------------------------------------------------------')
           print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
           # print('------------------------------------------------------------------')
           print('\n')

           
           of_d=np.empty(simulations)
           of_cs8=np.empty(simulations)
           of_cs5=np.empty(simulations)

           
           for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
               risk_vec=get_risk_vector()
               
               #Current Scenario
               # if ((lambda_1+lambda_2) == 1):
               #     current_node_list=list(range(0,100,1))
               #     current_node_list.reverse()
               # else:
               current_node_list=list(range(0,100,8))
               current_node_list.reverse()
                  
               of_cs8[sim]=calc_measof(current_node_list)
               
               current_node_list=list(range(0,100,5))
               current_node_list.reverse()
                  
               of_cs5[sim]=calc_measof(current_node_list)
               
               #Optimal Scenario
               risk_vec=np.sort(risk_vec)
               
               wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
               for i in range(0,pop):
                  for j in range(i+1,pop+1):
                        wt_D[i][j]=weight_D(i,j)            
        
               s = 0
        
               g_d = Graph(pop+1) 
      
               for i in range(0,pop):
                   for j in range(i+1,pop+1):
                       g_d.addEdge(i,j,wt_D[i][j])
               
               
               # of_d[sim]=Bellman(wt_D,20)[0]
               
               node_list=g_d.shortestPath(s) 
                 # part_d[sim]=len(node_list)
               of_d[sim]=calc_measof(node_list)


           results = results.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'OF_D':round(np.mean(of_d),4),
                                     'CI_OFD':1.96*np.std(of_d)/np.size(of_d),
                                     'OF_CS8':round(np.mean(of_cs8),4),
                                     'CI_OFCS8':1.96*np.std(of_cs8)/np.size(of_cs8),
                                     'OF_CS5':round(np.mean(of_cs5),4),
                                     'CI_OFCS5':1.96*np.std(of_cs5)/np.size(of_cs5)},ignore_index=True)
print("Time:")
print(time.time() - start_time) 


# #%% 7. SIMULATIONS
# # sys.stdout = open("FINALRESULTS.txt", "w")
# start_time=time.time()
# results = pd.DataFrame(columns= ['w1',
#                                  'w2',
#                                  'OF_D',
#                                  'CI_OFD',
#                                  'OF_CS8',
#                                  'CI_OFCS8',
#                                  'OF_CS5',
#                                  'CI_OFCS5']) 
           
# for w1 in tqdm(list(range(0, 12, 1)),desc='lambda loop',position=0,leave=True):
#    for w2 in list(range(0,1,1)):
#        if(w1+w2<=10):
#            lambda_1=w1/10
#            lambda_2=w2/10

#            print('\n')
#            # print('------------------------------------------------------------------')
#            print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
#            # print('------------------------------------------------------------------')
#            print('\n')

           
#            of_d=np.empty(simulations)
#            of_cs8=np.empty(simulations)
#            of_cs5=np.empty(simulations)

           
#            for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
#                risk_vec=get_risk_vector()
               
#                #Current Scenario
#                # if ((lambda_1+lambda_2) == 1):
#                #     current_node_list=list(range(0,100,1))
#                #     current_node_list.reverse()
#                # else:
#                current_node_list=list(range(0,100,8))
#                current_node_list.reverse()
                  
#                of_cs8[sim]=calc_measof(current_node_list)
               
#                current_node_list=list(range(0,100,5))
#                current_node_list.reverse()
                  
#                of_cs5[sim]=calc_measof(current_node_list)
               
#                #Optimal Scenario
#                risk_vec.sort()
               
#                wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
#                for i in range(0,pop):
#                   for j in range(i+1,pop+1):
#                         wt_D[i][j]=weight_D(i,j)            
        
#                s = 0
        
#                g_d = Graph(pop+1) 
      
#                for i in range(0,pop):
#                    for j in range(i+1,pop+1):
#                        g_d.addEdge(i,j,wt_D[i][j])
               
                
#                node_list=g_d.shortestPath(s) 
#                 # part_d[sim]=len(node_list)
#                of_d[sim]=calc_measof(node_list)


#            results = results.append({'w1':lambda_1,
#                                      'w2':lambda_2,
#                                      'OF_D':round(np.mean(of_d),4),
#                                      'CI_OFD':1.96*np.std(of_d)/np.size(of_d),
#                                      'OF_CS8':round(np.mean(of_cs8),4),
#                                      'CI_OFCS8':1.96*np.std(of_cs8)/np.size(of_cs8),
#                                      'OF_CS5':round(np.mean(of_cs5),4),
#                                      'CI_OFCS5':1.96*np.std(of_cs5)/np.size(of_cs5)},ignore_index=True)
# print("Time:")
# print(time.time() - start_time) 
#%% 8. OUTPUT TO CSV
results.to_csv('Covid_CSvsPS_100_B20_50sims.csv')      
    


# for w1 in tqdm(range(0, 12, 1),desc='lambda loop',position=0,leave=True):
#    for w2 in range(0,12,1):
#        if(w1+w2<=10):
#            lambda_1=w1/10
#            lambda_2=w2/10
           


#            print('\n')
#            print('------------------------------------------------------------------')
#            print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
#            print('------------------------------------------------------------------')
#            print('\n')
        
#            of_nd=np.empty(simulations)
#            et_nd=np.empty(simulations)
#            efp_nd=np.empty(simulations)
#            efn_nd=np.empty(simulations)
#            part_nd=np.empty(simulations)
           
#            of_d=np.empty(simulations)
#            et_d=np.empty(simulations)
#            efp_d=np.empty(simulations)
#            efn_d=np.empty(simulations)
#            part_d=np.empty(simulations)
           
#            of_cs=np.empty(simulations)
#            et_cs=np.empty(simulations)
#            efp_cs=np.empty(simulations)
#            efn_cs=np.empty(simulations)
           
#            for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
#                get_risk_vector()
               
#                #Current Scenario
#                current_node_list=list(range(0,100,16))
#                current_node_list.reverse()
               
#                of_cs[sim],et_cs[sim],efp_cs[sim],efn_cs[sim]=calc_meas(current_node_list)
               
#                #Optimal Scenario
#                risk_vec.sort()
#                g_nd = Graph(pop+1)
               
#                wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
#                for i in range(0,pop):
#                   for j in range(i+1,pop+1):
#                         wt_D[i][j]=weight_D(i,j)
            
#                wt_ND =[[None for i in range(pop+1)]for i in range(pop+1)]
#                for i in range(0,pop):
#                   for j in range(i+1,pop+1):
#                         wt_ND[i][j]=weight_ND(i,j)
                        
#                for i in range(0,pop):
#                    for j in range(i+1,pop+1):
#                        g_nd.addEdge(i,j,wt_ND[i][j])
        
#                s = 0
        
#                node_list=g_nd.shortestPath(s) 
#                part_nd[sim]=len(node_list)
#                of_nd[sim],et_nd[sim],efp_nd[sim],efn_nd[sim]=calc_meas(node_list)
        
#                g_d = Graph(pop+1) 
      
#                for i in range(0,pop):
#                    for j in range(i+1,pop+1):
#                        g_d.addEdge(i,j,wt_D[i][j])
                
#                node_list=g_d.shortestPath(s) 
#                part_d[sim]=len(node_list)
#                of_d[sim],et_d[sim],efp_d[sim],efn_d[sim]=calc_meas(node_list)
               

#            results = results.append({'w1':lambda_1,
#                                      'w2':lambda_2,
#                                      'OF_ND':round(np.mean(of_nd),4),
                                     
#                                      'OF_ND_HW':np.std(of_nd)*1.96/np.size(of_nd),
                                     
#                                      'OF_D':round(np.mean(of_d),4),
                                     
#                                      'OF_D_HW':np.std(of_d)*1.96/np.size(of_d),

                                     
#                                      'OF_CS':round(np.mean(of_cs),4),
                                     
#                                      'OF_CS_HW':np.std(of_cs)*1.96/np.size(of_cs),
                                     
#                                      'ET_ND':round(np.mean(et_nd),4),

#                                      'ET_ND_HW':np.std(et_nd)*1.96/np.size(et_nd),
                                     
#                                      'ET_D':round(np.mean(et_d),4),
                                     
#                                      'ET_D_HW':np.std(et_d)*1.96/np.size(et_d),

                                     
#                                      'ET_CS':round(np.mean(et_cs),4),
                                     
#                                      'ET_CS_HW':np.std(et_cs)*1.96/np.size(et_cs),

#                                      'EFP_ND':round(np.mean(efp_nd),4),
   
#                                      'EFP_ND_HW':np.std(efp_nd)*1.96/np.size(efp_nd),
                                     
                                     
#                                      'EFP_D':round(np.mean(efp_d),4),
                                     
#                                      'EFP_D_HW':np.std(efp_d)*1.96/np.size(efp_d),

                                     
#                                      'EFP_CS':round(np.mean(efp_cs),4),
                                     
#                                      'EFP_CS_HW':np.std(efp_cs)*1.96/np.size(efp_cs),
                                     
#                                      'EFN_ND':round(np.mean(efn_nd),4),

#                                      'EFN_ND_HW':np.std(efn_nd)*1.96/np.size(efn_nd),
                                     
#                                      'EFN_D':round(np.mean(efn_d),4),
                                     
#                                      'EFN_D_HW':np.std(efn_d)*1.96/np.size(efn_d),

                                     
#                                      'EFN_CS':round(np.mean(efn_cs),4),
                                     
#                                      'EFN_CS_HW':np.std(efn_cs)*1.96/np.size(efn_cs),
                                     
#                                      'PART_ND':round(np.mean(part_nd),4),

                                     
#                                      'PART_D':round(np.mean(part_d),4)},ignore_index=True)


#%% 7. SIMULATIONS
# sys.stdout = open("FINALRESULTS.txt", "w")
results = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'OF_ND',
                                 'OF_ND_B2',
                                 'OF_ND_B3',
                                 'OF_ND_B4',
                                 'OF_ND_B5',
                                 'OF_ND_B7',
                                 'OF_ND_B7_HW',
                                 'OF_D',
                                 'OF_D_B2',
                                 'OF_D_B3',
                                 'OF_D_B4',
                                 'OF_D_B5',
                                 'OF_D_B7',
                                 'OF_D_B7_HW',
                                 'OF_CS',
                                 'ET_ND',
                                 'ET_ND_B2',
                                 'ET_ND_B3',
                                 'ET_ND_B4',
                                 'ET_ND_B5',
                                 'ET_ND_B7',
                                 'ET_ND_B7_HW',
                                 'ET_D',
                                 'ET_D_B2',
                                 'ET_D_B3',
                                 'ET_D_B4',
                                 'ET_D_B5',
                                 'ET_D_B7',
                                 'ET_D_B7_HW',
                                 'ET_CS',
                                 'EFP_ND',
                                 'EFP_ND_B2',
                                 'EFP_ND_B3',
                                 'EFP_ND_B4',
                                 'EFP_ND_B5',
                                 'EFP_ND_B7',
                                 'EFP_ND_B7_HW',
                                 'EFP_D',
                                 'EFP_D_B2',
                                 'EFP_D_B3',
                                 'EFP_D_B4',
                                 'EFP_D_B5',
                                 'EFP_D_B7',
                                 'EFP_D_B7_HW',
                                 'EFP_CS',
                                 'EFN_ND',
                                 'EFN_ND_B2',
                                 'EFN_ND_B3',
                                 'EFN_ND_B4',
                                 'EFN_ND_B5',
                                 'EFN_ND_B7',
                                 'EFN_ND_B7_HW',
                                 'EFN_D',
                                 'EFN_D_B2',
                                 'EFN_D_B3',
                                 'EFN_D_B4',
                                 'EFN_D_B5',
                                 'EFN_D_B7',
                                 'EFN_D_B7_HW',
                                 'EFN_CS',
                                 'PART_ND',
                                 'PART_ND_B2',
                                 'PART_ND_B3',
                                 'PART_ND_B4',
                                 'PART_ND_B5',
                                 'PART_ND_B10',
                                 'PART_D',
                                 'PART_D_B2',
                                 'PART_D_B3',
                                 'PART_D_B4',
                                 'PART_D_B5',
                                 'PART_D_B10']) 
           
for w1 in tqdm(range(0, 12,10 ),desc='lambda loop',position=0,leave=True):
   for w2 in range(0,12,10):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10
           


           print('\n')
           print('------------------------------------------------------------------')
           print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
           print('------------------------------------------------------------------')
           print('\n')
        
           of_nd=np.empty(simulations)
           et_nd=np.empty(simulations)
           efp_nd=np.empty(simulations)
           efn_nd=np.empty(simulations)
           part_nd=np.empty(simulations)
           
           of_nd_b2=np.empty(simulations)
           et_nd_b2=np.empty(simulations)
           efp_nd_b2=np.empty(simulations)
           efn_nd_b2=np.empty(simulations)
           part_nd_b2=np.empty(simulations)
           
           of_nd_b3=np.empty(simulations)
           et_nd_b3=np.empty(simulations)
           efp_nd_b3=np.empty(simulations)
           efn_nd_b3=np.empty(simulations)
           part_nd_b3=np.empty(simulations)
           
           of_nd_b4=np.empty(simulations)
           et_nd_b4=np.empty(simulations)
           efp_nd_b4=np.empty(simulations)
           efn_nd_b4=np.empty(simulations)
           part_nd_b4=np.empty(simulations)
           
           of_nd_b5=np.empty(simulations)
           et_nd_b5=np.empty(simulations)
           efp_nd_b5=np.empty(simulations)
           efn_nd_b5=np.empty(simulations)
           part_nd_b5=np.empty(simulations)
           
           of_nd_b10=np.empty(simulations)
           et_nd_b10=np.empty(simulations)
           efp_nd_b10=np.empty(simulations)
           efn_nd_b10=np.empty(simulations)
           part_nd_b10=np.empty(simulations)
           
           of_d=np.empty(simulations)
           et_d=np.empty(simulations)
           efp_d=np.empty(simulations)
           efn_d=np.empty(simulations)
           part_d=np.empty(simulations)
           
           of_d_b2=np.empty(simulations)
           et_d_b2=np.empty(simulations)
           efp_d_b2=np.empty(simulations)
           efn_d_b2=np.empty(simulations)
           part_d_b2=np.empty(simulations)
           
           of_d_b3=np.empty(simulations)
           et_d_b3=np.empty(simulations)
           efp_d_b3=np.empty(simulations)
           efn_d_b3=np.empty(simulations)
           part_d_b3=np.empty(simulations)
           
           of_d_b4=np.empty(simulations)
           et_d_b4=np.empty(simulations)
           efp_d_b4=np.empty(simulations)
           efn_d_b4=np.empty(simulations)
           part_d_b4=np.empty(simulations)
           
           of_d_b5=np.empty(simulations)
           et_d_b5=np.empty(simulations)
           efp_d_b5=np.empty(simulations)
           efn_d_b5=np.empty(simulations)
           part_d_b5=np.empty(simulations)
           
           of_d_b10=np.empty(simulations)
           et_d_b10=np.empty(simulations)
           efp_d_b10=np.empty(simulations)
           efn_d_b10=np.empty(simulations)
           part_d_b10=np.empty(simulations)
           
           of_cs=np.empty(simulations)
           et_cs=np.empty(simulations)
           efp_cs=np.empty(simulations)
           efn_cs=np.empty(simulations)
           # part_cs=np.empty(simulations)
           
           for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
               get_risk_vector()
               
               #Current Scenario
               current_node_list=list(range(0,100,16))
               current_node_list.reverse()
               
               of_cs[sim],et_cs[sim],efp_cs[sim],efn_cs[sim]=calc_meas(current_node_list)
               
               #Optimal Scenario
               risk_vec.sort()
               g_nd = Graph(pop+1)
               
               wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
               for i in range(0,pop):
                  for j in range(i+1,pop+1):
                        wt_D[i][j]=weight_D(i,j)
            
               wt_ND =[[None for i in range(pop+1)]for i in range(pop+1)]
               for i in range(0,pop):
                  for j in range(i+1,pop+1):
                        wt_ND[i][j]=weight_ND(i,j)
                        
               for i in range(0,pop):
                   for j in range(i+1,pop+1):
                       g_nd.addEdge(i,j,wt_ND[i][j])
        
               s = 0
        
               node_list=g_nd.shortestPath(s) 
               part_nd[sim]=len(node_list)
               of_nd[sim],et_nd[sim],efp_nd[sim],efn_nd[sim]=calc_meas(node_list)
        
               g_d = Graph(pop+1) 
      
               for i in range(0,pop):
                   for j in range(i+1,pop+1):
                       g_d.addEdge(i,j,wt_D[i][j])
                
               node_list=g_d.shortestPath(s) 
               part_d[sim]=len(node_list)
               of_d[sim],et_d[sim],efp_d[sim],efn_d[sim]=calc_meas(node_list)
               
               of_nd_b2[sim],et_nd_b2[sim],efp_nd_b2[sim],efn_nd_b2[sim],part_nd_b2[sim]=Bellman(wt_ND,2)
               of_nd_b3[sim],et_nd_b3[sim],efp_nd_b3[sim],efn_nd_b3[sim],part_nd_b3[sim]=Bellman(wt_ND,3)
               of_nd_b4[sim],et_nd_b4[sim],efp_nd_b4[sim],efn_nd_b4[sim],part_nd_b4[sim]=Bellman(wt_ND,4)
               of_nd_b5[sim],et_nd_b5[sim],efp_nd_b5[sim],efn_nd_b5[sim],part_nd_b5[sim]=Bellman(wt_ND,5)
               of_nd_b10[sim],et_nd_b10[sim],efp_nd_b10[sim],efn_nd_b10[sim],part_nd_b10[sim]=Bellman(wt_ND,7)

               of_d_b2[sim],et_d_b2[sim],efp_d_b2[sim],efn_d_b2[sim],part_d_b2[sim]=Bellman(wt_D,2)
               of_d_b3[sim],et_d_b3[sim],efp_d_b3[sim],efn_d_b3[sim],part_d_b3[sim]=Bellman(wt_D,3)
               of_d_b4[sim],et_d_b4[sim],efp_d_b4[sim],efn_d_b4[sim],part_d_b4[sim]=Bellman(wt_D,4)
               of_d_b5[sim],et_d_b5[sim],efp_d_b5[sim],efn_d_b5[sim],part_d_b5[sim]=Bellman(wt_D,5)
               of_d_b10[sim],et_d_b10[sim],efp_d_b10[sim],efn_d_b10[sim],part_d_b10[sim]=Bellman(wt_D,7)

           results = results.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'OF_ND':round(np.mean(of_nd),4),
                                     'OF_ND_B2':round(np.mean(of_nd_b2),4),
                                     'OF_ND_B3':round(np.mean(of_nd_b3),4),
                                     'OF_ND_B4':round(np.mean(of_nd_b4),4),
                                     'OF_ND_B5':round(np.mean(of_nd_b5),4),
                                     'OF_ND_B7':round(np.mean(of_nd_b10),4),
                                     'OF_ND_B7_HW':np.std(of_nd_b10)*1.96/np.size(of_nd_b10),
                                     
                                     
                                     'OF_D':round(np.mean(of_d),4),
                                     'OF_D_B2':round(np.mean(of_d_b2),4),
                                     'OF_D_B3':round(np.mean(of_d_b3),4),
                                     'OF_D_B4':round(np.mean(of_d_b4),4),
                                     'OF_D_B5':round(np.mean(of_d_b5),4),
                                     'OF_D_B7':round(np.mean(of_d_b10),4),
                                     'OF_D_B7_HW':np.std(of_d_b10)*1.96/np.size(of_d_b10),
                                     
                                     'OF_CS':round(np.mean(of_cs),4),
                                     
                                     'ET_ND':round(np.mean(et_nd),4),
                                     'ET_ND_B2':round(np.mean(et_nd_b2),4),
                                     'ET_ND_B3':round(np.mean(et_nd_b3),4),
                                     'ET_ND_B4':round(np.mean(et_nd_b4),4),
                                     'ET_ND_B5':round(np.mean(et_nd_b5),4),
                                     'ET_ND_B7':round(np.mean(et_nd_b10),4),
                                     'ET_ND_B7_HW':np.std(et_nd_b10)*1.96/np.size(et_nd_b10),
                                     
                                     'ET_D':round(np.mean(et_d),4),
                                     'ET_D_B2':round(np.mean(et_d_b2),4),
                                     'ET_D_B3':round(np.mean(et_d_b3),4),
                                     'ET_D_B4':round(np.mean(et_d_b4),4),
                                     'ET_D_B5':round(np.mean(et_d_b5),4),
                                     'ET_D_B7':round(np.mean(et_d_b10),4),
                                     'ET_D_B7_HW':np.std(et_d_b10)*1.96/np.size(et_d_b10),
                                     
                                     'ET_CS':round(np.mean(et_cs),4),

                                     'EFP_ND':round(np.mean(efp_nd),4),
                                     'EFP_ND_B2':round(np.mean(efp_nd_b2),4),
                                     'EFP_ND_B3':round(np.mean(efp_nd_b3),4),
                                     'EFP_ND_B4':round(np.mean(efp_nd_b4),4),
                                     'EFP_ND_B5':round(np.mean(efp_nd_b5),4),
                                     'EFP_ND_B7':round(np.mean(efp_nd_b10),4),
                                     'EFP_ND_B7_HW':np.std(efp_nd_b10)*1.96/np.size(efp_nd_b10),
                                     
                                     'EFP_D':round(np.mean(efp_d),4),
                                     'EFP_D_B2':round(np.mean(efp_d_b2),4),
                                     'EFP_D_B3':round(np.mean(efp_d_b3),4),
                                     'EFP_D_B4':round(np.mean(efp_d_b4),4),
                                     'EFP_D_B5':round(np.mean(efp_d_b5),4),
                                     'EFP_D_B7':round(np.mean(efp_d_b10),4),
                                     'EFP_D_B7_HW':np.std(efp_d_b10)*1.96/np.size(efp_d_b10),
                                     
                                     'EFP_CS':round(np.mean(efp_cs),4),
                                     
                                     'EFN_ND':round(np.mean(efn_nd),4),
                                     'EFN_ND_B2':round(np.mean(efn_nd_b2),4),
                                     'EFN_ND_B3':round(np.mean(efn_nd_b3),4),
                                     'EFN_ND_B4':round(np.mean(efn_nd_b4),4),
                                     'EFN_ND_B5':round(np.mean(efn_nd_b5),4),
                                     'EFN_ND_B7':round(np.mean(efn_nd_b10),4),
                                     'EFN_ND_B7_HW':np.std(efn_nd_b10)*1.96/np.size(efn_nd_b10),
                                     
                                     'EFN_D':round(np.mean(efn_d),4),
                                     'EFN_D_B2':round(np.mean(efn_d_b2),4),
                                     'EFN_D_B3':round(np.mean(efn_d_b3),4),
                                     'EFN_D_B4':round(np.mean(efn_d_b4),4),
                                     'EFN_D_B5':round(np.mean(efn_d_b5),4),
                                     'EFN_D_B7':round(np.mean(efn_d_b10),4),
                                     'EFN_D_B7_HW':np.std(efn_d_b10)*1.96/np.size(efn_d_b10),
                                     
                                     'EFN_CS':round(np.mean(efn_cs),4),
                                     
                                     'PART_ND':round(np.mean(part_nd),4),
                                     'PART_ND_B2':round(np.mean(part_nd_b2),4),
                                     'PART_ND_B3':round(np.mean(part_nd_b3),4),
                                     'PART_ND_B4':round(np.mean(part_nd_b4),4),
                                     'PART_ND_B5':round(np.mean(part_nd_b5),4),
                                     'PART_ND_B10':round(np.mean(part_nd_b10),4),
                                     
                                     'PART_D':round(np.mean(part_d),4),
                                     'PART_D_B2':round(np.mean(part_d_b2),4),
                                     'PART_D_B3':round(np.mean(part_d_b3),4),
                                     'PART_D_B4':round(np.mean(part_d_b4),4),
                                     'PART_D_B5':round(np.mean(part_d_b5),4),
                                     'PART_D_B10':round(np.mean(part_d_b10),4)},ignore_index=True)
           
           # results = results.append({'w1':lambda_1,'w2':lambda_2,'OF_ND':round(np.mean(of_nd),4),'OF_D':round(np.mean(of_d),4),'Change_OF':(round(np.mean(of_d),4)-round(np.mean(of_nd),4))*100/round(np.mean(of_nd),4),'OF_B5':round(np.mean(of_b5),4),'OF_B10':round(np.mean(of_b10),4),'ET_ND':round(np.mean(et_nd),4),'ET_D':round(np.mean(et_d),4),'Change_ET':(round(np.mean(et_d),4)-round(np.mean(et_nd),4))*100/round(np.mean(et_nd),4),'ET_B5':round(np.mean(et_b5),4),'ET_B10':round(np.mean(et_b10),4),'EFP_ND':round(np.mean(efp_nd),4),'EFP_D':round(np.mean(efp_d),4),'Change_EFP':(round(np.mean(efp_d),4)-round(np.mean(efp_nd),4))*100/round(np.mean(efp_nd),4),'EFP_B5':round(np.mean(efp_b5),4),'EFP_B10':round(np.mean(efp_b10),4),'EFN_ND':round(np.mean(efn_nd),4),'EFN_D':round(np.mean(efn_d),4),'Change_EFN':(round(np.mean(efn_d),4)-round(np.mean(efn_nd),4))*100/round(np.mean(efn_nd),4),'EFN_B5':round(np.mean(efn_b5),4),'EFN_B10':round(np.mean(efn_b10),4),},ignore_index=True)
           # results = results.append({'w1':lambda_1,'w2':lambda_2,'OF_ND':round(np.mean(of_nd),4),'OF_D':round(np.mean(of_d),4),'Change_OF':(round(np.mean(of_d),4)-round(np.mean(of_nd),4))*100/round(np.mean(of_nd),4),'ET_ND':round(np.mean(et_nd),4),'ET_D':round(np.mean(et_d),4),'Change_ET':(round(np.mean(et_d),4)-round(np.mean(et_nd),4))*100/round(np.mean(et_nd),4),'EFP_ND':round(np.mean(efp_nd),4),'EFP_D':round(np.mean(efp_d),4),'Change_EFP':(round(np.mean(efp_d),4)-round(np.mean(efp_nd),4))*100/round(np.mean(efp_nd),4),'EFN_ND':round(np.mean(efn_nd),4),'EFN_D':round(np.mean(efn_d),4),'Change_EFN':(round(np.mean(efn_d),4)-round(np.mean(efn_nd),4))*100/round(np.mean(efn_nd),4)},ignore_index=True)
           # #PRINT RESULTS
           # print("OF_ND:\n")
           # interval_of_nd = stats.t.ppf(1.0 - (CI / 2.0),of_nd.size-1) * (np.std(of_nd)/ np.sqrt(of_nd.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(of_nd),4))+"+-"+str(round(interval_of_nd,4)))
           
           # print("\nET_ND:\n")
           # interval_et_nd = stats.t.ppf(1.0 - (CI / 2.0),et_nd.size-1) * (np.std(et_nd)/ np.sqrt(et_nd.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(et_nd),4))+"+-"+str(round(interval_et_nd,4)))
           
           # print("\nEFP_ND:\n")
           # interval_efp_nd = stats.t.ppf(1.0 - (CI / 2.0),efp_nd.size-1) * (np.std(efp_nd)/ np.sqrt(efp_nd.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(efp_nd),4))+"+-"+str(round(interval_efp_nd,4)))
          
           # print("\nEFN_ND:\n")
           # interval_efn_nd = stats.t.ppf(1.0 - (CI / 2.0),efn_nd.size-1) * (np.std(efn_nd)/ np.sqrt(efn_nd.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(efn_nd),4))+"+-"+str(round(interval_efn_nd,4)))
           
           # print('\n')
           # print('--------------------------------')
           # print('\n')
           # print("OF_D:\n")
           # interval_of_d = stats.t.ppf(1.0 - (CI / 2.0),of_d.size-1) * (np.std(of_d)/ np.sqrt(of_d.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(of_d),4))+"+-"+str(round(interval_of_d,4)))
           
           # print("\nET_D:\n")
           # interval_et_d = stats.t.ppf(1.0 - (CI / 2.0),et_d.size-1) * (np.std(et_d)/ np.sqrt(et_d.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(et_d),4))+"+-"+str(round(interval_et_d,4)))
           
           # print("\nEFP_D:\n")
           # interval_efp_d = stats.t.ppf(1.0 - (CI / 2.0),efp_d.size-1) * (np.std(efp_d)/ np.sqrt(efp_d.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(efp_d),4))+"+-"+str(round(interval_efp_d,4)))
           
           # print("\nEFN_D:\n")
           # interval_efn_d = stats.t.ppf(1.0 - (CI / 2.0),efn_d.size-1) * (np.std(efn_d)/ np.sqrt(efn_d.size))
           # # ci = (a.mean() - interval, a.mean() + interval)
           # print(str(round(np.mean(efn_d),4))+"+-"+str(round(interval_efn_d,4)))
           
           
# sys.stdout.close()
#%% 8. OUTPUT TO CSV
# of_d.to_csv('OF_D.csv')
# of_cs.to_csv('OF_CS.csv')
results.to_csv('qwertyuiop.csv')