# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%% 
# 1. IMPORT LIBRARIES
import matplotlib.pyplot as plt
import seaborn as sns
import random
from sympy import bell
from tqdm import tqdm
import numpy as np 
# import scipy.stats as stats
# import sys
import pandas as pd
import time
#%% 
# 2. RISK VECTOR REALIZATION
'''
LIST OF FUNCTIONS 
FOR GENERATING RISK VECTOR
'''

def coin_flip(p):
    toss=random.uniform(0, 1)
    if toss>p:
        return 0
    else:
        return 1
    

def get_status_vec(risk_vec):
    return [coin_flip(i) for i in risk_vec]    
    
# COVID19 risk distribution   
def get_risk_vector(k):
    # '''
    # fill risk_vec
    # using coin flips
    # sort risk_vec
    # '''
    # for i in range(0,pop):
    #   risk_vec[i]=round(coin_flip()/100,4)
    # elements = [0.198,0.138,0.46,0.395,0.126]
    elements = [0.0198,0.0138,0.046,0.0395,0.0126]
    elements = [i*k for i in elements]
    probabilities = [0.025,0.3715,0.2685,0.0031,0.3319]
    a=np.random.choice(elements, pop, p=probabilities)
    return a

# HBV risk distribution
# def get_risk_vector():
#     elements = [0.0009,0.0048,0.0107,0.0081,0.0042,0.0018]
#     probabilities = [0.16339801983321,0.171090634447457,0.15443199329378,0.160180315971968,0.158686396851769,0.192212639601815]
#     a=np.random.choice(elements, pop, p=probabilities)
#     return a


#%% 
# 3. Debugging Block
# abc=get_risk_vector()
# print(abc)
# print(len(abc))
#%% 
# 3B. Dorfman Simulation Block
def partition(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset 
        yield [ [ first ] ] + smaller
        
        
def DorfmanSim(rv,part):
    sv=get_status_vec(rv)
    # print(sv)
    exp_tests=0
    exp_fn=0
    exp_fp=0
    for grp in part:
        k=0
        for i in grp:
            k+=sv[i]
        sens=gen_sens(len(grp), k)
        gt_result=int(np.random.choice([0,1],1,p=[1-sens,sens]))
        if gt_result==0:
            exp_tests+=1
            exp_fn+=k
        if gt_result==1:
            if len(grp)==1:
                 exp_tests+=1
            else:     
                exp_tests+=1+len(grp)
            for i in grp: 
                it_result=int(np.random.choice([0,1],1,p=[1-gen_sens(1,sv[i]),gen_sens(1,sv[i])]))
                if (it_result-sv[i])==-1:
                    exp_fn+=1
                if (it_result-sv[i])==1:
                    exp_fp+=1
    return exp_tests,exp_fn,exp_fp,lambda_1*exp_fn + lambda_2*exp_fp + (1-lambda_1-lambda_2)*exp_tests
# #%%
# simulations=1000
# t=np.empty(simulations)
# fn=np.empty(simulations)
# fp=np.empty(simulations)
# for i in range(simulations):
#     rv = np.array([1 for i in range(7)])
#     t[i],fn[i],fp[i]= DorfmanSim(rv,[[0,1,2,3,4,5,6]])

#%% 
# 4. DILUTION CASE PERFORMANCE MEASURES

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
# 4.c U-RDDS PERFORMANCE MEASURES
'''
LIST OF 
PERFORMANCE MEASURES
'''
def risk_sum(Node_1,Node_2):
    s = 0 
    for ele in range(Node_1,Node_2):
      s = s + risk_vec[ele]
    return s  

def gen_sens(n,k):  
    return sens(n,k)

def gen_sensBM(n,k):
  '''
  Input- 
    n(int): number of people in group
    k(int): number of infected people

  Output- 
    returns sensitivity(float) from dilution equation
  '''
  # return 1-spec*alpha**(k/(n**gamma))
  # return -1.36+2.359*((k/n)**0.02)
  # return 0.984*(k/n)**0.0738
  # return 0.993*(k/n)**0.08888
  return 1.00*(k/n)**0.0507

def Exp_FalseNeg_UB(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return (1-gen_sens(1, 1))*risk_vec[Node_1]
  else:
    sum = risk_sum(Node_1,Node_2)
    return gen_sens(1, 1)*(1-gen_sens(Node_2-Node_1, sum))*sum
    
def Exp_FalsePos_UB(Node_1,Node_2):
  if Node_2==Node_1+1: 
    return (1-spec)*(1-risk_vec[Node_1])
  else:
    sum = risk_sum(Node_1,Node_2)
    return (1-spec)*(Node_2-Node_1-sum)*gen_sens(Node_2-Node_1, sum)

def Exp_Test_Num_UB(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return 1
  else:
    sum = risk_sum(Node_1,Node_2)
    return 1+(Node_2-Node_1)*gen_sens(Node_2-Node_1, sum)
    
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_UB(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_UB(Node_1,Node_2))+(lambda_2*Exp_FalsePos_UB(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_UB(Node_1,Node_2))
  return weight

def calc_measub(node_list):

    of=0
    of+=weight_UB(node_list[0],pop)
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
        of+=weight_UB(node_list[ele+1],node_list[ele])
    return of


def get_urdds_part():
    wt_UB =[[None for i in range(pop+1)]for i in range(pop+1)]
    for i in range(0,pop):
        for j in range(i+1,pop+1):
              wt_UB[i][j]=weight_UB(i,j)              
    s = 0
    g_ub = Graph(pop+1)  
    for i in range(0,pop):
        for j in range(i+1,pop+1):
            g_ub.addEdge(i, j, wt_UB[i][j])
    node_list_ub=g_ub.shortestPath(s)
    return node_list_ub,calc_measub(node_list_ub)

def modify_part(part):
    nl_rev=part[::-1]
    output=[]
    for i in range(0, len(nl_rev)):
        if (i==len(nl_rev)-1):
            output.append(list(range(nl_rev[i],pop)))
        else:
            output.append(list(range(nl_rev[i],nl_rev[i+1])))    
    # for ele in range(0,len(output)):
    #     for pat in range(0,len(output[ele])): 
    #         output[ele][pat]=perm_vec[output[ele][pat]]
    # output.sort()        
    return output     

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


def get_optimal_part():
    wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
    for i in range(0,pop):
        for j in range(i+1,pop+1):
              wt_D[i][j]=weight_D(i,j)              
    s = 0
    g_d = Graph(pop+1)  
    for i in range(0,pop):
        for j in range(i+1,pop+1):
            g_d.addEdge(i, j, wt_D[i][j])
    node_list_d=g_d.shortestPath(s)
    return node_list_d,calc_measof(node_list_d)


#%% 
# 6. ESTABLISH CONSTANT PARAMETERS
#population
pop = 100

# COVID Parameters
spec = 1
alpha= 0.012595
gamma= 0.257138
const_sens = sens(1,1)

# # HBV Parameters
# spec = 0.99
# alpha= 0.0514
# gamma= 0.2416
# const_sens = sens(1,1)

#Confidence Interval
CI = 0.05

#Simulations
simulations=100
DSim=100
lambda_1=0
lambda_2=0
#%% 
# 7 SIMULATIONS- Special sensitivity function
# ET_CP=np.empty(simulations)
ub_part=[]
opti_part=[]
ET_URDDS=np.empty(simulations)
opti=np.empty(simulations)
for sim in tqdm(range(0,simulations)):
    uns_risk_vec=get_risk_vector(1)
    risk_vec=np.sort(uns_risk_vec)
    urdds_part=get_urdds_part()[0]
    ub_part.append(urdds_part)
    ET_URDDS_UB=get_urdds_part()[1]
    opti_part.append(get_optimal_part()[0])
    opti[sim]=get_optimal_part()[1]
    ET_URDDS[sim]=calc_measof(urdds_part)

    # print('---')   
    
print('Percentage Difference: '+str(round((np.mean(ET_URDDS)-np.mean(opti))/np.mean(opti)*100,3))+'%' )    
    
#%% 
# 7 SIMULATIONS- Special sensitivity function
# ET_CP=np.empty(simulations)
ub_part=[]
opti_part=[]
ET_URDDS=np.empty(simulations)
opti=np.empty(simulations)
for sim in tqdm(range(0,simulations)):
    uns_risk_vec=get_risk_vector(1)
    risk_vec=np.sort(uns_risk_vec)
    urdds_part=get_urdds_part()[0]
    ub_part.append(urdds_part)
    ET_URDDS_UB=get_urdds_part()[1]
    opti_part.append(get_optimal_part()[0])
    opti[sim]=get_optimal_part()[1]
    ET_URDDS[sim]=calc_measof(urdds_part)

    # print('---')   
    
print('Percentage Difference: '+str(round((np.mean(ET_URDDS)-np.mean(opti))/np.mean(opti)*100,3))+'%' )    
    
#%% 
# 7B. PLOT HISTOGRAM
sns.distplot(opti, hist=True, kde=False, 
             bins=5, color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='GO')

mean_opti=np.mean(opti)
mean_ub=np.mean(ET_URDDS)
plt.axvline(mean_opti, c= 'blue', ls= '--', label='Mean (GO)')
sns.distplot(ET_URDDS, hist=True, kde=False, 
             bins=5, color = 'darkred', 
             hist_kws={'edgecolor':'red'},
             kde_kws={'linewidth': 4}, label='URDDS',axlabel='ET')

# ax_hist.axvline(mean_d, color='blue', linestyle='--')   
# ax_hist.axvline(mean_cs, color='r', linestyle='--')

plt.ylabel("Frequency")
plt.axvline(mean_ub, c= 'red', ls='--', label='Mean (URDDS)')

# plt.xlim(0.29, 0.48)
# plt.ylim(0, 50)
plt.title('Percentage Difference: '+str(round((np.mean(ET_URDDS)-np.mean(opti))/np.mean(opti)*100,3))+'%' ) 

plt.legend()
plt.show()
# sns.plt.show()

#%% 7.a SIMULATIONS- wo perm vec
# ET_CP=np.empty(simulations)
# def game(km):
ub_part=[]
opti_part=[]
risks=[]
ET_URDDS=np.empty(simulations)
opti=np.empty(simulations)
# for w1 in tqdm(list(range(0, 120, 25)),desc='lambda loop',position=0,leave=True):
for w1 in range(0, 120, 25):  
    lambda_1=w1/100
    # for sim in tqdm(list(range(0,simulations)),position=0,leave=True):
    for sim in range(0,simulations):
        uns_risk_vec=get_risk_vector(1)
        risk_vec=np.sort(uns_risk_vec)
        risks.append(risk_vec) 
        urdds_part=get_urdds_part()[0]
        ET_URDDS_UB=get_urdds_part()[1]
        urdds_part_mod=modify_part(urdds_part)
        ub_part.append(urdds_part)
        for ele in range(0,len(urdds_part_mod)):
            urdds_part_mod[ele].sort()
        ET_allpart=np.empty((bell(pop)))
        parts=[]

        running_part=0
        for n, p in (enumerate(partition(list(range(0,pop))),1)):
            # print(n,p)
            p.sort()
            # if (n%1000)==0:
            #     print (n)
            parts= np.append(parts,str(p))
            ET=np.empty(DSim)
            for dg in range(0,DSim):
                ET[dg]=DorfmanSim(risk_vec,p)[3]
            ET_allpart[n-1]=round(np.mean(ET),4) 
            # if n==1:
            #     running_part=p
            # else: 
            #     if(ET_allpart[n-1]<=ET_allpart[n-2]):
            #         running_part=p
            # print(n,p,ET_allpart[n-1])
            if (urdds_part_mod==p):
                ET_URDDS[sim]=round(np.mean(ET),4)
                # print(urdds_part_mod)
                # print(p)
                # print('URDDS_OF: '+str(round(np.mean(ET),4) ))
        opti[sim]=np.min(ET_allpart)
        opti_part.append(p)
        # print(gen_sensBM(5,sum(risk_vec)))
        # print('riskvec'+str(risk_vec))
        # print('GO'+str(opti[sim]))
        # print('uRDDS'+str(ET_URDDS[sim]))
        
        # print(sim)
        # print('---')    
    km=1            
    print('MF:'+str(km)+'. '+'L1: '+ str(lambda_1)+'. Percentage Difference: '+str(round((np.mean(ET_URDDS)-np.mean(opti))/np.mean(opti)*100,3))+'%' )    


# for i in range(1,6):
#   game(i)
# %%

# %%
