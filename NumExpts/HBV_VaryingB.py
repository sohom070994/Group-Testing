# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%% 1. IMPORT LIBRARIES
import seaborn as sns
# from matplotlib.lines import Line2D
import math
import time 
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
import numpy as np 
# import scipy.stats as stats
# import sys
import pandas as pd
#%% 2. RISK VECTOR REALIZATION
'''
LIST OF FUNCTIONS 
FOR GENERATING RISK VECTOR
'''


def coin_flip():

    # randnum = round(random.uniform(0,100),2)
    # if 0.00 <= randnum <= 16.34:
    #   return 0.09*3
    # elif randnum <= 33.45:
    #   return 0.48*3
    # elif randnum <= 48.89:
    #   return 1.07*3
    # elif randnum <= 64.91:
    #   return 0.81*3
    # elif randnum <= 80.78:
    #   return 0.42*3
    # else:
    #   return 0.18*3
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

  # risk_vec.sort()
#%% 3. Debugging Block
# get_risk_vector()
# print(risk_vec)
# print(len(risk_vec))
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

def Exp_FalseNeg_ND_Mod(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return (1-const_sens)*risk_vec[Node_1]
  else:
    wt = sum(risk_vec[Node_1:Node_2])
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

# def Exp_FalsePos_ND_Mod(Node_1,Node_2):
#   if Node_2==Node_1+1: 
#     return (1-spec)*(1-risk_vec[Node_1])
#   else:
#     summ = Node_2-Node_1-sum(risk_vec[Node_1:Node_2])
#     wt = math.prod([(1-risk_vec[ele]) for ele in range(Node_1,Node_2)])
#     return ((1-spec)*const_sens*sum)-((Node_2-Node_1)*(1-spec)*(const_sens+spec-1)*(wt))


def Exp_Test_Num_ND(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return 1
  else:
    wt = 1
    for ele in range(Node_1,Node_2):
      wt = wt*(1-risk_vec[ele])
    return 1+(Node_2-Node_1)*(const_sens-(const_sens+spec-1)*(wt))
 
# def Exp_Test_Num_ND_Mod(Node_1,Node_2): 
#   if Node_2==Node_1+1: 
#     return 1
#   else:
#     wt = math.prod([1-risk_vec[ele] for ele in range(Node_1,Node_2)])
#     return 1+(Node_2-Node_1)*(const_sens-(const_sens+spec-1)*(wt))    

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_ND_old(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_ND(Node_1,Node_2))+(lambda_2*Exp_FalsePos_ND(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_ND(Node_1,Node_2))
  return weight

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

# def Exp_FalseNeg_D_Mod(Node_1,Node_2): 
#     if Node_2==Node_1+1: 
#        return (1-sens(1,1))*risk_vec[Node_1]
#     else:
#         sum_1 = sum(risk_vec[Node_1:Node_2])
#         sum_2 = sum([((risk_vec[ele])/(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])) for ele in range(Node_1,Node_2)]) 
#         wt = math.prod([(1-alpha_pr(Node_1,Node_2)*risk_vec[ele]) for ele in range(Node_1,Node_2)])
#         return  (1-sens(1,1))*sum_1+sens(1,1)*spec*alpha**(1/((Node_2-Node_1)**gamma))*wt*sum_2


    
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
    
    
# def Exp_FalsePos_D_Mod(Node_1,Node_2):
#     if Node_2==Node_1+1: 
#         return (1-spec)*(1-risk_vec[Node_1])
#     else:
#         sum_1 = Node_2-Node_1-sum(risk_vec[Node_1:Node_2]) 
#         sum_2 = sum([((1-risk_vec[ele])/(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])) for ele in range(Node_1,Node_2)])
#         wt = math.prod([(1-alpha_pr(Node_1,Node_2)*risk_vec[ele]) for ele in range(Node_1,Node_2)])
#         return ((1-spec)*sum_1)-((1-spec)*spec*wt*sum_2)    

def Exp_Test_Num_D(Node_1,Node_2): 
    if Node_2==Node_1+1: 
        return 1
    else:
        wt = 1
        for ele in range(Node_1,Node_2):
          wt = wt*(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])
        return 1+(Node_2-Node_1)*(1-spec*(wt))

# def Exp_Test_Num_D_Mod(Node_1,Node_2): 
#     if Node_2==Node_1+1: 
#       return 1
#     else:
#         wt = math.prod([(1-alpha_pr(Node_1,Node_2)*risk_vec[ele]) for ele in range(Node_1,Node_2)])
#         return 1+(Node_2-Node_1)*(1-spec*(wt))
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight_D(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_D(Node_1,Node_2))+(lambda_2*Exp_FalsePos_D(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_D(Node_1,Node_2))
  return weight

def weight_D_old(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg_D(Node_1,Node_2))+(lambda_2*Exp_FalsePos_D(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num_D_Mod(Node_1,Node_2))
  return weight
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
    # return of,et,efp,efn,len(node_list)
    return of,len(node_list)

def Bellmanmod(m,cost,B):
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
risk_vec = [0] * pop

# #Defining weights for Objective Function
# lambda_1 = 1.0 #weight on EFN
# lambda_2 = 0.0 #weight on EFP

#Parameters
spec = 0.99
alpha= 0.0514
gamma= 0.2416
const_sens = sens(1,1)

# #Confidence Interval
# CI = 0.05

#Simulations
simulations=100
# #%% ONE SIM
# start_time=time.time()
# lambda_1=0
# lambda_2=0
# get_risk_vector()
# g_d = Graph(pop+1)
# risk_vec.sort()
# for i in range(0,pop):
#     for j in range(i+1,pop+1):
#         g_d.addEdge(i,j,weight_D(i,j))
# print("--- %s seconds ---" % (time.time() - start_time))
# node_list=g_d.shortestPath(s) 
# print("--- %s seconds ---" % (time.time() - start_time))
w1=10
w2=0




               



# wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
# for i in range(0,pop):
#     for j in range(i+1,pop+1):
#           wt_D[i][j]=weight_D(i,j)

# for i in range(0,pop):
#     for j in range(i+1,pop+1):
#         g_d.addEdge(i,j,wt_D[i][j])
# print("--- %s seconds ---" % (time.time() - start_time))
# of_md_b5,et_md_b5,efp_md_b5,efn_md_b5,part_md_b5=Bellmanmod(21,wt_D,15)
# of_d_b5,et_d_b5,efp_d_b5,efn_d_b5,part_d_b5=Bellman(wt_D,15)
#%% 7. SIMULATIONS
# sys.stdout = open("FINALRESULTS.txt", "w")
results = pd.DataFrame(columns= ['B','OF_ND_B','OF_D_B','OF_CP','PART_ND_B','PART_D_B']) 


for b in tqdm([7,8,9],desc='lambda loop',position=0,leave=True):
    pop = 100
    lambda_1=w1/10
    lambda_2=w2/10
    print('\n')
    print('------------------------------------------------------------------')
    print('B :'+str(b)+'\n')
    print('------------------------------------------------------------------')
    print('\n')
    of_nd_b=np.empty(simulations)
    part_nd_b=np.empty(simulations)
    of_d_b=np.empty(simulations)
    part_d_b=np.empty(simulations)
    of_cs=np.empty(simulations)
    for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
        print('SIM :'+str(sim)+'\n')
        get_risk_vector()
        
        #Current Scenario
        current_node_list=[96,80,64,49,32,16,0]
        of_cs[sim]=calc_meas(current_node_list)[0]
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
        
        # node_list=g_nd.shortestPath(s) 
        # part_nd_b[sim]=len(node_list)
        # of_nd_b[sim]=calc_meas(node_list)[0]
        
        g_d = Graph(pop+1) 
   
        for i in range(0,pop):
            for j in range(i+1,pop+1):
                g_d.addEdge(i,j,wt_D[i][j])
         
        # node_list=g_d.shortestPath(s) 
        # part_d_b[sim]=len(node_list)
        # of_d_b[sim]=calc_meas(node_list)[0]
        
        of_nd_b[sim],part_nd_b[sim]=Bellman(wt_ND,b)
        of_d_b[sim],part_d_b[sim]=Bellman(wt_D,b)
    
    results = results.append({'B':b,
                                 'OF_ND_B':round(np.mean(of_nd_b),4),
                                 
                                 'OF_D_B':round(np.mean(of_d_b),4),
                                 
                                 
                                 'OF_CP':round(np.mean(of_cs),4),
                                 
                                 
                                 'PART_ND_B':round(np.mean(part_nd_b),4),
                                 'PART_D_B':round(np.mean(part_d_b),4),},ignore_index=True)
           
# for w1 in tqdm(list(range(0, 12, 2))+[9.5],desc='lambda loop',position=0,leave=True):
#    for w2 in list(range(0,12,2))+[0.25]:
#        # if(w1!=9.6)&(w2==0.2):
#        #     continue
#        if (w1==9.5)&(w2!=0.25):
#             continue
#        # if not (((w1==10)&(w2==0))|((w1==0)&(w2==10))):
#        #      continue
#        if(w1+w2<=10):
#            lambda_1=w1/10
#            lambda_2=w2/10
           
           


#            # print('\n')
#            # print('------------------------------------------------------------------')
#            print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
#            # print('------------------------------------------------------------------')
#            # print('\n')
        
#            of_nd=np.empty(simulations)
#            et_nd=np.empty(simulations)
#            efp_nd=np.empty(simulations)
#            efn_nd=np.empty(simulations)
#            part_nd=np.empty(simulations)
           
#            of_nd_b7=np.empty(simulations)
#            et_nd_b7=np.empty(simulations)
#            efp_nd_b7=np.empty(simulations)
#            efn_nd_b7=np.empty(simulations)
#            part_nd_b7=np.empty(simulations)          
           
#            of_d=np.empty(simulations)
#            et_d=np.empty(simulations)
#            efp_d=np.empty(simulations)
#            efn_d=np.empty(simulations)
#            part_d=np.empty(simulations)
           
#            of_d_b7=np.empty(simulations)
#            et_d_b7=np.empty(simulations)
#            efp_d_b7=np.empty(simulations)
#            efn_d_b7=np.empty(simulations)
#            part_d_b7=np.empty(simulations)
           
#            of_cs=np.empty(simulations)
#            et_cs=np.empty(simulations)
#            efp_cs=np.empty(simulations)
#            efn_cs=np.empty(simulations)
#            # part_cs=np.empty(simulations)
           
#            for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
#                get_risk_vector()
               
#                #Current Scenario
#                current_node_list=[96,80,64,49,32,16,0]
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
               
#                of_nd_b7[sim],et_nd_b7[sim],efp_nd_b7[sim],efn_nd_b7[sim],part_nd_b7[sim]=Bellman(wt_ND,8)

#                of_d_b7[sim],et_d_b7[sim],efp_d_b7[sim],efn_d_b7[sim],part_d_b7[sim]=Bellman(wt_D,8)

#            results = results.append({'w1':lambda_1,
#                                      'w2':lambda_2,
#                                      'OF_ND':round(np.mean(of_nd),4),
#                                      'OF_ND_HW':1.96*np.std(of_nd)/np.size(of_nd),
#                                      'OF_ND_B7':round(np.mean(of_nd_b7),4),                                     
#                                      'OF_ND_B7_HW':1.96*np.std(of_nd_b7)/np.size(of_nd_b7),
                                     
#                                      'OF_D':round(np.mean(of_d),4),
#                                      'OF_D_HW':1.96*np.std(of_d)/np.size(of_d),
#                                      'OF_D_B7':round(np.mean(of_d_b7),4),
#                                      'OF_D_B7_HW':1.96*np.std(of_d_b7)/np.size(of_d_b7),
                                     
                                     
#                                      'OF_CS':round(np.mean(of_cs),4),
#                                      'OF_CS_HW':1.96*np.std(of_cs)/np.size(of_cs),
                                     
                                     
#                                      'ET_ND':round(np.mean(et_nd),4),
#                                      'ET_ND_B7':round(np.mean(et_nd_b7),4),
                                     
#                                      'ET_D':round(np.mean(et_d),4),
#                                      'ET_D_B7':round(np.mean(et_d_b7),4),
                                     
#                                      'ET_CS':round(np.mean(et_cs),4),

#                                      'EFP_ND':round(np.mean(efp_nd),4),
#                                      'EFP_ND_B7':round(np.mean(efp_nd_b7),4),
                                     
#                                      'EFP_D':round(np.mean(efp_d),4),
#                                      'EFP_D_B7':round(np.mean(efp_d_b7),4),
                                     
#                                      'EFP_CS':round(np.mean(efp_cs),4),
                                     
#                                      'EFN_ND':round(np.mean(efn_nd),4),
#                                      'EFN_ND_B7':round(np.mean(efn_nd_b7),4),
                                     
#                                      'EFN_D':round(np.mean(efn_d),4),
#                                      'EFN_D_B7':round(np.mean(efn_d_b7),4),
                                     
#                                      'EFN_CS':round(np.mean(efn_cs),4),
                                     
#                                      'PART_ND':round(np.mean(part_nd),4),
#                                      'PART_ND_B7':round(np.mean(part_nd_b7),4),
                                     
#                                      'PART_D':round(np.mean(part_d),4),
#                                      'PART_D_B7':round(np.mean(part_d_b7),4),},ignore_index=True)
           
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

# pd.DataFrame(of_d).to_csv("OF_D.csv")
# pd.DataFrame(of_cs).to_csv("OF_CS.csv")
results.to_csv('HBVTable_Bvar789.csv')

#%% 9. PLOT HISTOGRAM
sns.distplot(of_d, hist=True, kde=True, 
             bins=20, color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4},label='Proposed')

mean_d=of_d.mean()
mean_cs=of_cs.mean()
plt.axvline(mean_d, c= 'blue', ls= '--', label='Mean (Proposed)')
sns.distplot(of_cs, hist=True, kde=True, 
             bins=20, color = 'darkred', 
             hist_kws={'edgecolor':'red'},
             kde_kws={'linewidth': 4}, label='Current',axlabel='Objective Function')

# ax_hist.axvline(mean_d, color='blue', linestyle='--')   
# ax_hist.axvline(mean_cs, color='r', linestyle='--')

plt.ylabel("Frequency")
plt.axvline(mean_cs, c= 'red', ls='--', label='Mean (Current)')

plt.xlim(0.29, 0.48)
plt.ylim(0, 50)
plt.title("Current Scheme vs Proposed Scheme, "+r'$\lambda_1$'+"="+str(lambda_1)+", "+r'$\lambda_2$'+"="+str(lambda_2)) 

plt.legend()
plt.show()
# sns.plt.show()
#%% 10. Calibration
def sen(al,gm,x,sp,n):
  '''
  Input- 
    n(int): number of people in group
    k(int): number of infected people

  Output- 
    returns sensitivity(float) from dilution equation
  '''
  return 1-sp*al**(x/(n**gm))

# al=1
# gm=1
# sp=0.8
# n=10
a=list(range(0,100001,1))
x=[i/10000 for i in a]
y1=[sen(0.9,0.5,a,0.95,10) for a in x ]
y2=[sen(0.2,1.5,a,0.95,10) for a in x ]
y3=[sen(0.05,0.9,a,0.95,10) for a in x ]
y4=[sen(0.001,0.2,a,0.95,10) for a in x ]
y5=[sen(0.8,0.2,a,0.95,10) for a in x ]
plt.xlim(-0.1,10.1)
plt.ylim(-0.1, 1.2)
plt.scatter(x, y1,s=1,marker=',')
plt.scatter(x, y2,s=1,marker=',')
plt.scatter(x, y3,s=1,marker=',')
plt.scatter(x, y4,s=1,marker=',')
plt.scatter(x, y5,s=1,marker=',')
plt.xlabel("Infected individuals")
plt.ylabel("Sensitivity")
# legend_elements = [Line2D([0], [0], color='b', lw=4, label='Line')]
# ax.legend(handles=legend_elements)


plt.axhline(1.0, c= 'red', ls='--')
plt.title("Sensitivity functions for group size, n= 10") 

# plt.legend()
# ax = sns.regplot(x=x, y=y1,
#                  # scatter_kws={"s": 80},
#                  order=100, ci=None, scatter=False)

plt.show()