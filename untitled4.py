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

#%% 1. IMPORT LIBRARIES
# import matplotlib.pyplot as plt
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

  risk_vec.sort()
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

def Bellman(B):   
    d=[[None for i in range(pop+1)] for i in range(B+1)]
    d[0]=[float("inf") for i in range(pop+1)]
    d[0][0]=0.0
    pathNode=[[None for i in range(pop+1)] for i in range(B+1)]
    for k in range(B):
        for i in range(pop+1):
            tempMin=float("inf")
            tempIndex=-1.0
            for j in range(i):
                if d[k][j]+wt_D[j][i]<tempMin:
                    tempMin=d[k][j]+wt_D[j][i]
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

    # of=0
    # of+=weight_D(node_list[0],pop)
    et=0
    et+=Exp_Test_Num_D(node_list[0],pop)
    efp=0
    efp+=Exp_FalsePos_D(node_list[0],pop)
    efn=0
    efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        # of+=weight_D(node_list[ele+1],node_list[ele])
        et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])  
    
    print(node_list)  
    print(d[B][pop])
    # return d[B][pop],et,efp,efn,len(node_list)


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

#Confidence Interval
CI = 0.05

#Simulations
simulations=1





#%% 7. SIMULATIONS
# sys.stdout = open("FINALRESULTS.txt", "w")
results = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'OF_ND',
                                 'OF_D',
                                 'OF_B5',
                                 'OF_B10',
                                 'ET_ND',
                                 'ET_D',
                                 'ET_B5',
                                 'ET_B10',
                                 'EFP_ND',
                                 'EFP_D',
                                 'EFP_B5',
                                 'EFP_B10',
                                 'EFN_ND',
                                 'EFN_D',
                                 'EFN_B5',
                                 'EFN_B10',
                                 'PART_ND',
                                 'PART_D',
                                 'PART_B5',
                                 'PART_B10']) 
           
for w1 in tqdm(range(0, 12, 2),desc='lambda loop',position=0,leave=True):
   for w2 in range(0,12,2):
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
           
           of_d=np.empty(simulations)
           et_d=np.empty(simulations)
           efp_d=np.empty(simulations)
           efn_d=np.empty(simulations)
           part_d=np.empty(simulations)
           
           of_b5=np.empty(simulations)
           et_b5=np.empty(simulations)
           efp_b5=np.empty(simulations)
           efn_b5=np.empty(simulations)
           part_b5=np.empty(simulations)
           
           of_b10=np.empty(simulations)
           et_b10=np.empty(simulations)
           efp_b10=np.empty(simulations)
           efn_b10=np.empty(simulations)
           part_b10=np.empty(simulations)
           for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
               get_risk_vector()
               
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
               
               of_b5[sim],et_b5[sim],efp_b5[sim],efn_b5[sim],part_b5[sim]=Bellman(5)
               of_b10[sim],et_b10[sim],efp_b10[sim],efn_b10[sim],part_b10[sim]=Bellman(10)

           results = results.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'OF_ND':round(np.mean(of_nd),4),
                                     'OF_D':round(np.mean(of_d),4),
                                     'OF_B5':round(np.mean(of_b5),4),
                                     'OF_B10':round(np.mean(of_b10),4),
                                     'ET_ND':round(np.mean(et_nd),4),
                                     'ET_D':round(np.mean(et_d),4),
                                     'ET_B5':round(np.mean(et_b5),4),
                                     'ET_B10':round(np.mean(et_b10),4),
                                     'EFP_ND':round(np.mean(efp_nd),4),
                                     'EFP_D':round(np.mean(efp_d),4),
                                     'EFP_B5':round(np.mean(efp_b5),4),
                                     'EFP_B10':round(np.mean(efp_b10),4),
                                     'EFN_ND':round(np.mean(efn_nd),4),
                                     'EFN_D':round(np.mean(efn_d),4),
                                     'EFN_B5':round(np.mean(efn_b5),4),
                                     'EFN_B10':round(np.mean(efn_b10),4),
                                     'PART_ND':round(np.mean(part_nd),4),
                                     "PART_D":round(np.mean(part_d),4),
                                     "PART_B5":round(np.mean(part_b5),4),
                                     "PART_B10":round(np.mean(part_b10),4)},ignore_index=True)
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

results.to_csv('Bellman1SimulationCorrect_Modified.csv')