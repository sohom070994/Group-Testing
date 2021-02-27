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
simulations=500

#%% 8. VARYING B SIMULATIONS
lambda_1=1
lambda_2=0
simulations = 400
results = pd.DataFrame(columns= ['OF_ND',
                                 'OF_ND_B5',
                                 'OF_ND_B10',
                                 'OF_ND_B15',
                                 'OF_ND_B16',    
                                 'OF_ND_B17',    
                                 'OF_ND_B18',   
                                 'OF_ND_B19',                                 
                                 'OF_ND_B20',
                                 'OF_ND_B30',
                                 'OF_ND_B40',
                                 'OF_ND_B50',
                                 'OF_ND_B60',
                                 'OF_ND_B75',
                                 'PART_ND',
                                 'PART_ND_B5',
                                 'PART_ND_B10',
                                 'PART_ND_B15',
                                 'PART_ND_B16',    
                                 'PART_ND_B17',    
                                 'PART_ND_B18',   
                                 'PART_ND_B19',                                 
                                 'PART_ND_B20',
                                 'PART_ND_B30',
                                 'PART_ND_B40',
                                 'PART_ND_B50',
                                 'PART_ND_B60',
                                 'PART_ND_B75',
                                 'OF_D',
                                 'OF_D_B5',
                                 'OF_D_B10',
                                 'OF_D_B15',
                                 'OF_D_B16',    
                                 'OF_D_B17',    
                                 'OF_D_B18',   
                                 'OF_D_B19',                                 
                                 'OF_D_B20',
                                 'OF_D_B30',
                                 'OF_D_B40',
                                 'OF_D_B50',
                                 'OF_D_B60',
                                 'OF_D_B75',
                                 'PART_D',
                                 'PART_D_B5',
                                 'PART_D_B10',
                                 'PART_D_B15',
                                 'PART_D_B16',    
                                 'PART_D_B17',    
                                 'PART_D_B18',   
                                 'PART_D_B19',
                                 'PART_D_B20',
                                 'PART_D_B30',
                                 'PART_D_B40',
                                 'PART_D_B50',
                                 'PART_D_B60',
                                 'PART_D_B75',
                                 'DIFF',
                                 'DIFF_B5',
                                 'DIFF_B10',
                                 'DIFF_B15',
                                 'DIFF_B16',    
                                 'DIFF_B17',    
                                 'DIFF_B18',   
                                 'DIFF_B19',
                                 'DIFF_B20',
                                 'DIFF_B30',
                                 'DIFF_B40',
                                 'DIFF_B50',
                                 'DIFF_B60',
                                 'DIFF_B75'
                                 ])

of_nd=np.empty(simulations)
of_nd_b5=np.empty(simulations)
of_nd_b10=np.empty(simulations)
of_nd_b15=np.empty(simulations)
of_nd_b16=np.empty(simulations)
of_nd_b17=np.empty(simulations)
of_nd_b18=np.empty(simulations)
of_nd_b19=np.empty(simulations)
of_nd_b20=np.empty(simulations)
of_nd_b30=np.empty(simulations)
of_nd_b40=np.empty(simulations)
of_nd_b50=np.empty(simulations)
of_nd_b60=np.empty(simulations)
of_nd_b75=np.empty(simulations)

part_d=np.empty(simulations)
part_d_b5=np.empty(simulations)
part_d_b10=np.empty(simulations)
part_d_b15=np.empty(simulations)
part_d_b16=np.empty(simulations)
part_d_b17=np.empty(simulations)
part_d_b18=np.empty(simulations)
part_d_b19=np.empty(simulations)
part_d_b20=np.empty(simulations)
part_d_b30=np.empty(simulations)
part_d_b40=np.empty(simulations)
part_d_b50=np.empty(simulations)
part_d_b60=np.empty(simulations)
part_d_b75=np.empty(simulations)

of_d=np.empty(simulations)
of_d_b5=np.empty(simulations)
of_d_b10=np.empty(simulations)
of_d_b15=np.empty(simulations)
of_d_b16=np.empty(simulations)
of_d_b17=np.empty(simulations)
of_d_b18=np.empty(simulations)
of_d_b19=np.empty(simulations)
of_d_b20=np.empty(simulations)
of_d_b30=np.empty(simulations)
of_d_b40=np.empty(simulations)
of_d_b50=np.empty(simulations)
of_d_b60=np.empty(simulations)
of_d_b75=np.empty(simulations)

part_nd=np.empty(simulations)
part_nd_b5=np.empty(simulations)
part_nd_b10=np.empty(simulations)
part_nd_b15=np.empty(simulations)
part_nd_b16=np.empty(simulations)
part_nd_b17=np.empty(simulations)
part_nd_b18=np.empty(simulations)
part_nd_b19=np.empty(simulations)
part_nd_b20=np.empty(simulations)
part_nd_b30=np.empty(simulations)
part_nd_b40=np.empty(simulations)
part_nd_b50=np.empty(simulations)
part_nd_b60=np.empty(simulations)
part_nd_b75=np.empty(simulations)


diff=np.empty(simulations)
diff_b5=np.empty(simulations)
diff_b10=np.empty(simulations)
diff_b15=np.empty(simulations)
diff_b16=np.empty(simulations)
diff_b17=np.empty(simulations)
diff_b18=np.empty(simulations)
diff_b19=np.empty(simulations)
diff_b20=np.empty(simulations)
diff_b30=np.empty(simulations)
diff_b40=np.empty(simulations)
diff_b50=np.empty(simulations)
diff_b60=np.empty(simulations)
diff_b75=np.empty(simulations)


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
               of_nd[sim]=calc_measof(node_list)
        
               g_d = Graph(pop+1) 
      
               for i in range(0,pop):
                   for j in range(i+1,pop+1):
                       g_d.addEdge(i,j,wt_D[i][j])
                
               node_list=g_d.shortestPath(s) 
               part_d[sim]=len(node_list)
               of_d[sim]=calc_measof(node_list)
               
               
               


               of_nd_b5[sim],part_nd_b5[sim]=Bellmanof(wt_ND,5)
               of_nd_b10[sim],part_nd_b10[sim]=Bellmanof(wt_ND,10)
               of_nd_b15[sim],part_nd_b15[sim]=Bellmanof(wt_ND,15)
               of_nd_b16[sim],part_nd_b16[sim]=Bellmanof(wt_ND,16)
               of_nd_b17[sim],part_nd_b17[sim]=Bellmanof(wt_ND,17)
               of_nd_b18[sim],part_nd_b18[sim]=Bellmanof(wt_ND,18)
               of_nd_b19[sim],part_nd_b19[sim]=Bellmanof(wt_ND,19)
               of_nd_b20[sim],part_nd_b20[sim]=Bellmanof(wt_ND,20)     
               of_nd_b30[sim],part_nd_b30[sim]=Bellmanof(wt_ND,30)   
               of_nd_b40[sim],part_nd_b40[sim]=Bellmanof(wt_ND,40)    
               of_nd_b50[sim],part_nd_b50[sim]=Bellmanof(wt_ND,50)      
               of_nd_b60[sim],part_nd_b60[sim]=Bellmanof(wt_ND,60) 
               of_nd_b75[sim],part_nd_b75[sim]=Bellmanof(wt_ND,75)    
               
               of_d_b5[sim],part_d_b5[sim]=Bellmanof(wt_D,5)
               of_d_b10[sim],part_d_b10[sim]=Bellmanof(wt_D,10)
               of_d_b15[sim],part_d_b15[sim]=Bellmanof(wt_D,15)
               of_d_b16[sim],part_d_b16[sim]=Bellmanof(wt_D,16)
               of_d_b17[sim],part_d_b17[sim]=Bellmanof(wt_D,17)
               of_d_b18[sim],part_d_b18[sim]=Bellmanof(wt_D,18)
               of_d_b19[sim],part_d_b19[sim]=Bellmanof(wt_D,19)
               of_d_b20[sim],part_d_b20[sim]=Bellmanof(wt_D,20)     
               of_d_b30[sim],part_d_b30[sim]=Bellmanof(wt_D,30)   
               of_d_b40[sim],part_d_b40[sim]=Bellmanof(wt_D,40)    
               of_d_b50[sim],part_d_b50[sim]=Bellmanof(wt_D,50)      
               of_d_b60[sim],part_d_b60[sim]=Bellmanof(wt_D,60) 
               of_d_b75[sim],part_d_b75[sim]=Bellmanof(wt_D,75)  
               
               diff=(of_d[sim]-of_nd[sim])*100/of_nd[sim]
               diff_b5=(of_d_b5[sim]-of_nd_b5[sim])*100/of_nd_b5[sim]
               diff_b10=(of_d_b10[sim]-of_nd_b10[sim])*100/of_nd_b10[sim] 
               diff_b15=(of_d_b15[sim]-of_nd_b15[sim])*100/of_nd_b15[sim] 
               diff_b16=(of_d_b16[sim]-of_nd_b16[sim])*100/of_nd_b16[sim]
               diff_b17=(of_d_b17[sim]-of_nd_b17[sim])*100/of_nd_b17[sim]
               diff_b18=(of_d_b18[sim]-of_nd_b18[sim])*100/of_nd_b18[sim]
               diff_b19=(of_d_b19[sim]-of_nd_b19[sim])*100/of_nd_b19[sim]
               diff_b20=(of_d_b20[sim]-of_nd_b20[sim])*100/of_nd_b20[sim]   
               diff_b30=(of_d_b30[sim]-of_nd_b30[sim])*100/of_nd_b30[sim]   
               diff_b40=(of_d_b40[sim]-of_nd_b40[sim])*100/of_nd_b40[sim]  
               diff_b50=(of_d_b50[sim]-of_nd_b50[sim])*100/of_nd_b50[sim]     
               diff_b60=(of_d_b60[sim]-of_nd_b60[sim])*100/of_nd_b60[sim]     
               diff_b75=(of_d_b75[sim]-of_nd_b75[sim])*100/of_nd_b75[sim]               
 
              
results = results.append({           
                                 'OF_ND':round(np.mean(of_nd),4),
                                 'OF_ND_B5':round(np.mean(of_nd_b5),4),
                                 'OF_ND_B10':round(np.mean(of_nd_b10),4),    
                                 'OF_ND_B15':round(np.mean(of_nd_b15),4),    
                                 'OF_ND_B16':round(np.mean(of_nd_b16),4),
                                 'OF_ND_B17':round(np.mean(of_nd_b17),4),
                                 'OF_ND_B18':round(np.mean(of_nd_b18),4),
                                 'OF_ND_B19':round(np.mean(of_nd_b19),4),
                                 'OF_ND_B20':round(np.mean(of_nd_b20),4), 
                                 'OF_ND_B30':round(np.mean(of_nd_b30),4),  
                                 'OF_ND_B40':round(np.mean(of_nd_b40),4), 
                                 'OF_ND_B50':round(np.mean(of_nd_b50),4),   
                                 'OF_ND_B60':round(np.mean(of_nd_b60),4),  
                                 'OF_ND_B75':round(np.mean(of_nd_b75),4),
                                 
                                 'PART_ND':round(np.mean(part_nd),4),
                                 'PART_ND_B5':round(np.mean(part_nd_b5),4),
                                 'PART_ND_B10':round(np.mean(part_nd_b10),4),    
                                 'PART_ND_B15':round(np.mean(part_nd_b15),4),  
                                 'PART_ND_B16':round(np.mean(part_nd_b16),4),
                                 'PART_ND_B17':round(np.mean(part_nd_b17),4),
                                 'PART_ND_B18':round(np.mean(part_nd_b18),4),
                                 'PART_ND_B19':round(np.mean(part_nd_b19),4),
                                 'PART_ND_B20':round(np.mean(part_nd_b20),4), 
                                 'PART_ND_B30':round(np.mean(part_nd_b30),4),  
                                 'PART_ND_B40':round(np.mean(part_nd_b40),4), 
                                 'PART_ND_B50':round(np.mean(part_nd_b50),4),   
                                 'PART_ND_B60':round(np.mean(part_nd_b60),4),  
                                 'PART_ND_B75':round(np.mean(part_nd_b75),4),                                 
                                 
                                 'OF_D':round(np.mean(of_d),4),
                                 'OF_D_B5':round(np.mean(of_d_b5),4),
                                 'OF_D_B10':round(np.mean(of_d_b10),4),    
                                 'OF_D_B15':round(np.mean(of_d_b15),4),    
                                 'OF_D_B16':round(np.mean(of_d_b16),4), 
                                 'OF_D_B17':round(np.mean(of_d_b17),4), 
                                 'OF_D_B18':round(np.mean(of_d_b18),4), 
                                 'OF_D_B19':round(np.mean(of_d_b19),4), 
                                 'OF_D_B20':round(np.mean(of_d_b20),4), 
                                 'OF_D_B30':round(np.mean(of_d_b30),4),  
                                 'OF_D_B40':round(np.mean(of_d_b40),4), 
                                 'OF_D_B50':round(np.mean(of_d_b50),4),   
                                 'OF_D_B60':round(np.mean(of_d_b60),4),  
                                 'OF_D_B75':round(np.mean(of_d_b75),4), 
                                 
                                 'PART_D':round(np.mean(part_d),4),
                                 'PART_D_B5':round(np.mean(part_d_b5),4),
                                 'PART_D_B10':round(np.mean(part_d_b10),4),    
                                 'PART_D_B15':round(np.mean(part_d_b15),4), 
                                 'PART_D_B16':round(np.mean(part_d_b16),4),
                                 'PART_D_B17':round(np.mean(part_d_b17),4),
                                 'PART_D_B18':round(np.mean(part_d_b18),4),
                                 'PART_D_B19':round(np.mean(part_d_b19),4),
                                 'PART_D_B20':round(np.mean(part_d_b20),4), 
                                 'PART_D_B30':round(np.mean(part_d_b30),4),  
                                 'PART_D_B40':round(np.mean(part_d_b40),4), 
                                 'PART_D_B50':round(np.mean(part_d_b50),4),   
                                 'PART_D_B60':round(np.mean(part_d_b60),4),  
                                 'PART_D_B75':round(np.mean(part_d_b75),4), 
                                 
                                 'DIFF':round(np.mean(diff),4),
                                 'DIFF_B5':round(np.mean(diff_b5),4),
                                 'DIFF_B10':round(np.mean(diff_b10),4),    
                                 'DIFF_B15':round(np.mean(diff_b15),4),  
                                 'DIFF_B16':round(np.mean(diff_b16),4),
                                 'DIFF_B17':round(np.mean(diff_b17),4),
                                 'DIFF_B18':round(np.mean(diff_b18),4),
                                 'DIFF_B19':round(np.mean(diff_b19),4),
                                 'DIFF_B20':round(np.mean(diff_b20),4), 
                                 'DIFF_B30':round(np.mean(diff_b30),4),  
                                 'DIFF_B40':round(np.mean(diff_b40),4), 
                                 'DIFF_B50':round(np.mean(diff_b50),4),   
                                 'DIFF_B60':round(np.mean(diff_b60),4),  
                                 'DIFF_B75':round(np.mean(diff_b75),4)
    
                                },ignore_index=True)


# #%% 7. SIMULATIONS
# # sys.stdout = open("FINALRESULTS.txt", "w")
# results = pd.DataFrame(columns= ['w1',
#                                  'w2',
#                                  'OF_ND',
#                                  'OF_ND_B2',
#                                  'OF_ND_B3',
#                                  'OF_ND_B4',
#                                  'OF_ND_B5',
#                                  'OF_ND_B10',
#                                  'OF_D',
#                                  'OF_D_B2',
#                                  'OF_D_B3',
#                                  'OF_D_B4',
#                                  'OF_D_B5',
#                                  'OF_D_B10',
#                                  'ET_ND',
#                                  'ET_ND_B2',
#                                  'ET_ND_B3',
#                                  'ET_ND_B4',
#                                  'ET_ND_B5',
#                                  'ET_ND_B10',
#                                  'ET_D',
#                                  'ET_D_B2',
#                                  'ET_D_B3',
#                                  'ET_D_B4',
#                                  'ET_D_B5',
#                                  'ET_D_B10',
#                                  'EFP_ND',
#                                  'EFP_ND_B2',
#                                  'EFP_ND_B3',
#                                  'EFP_ND_B4',
#                                  'EFP_ND_B5',
#                                  'EFP_ND_B10',
#                                  'EFP_D',
#                                  'EFP_D_B2',
#                                  'EFP_D_B3',
#                                  'EFP_D_B4',
#                                  'EFP_D_B5',
#                                  'EFP_D_B10',
#                                  'EFN_ND',
#                                  'EFN_ND_B2',
#                                  'EFN_ND_B3',
#                                  'EFN_ND_B4',
#                                  'EFN_ND_B5',
#                                  'EFN_ND_B10',
#                                  'EFN_D',
#                                  'EFN_D_B2',
#                                  'EFN_D_B3',
#                                  'EFN_D_B4',
#                                  'EFN_D_B5',
#                                  'EFN_D_B10',
#                                  'PART_ND',
#                                  'PART_ND_B2',
#                                  'PART_ND_B3',
#                                  'PART_ND_B4',
#                                  'PART_ND_B5',
#                                  'PART_ND_B10',
#                                  'PART_D',
#                                  'PART_D_B2',
#                                  'PART_D_B3',
#                                  'PART_D_B4',
#                                  'PART_D_B5',
#                                  'PART_D_B10']) 
           
# for w1 in tqdm(range(0, 12, 2),desc='lambda loop',position=0,leave=True):
#    for w2 in range(0,12,2):
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
           
#            of_nd_b2=np.empty(simulations)
#            et_nd_b2=np.empty(simulations)
#            efp_nd_b2=np.empty(simulations)
#            efn_nd_b2=np.empty(simulations)
#            part_nd_b2=np.empty(simulations)
           
#            of_nd_b3=np.empty(simulations)
#            et_nd_b3=np.empty(simulations)
#            efp_nd_b3=np.empty(simulations)
#            efn_nd_b3=np.empty(simulations)
#            part_nd_b3=np.empty(simulations)
           
#            of_nd_b4=np.empty(simulations)
#            et_nd_b4=np.empty(simulations)
#            efp_nd_b4=np.empty(simulations)
#            efn_nd_b4=np.empty(simulations)
#            part_nd_b4=np.empty(simulations)
           
#            of_nd_b5=np.empty(simulations)
#            et_nd_b5=np.empty(simulations)
#            efp_nd_b5=np.empty(simulations)
#            efn_nd_b5=np.empty(simulations)
#            part_nd_b5=np.empty(simulations)
           
#            of_nd_b10=np.empty(simulations)
#            et_nd_b10=np.empty(simulations)
#            efp_nd_b10=np.empty(simulations)
#            efn_nd_b10=np.empty(simulations)
#            part_nd_b10=np.empty(simulations)
           
#            of_d=np.empty(simulations)
#            et_d=np.empty(simulations)
#            efp_d=np.empty(simulations)
#            efn_d=np.empty(simulations)
#            part_d=np.empty(simulations)
           
#            of_d_b2=np.empty(simulations)
#            et_d_b2=np.empty(simulations)
#            efp_d_b2=np.empty(simulations)
#            efn_d_b2=np.empty(simulations)
#            part_d_b2=np.empty(simulations)
           
#            of_d_b3=np.empty(simulations)
#            et_d_b3=np.empty(simulations)
#            efp_d_b3=np.empty(simulations)
#            efn_d_b3=np.empty(simulations)
#            part_d_b3=np.empty(simulations)
           
#            of_d_b4=np.empty(simulations)
#            et_d_b4=np.empty(simulations)
#            efp_d_b4=np.empty(simulations)
#            efn_d_b4=np.empty(simulations)
#            part_d_b4=np.empty(simulations)
           
#            of_d_b5=np.empty(simulations)
#            et_d_b5=np.empty(simulations)
#            efp_d_b5=np.empty(simulations)
#            efn_d_b5=np.empty(simulations)
#            part_d_b5=np.empty(simulations)
           
#            of_d_b10=np.empty(simulations)
#            et_d_b10=np.empty(simulations)
#            efp_d_b10=np.empty(simulations)
#            efn_d_b10=np.empty(simulations)
#            part_d_b10=np.empty(simulations)
           
#            for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
#                get_risk_vector()
               
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
               
#                of_nd_b2[sim],et_nd_b2[sim],efp_nd_b2[sim],efn_nd_b2[sim],part_nd_b2[sim]=Bellman(wt_ND,2)
#                of_nd_b3[sim],et_nd_b3[sim],efp_nd_b3[sim],efn_nd_b3[sim],part_nd_b3[sim]=Bellman(wt_ND,3)
#                of_nd_b4[sim],et_nd_b4[sim],efp_nd_b4[sim],efn_nd_b4[sim],part_nd_b4[sim]=Bellman(wt_ND,4)
#                of_nd_b5[sim],et_nd_b5[sim],efp_nd_b5[sim],efn_nd_b5[sim],part_nd_b5[sim]=Bellman(wt_ND,5)
#                of_nd_b10[sim],et_nd_b10[sim],efp_nd_b10[sim],efn_nd_b10[sim],part_nd_b10[sim]=Bellman(wt_ND,10)

#                of_d_b2[sim],et_d_b2[sim],efp_d_b2[sim],efn_d_b2[sim],part_d_b2[sim]=Bellman(wt_D,2)
#                of_d_b3[sim],et_d_b3[sim],efp_d_b3[sim],efn_d_b3[sim],part_d_b3[sim]=Bellman(wt_D,3)
#                of_d_b4[sim],et_d_b4[sim],efp_d_b4[sim],efn_d_b4[sim],part_d_b4[sim]=Bellman(wt_D,4)
#                of_d_b5[sim],et_d_b5[sim],efp_d_b5[sim],efn_d_b5[sim],part_d_b5[sim]=Bellman(wt_D,5)
#                of_d_b10[sim],et_d_b10[sim],efp_d_b10[sim],efn_d_b10[sim],part_d_b10[sim]=Bellman(wt_D,10)

#            results = results.append({'w1':lambda_1,
#                                      'w2':lambda_2,
#                                      'OF_ND':round(np.mean(of_nd),4),
#                                      'OF_ND_B2':round(np.mean(of_nd_b2),4),
#                                      'OF_ND_B3':round(np.mean(of_nd_b3),4),
#                                      'OF_ND_B4':round(np.mean(of_nd_b4),4),
#                                      'OF_ND_B5':round(np.mean(of_nd_b5),4),
#                                      'OF_ND_B10':round(np.mean(of_nd_b10),4),
                                     
#                                      'OF_D':round(np.mean(of_d),4),
#                                      'OF_D_B2':round(np.mean(of_d_b2),4),
#                                      'OF_D_B3':round(np.mean(of_d_b3),4),
#                                      'OF_D_B4':round(np.mean(of_d_b4),4),
#                                      'OF_D_B5':round(np.mean(of_d_b5),4),
#                                      'OF_D_B10':round(np.mean(of_d_b10),4),
                                     
#                                      'ET_ND':round(np.mean(et_nd),4),
#                                      'ET_ND_B2':round(np.mean(et_nd_b2),4),
#                                      'ET_ND_B3':round(np.mean(et_nd_b3),4),
#                                      'ET_ND_B4':round(np.mean(et_nd_b4),4),
#                                      'ET_ND_B5':round(np.mean(et_nd_b5),4),
#                                      'ET_ND_B10':round(np.mean(et_nd_b10),4),
                                     
#                                      'ET_D':round(np.mean(et_d),4),
#                                      'ET_D_B2':round(np.mean(et_d_b2),4),
#                                      'ET_D_B3':round(np.mean(et_d_b3),4),
#                                      'ET_D_B4':round(np.mean(et_d_b4),4),
#                                      'ET_D_B5':round(np.mean(et_d_b5),4),
#                                      'ET_D_B10':round(np.mean(et_d_b10),4),

#                                      'EFP_ND':round(np.mean(efp_nd),4),
#                                      'EFP_ND_B2':round(np.mean(efp_nd_b2),4),
#                                      'EFP_ND_B3':round(np.mean(efp_nd_b3),4),
#                                      'EFP_ND_B4':round(np.mean(efp_nd_b4),4),
#                                      'EFP_ND_B5':round(np.mean(efp_nd_b5),4),
#                                      'EFP_ND_B10':round(np.mean(efp_nd_b10),4),
                                     
#                                      'EFP_D':round(np.mean(efp_d),4),
#                                      'EFP_D_B2':round(np.mean(efp_d_b2),4),
#                                      'EFP_D_B3':round(np.mean(efp_d_b3),4),
#                                      'EFP_D_B4':round(np.mean(efp_d_b4),4),
#                                      'EFP_D_B5':round(np.mean(efp_d_b5),4),
#                                      'EFP_D_B10':round(np.mean(efp_d_b10),4),
                                     
#                                      'EFN_ND':round(np.mean(efn_nd),4),
#                                      'EFN_ND_B2':round(np.mean(efn_nd_b2),4),
#                                      'EFN_ND_B3':round(np.mean(efn_nd_b3),4),
#                                      'EFN_ND_B4':round(np.mean(efn_nd_b4),4),
#                                      'EFN_ND_B5':round(np.mean(efn_nd_b5),4),
#                                      'EFN_ND_B10':round(np.mean(efn_nd_b10),4),
                                     
#                                      'EFN_D':round(np.mean(efn_d),4),
#                                      'EFN_D_B2':round(np.mean(efn_d_b2),4),
#                                      'EFN_D_B3':round(np.mean(efn_d_b3),4),
#                                      'EFN_D_B4':round(np.mean(efn_d_b4),4),
#                                      'EFN_D_B5':round(np.mean(efn_d_b5),4),
#                                      'EFN_D_B10':round(np.mean(efn_d_b10),4),
                                     
#                                      'PART_ND':round(np.mean(part_nd),4),
#                                      'PART_ND_B2':round(np.mean(part_nd_b2),4),
#                                      'PART_ND_B3':round(np.mean(part_nd_b3),4),
#                                      'PART_ND_B4':round(np.mean(part_nd_b4),4),
#                                      'PART_ND_B5':round(np.mean(part_nd_b5),4),
#                                      'PART_ND_B10':round(np.mean(part_nd_b10),4),
                                     
#                                      'PART_D':round(np.mean(part_d),4),
#                                      'PART_D_B2':round(np.mean(part_d_b2),4),
#                                      'PART_D_B3':round(np.mean(part_d_b3),4),
#                                      'PART_D_B4':round(np.mean(part_d_b4),4),
#                                      'PART_D_B5':round(np.mean(part_d_b5),4),
#                                      'PART_D_B10':round(np.mean(part_d_b10),4)},ignore_index=True)
           
#            # results = results.append({'w1':lambda_1,'w2':lambda_2,'OF_ND':round(np.mean(of_nd),4),'OF_D':round(np.mean(of_d),4),'Change_OF':(round(np.mean(of_d),4)-round(np.mean(of_nd),4))*100/round(np.mean(of_nd),4),'OF_B5':round(np.mean(of_b5),4),'OF_B10':round(np.mean(of_b10),4),'ET_ND':round(np.mean(et_nd),4),'ET_D':round(np.mean(et_d),4),'Change_ET':(round(np.mean(et_d),4)-round(np.mean(et_nd),4))*100/round(np.mean(et_nd),4),'ET_B5':round(np.mean(et_b5),4),'ET_B10':round(np.mean(et_b10),4),'EFP_ND':round(np.mean(efp_nd),4),'EFP_D':round(np.mean(efp_d),4),'Change_EFP':(round(np.mean(efp_d),4)-round(np.mean(efp_nd),4))*100/round(np.mean(efp_nd),4),'EFP_B5':round(np.mean(efp_b5),4),'EFP_B10':round(np.mean(efp_b10),4),'EFN_ND':round(np.mean(efn_nd),4),'EFN_D':round(np.mean(efn_d),4),'Change_EFN':(round(np.mean(efn_d),4)-round(np.mean(efn_nd),4))*100/round(np.mean(efn_nd),4),'EFN_B5':round(np.mean(efn_b5),4),'EFN_B10':round(np.mean(efn_b10),4),},ignore_index=True)
#            # results = results.append({'w1':lambda_1,'w2':lambda_2,'OF_ND':round(np.mean(of_nd),4),'OF_D':round(np.mean(of_d),4),'Change_OF':(round(np.mean(of_d),4)-round(np.mean(of_nd),4))*100/round(np.mean(of_nd),4),'ET_ND':round(np.mean(et_nd),4),'ET_D':round(np.mean(et_d),4),'Change_ET':(round(np.mean(et_d),4)-round(np.mean(et_nd),4))*100/round(np.mean(et_nd),4),'EFP_ND':round(np.mean(efp_nd),4),'EFP_D':round(np.mean(efp_d),4),'Change_EFP':(round(np.mean(efp_d),4)-round(np.mean(efp_nd),4))*100/round(np.mean(efp_nd),4),'EFN_ND':round(np.mean(efn_nd),4),'EFN_D':round(np.mean(efn_d),4),'Change_EFN':(round(np.mean(efn_d),4)-round(np.mean(efn_nd),4))*100/round(np.mean(efn_nd),4)},ignore_index=True)
#            # #PRINT RESULTS
#            # print("OF_ND:\n")
#            # interval_of_nd = stats.t.ppf(1.0 - (CI / 2.0),of_nd.size-1) * (np.std(of_nd)/ np.sqrt(of_nd.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(of_nd),4))+"+-"+str(round(interval_of_nd,4)))
           
#            # print("\nET_ND:\n")
#            # interval_et_nd = stats.t.ppf(1.0 - (CI / 2.0),et_nd.size-1) * (np.std(et_nd)/ np.sqrt(et_nd.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(et_nd),4))+"+-"+str(round(interval_et_nd,4)))
           
#            # print("\nEFP_ND:\n")
#            # interval_efp_nd = stats.t.ppf(1.0 - (CI / 2.0),efp_nd.size-1) * (np.std(efp_nd)/ np.sqrt(efp_nd.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(efp_nd),4))+"+-"+str(round(interval_efp_nd,4)))
          
#            # print("\nEFN_ND:\n")
#            # interval_efn_nd = stats.t.ppf(1.0 - (CI / 2.0),efn_nd.size-1) * (np.std(efn_nd)/ np.sqrt(efn_nd.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(efn_nd),4))+"+-"+str(round(interval_efn_nd,4)))
           
#            # print('\n')
#            # print('--------------------------------')
#            # print('\n')
#            # print("OF_D:\n")
#            # interval_of_d = stats.t.ppf(1.0 - (CI / 2.0),of_d.size-1) * (np.std(of_d)/ np.sqrt(of_d.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(of_d),4))+"+-"+str(round(interval_of_d,4)))
           
#            # print("\nET_D:\n")
#            # interval_et_d = stats.t.ppf(1.0 - (CI / 2.0),et_d.size-1) * (np.std(et_d)/ np.sqrt(et_d.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(et_d),4))+"+-"+str(round(interval_et_d,4)))
           
#            # print("\nEFP_D:\n")
#            # interval_efp_d = stats.t.ppf(1.0 - (CI / 2.0),efp_d.size-1) * (np.std(efp_d)/ np.sqrt(efp_d.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(efp_d),4))+"+-"+str(round(interval_efp_d,4)))
           
#            # print("\nEFN_D:\n")
#            # interval_efn_d = stats.t.ppf(1.0 - (CI / 2.0),efn_d.size-1) * (np.std(efn_d)/ np.sqrt(efn_d.size))
#            # # ci = (a.mean() - interval, a.mean() + interval)
#            # print(str(round(np.mean(efn_d),4))+"+-"+str(round(interval_efn_d,4)))
           
           
# # sys.stdout.close()
#%% 8. OUTPUT TO CSV

results.to_csv('Bvar_100_400sims_16to192.csv'
#%%             
budget=pd.read_csv("C:/Users/Sohom/Desktop/Group Testing/Budget.csv")  

x=budget["Budget"]
y=budget["Difference"]
plt.xlabel("Maximum allowable groups")
plt.ylabel("Percent Improvement")
# plot the data itself
# pylab.plot(x,y,'o')
# plt.axhline(max(y), c= 'red', ls='--')
plt.xlim(-5,105)
plt.ylim(0, 30)
plt.title("Improvement of proposed scheme over existing results") 
# calc the trendline
z = numpy.polyfit(x, y, 15)
p = numpy.poly1d(z)
pylab.plot(x,p(x),"b--")
# the line equation:
# print "y=%.6fx+(%.6f)"%(z[0],z[1])