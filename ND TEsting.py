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
#%% 1. Import Libraries
import matplotlib.pyplot as plt
import random,math 
from tqdm import tqdm
import numpy as np 
import scipy.stats as stats

#%% 2. Establish constants
#population
pop = 100
# risk_vec = [0] * pop


#Defining weights for Objective Function
lambda_1 = 0.96
lambda_2 = 0.02
# sens = 0.95
sens = 0.949114
# spec = 0.95
spec = 0.99
#confidence interval
CI = 0.05

#simulations
simulations=1

# partition = [9999] * (pop+1)
#%% 3. Functions
'''
LIST OF FUNCTIONS 
FOR GENERATING RISK VECTOR
'''


def coin_flip():
  '''
  get random float between 0 to 100 
  round to 2 digits
  check range and return proper risk
  '''
  randnum = round(random.uniform(0,100),2)
  if 0.00 <= randnum <= 18.75:
    return 0.00
  elif randnum <= 32.02:
    return 0.09*3
  elif randnum <= 45.93:
    return 0.48*3
  elif randnum <= 58.47:
    return 1.07*3
  elif randnum <= 71.49:
    return 0.81*3
  elif randnum <= 84.38:
    return 0.42*3
  else:
    return 0.18*3

def get_risk_vector():
  '''
  fill risk_vec
  using coin flips
  sort risk_vec
  '''
  for i in range(0,pop):
    risk_vec[i]=round(coin_flip()/100,4)

  risk_vec.sort()
#%% 4. Debugging Block
# get_risk_vector()
# print(risk_vec)
# print(len(risk_vec))
#%% 5. Performance Measures
'''
LIST OF 
PERFORMANCE MEASURES
'''

def Exp_FalseNeg(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return (1-sens)*risk_vec[Node_1]
  else:
    wt = 0
    for ele in range(Node_1,Node_2):
      wt = wt + risk_vec[ele]
    return (1-sens**2)*wt  
    
def Exp_FalsePos(Node_1,Node_2):
  if Node_2==Node_1+1: 
    return (1-spec)*(1-risk_vec[Node_1])
  else:
    sum = 0 
    for ele in range(Node_1,Node_2):
      sum = sum + (1-risk_vec[ele])
    wt = 1
    for ele in range(Node_1,Node_2):
      wt = wt*(1-risk_vec[ele]) 
    return ((1-spec)*sens*sum)-((Node_2-Node_1)*(1-spec)*(sens+spec-1)*(wt))

def Exp_Test_Num(Node_1,Node_2): 
  if Node_2==Node_1+1: 
    return 1
  else:
    wt = 1
    for ele in range(Node_1,Node_2):
      wt = wt*(1-risk_vec[ele])
    return 1+(Node_2-Node_1)*(sens-(sens+spec-1)*(wt))
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg(Node_1,Node_2))+(lambda_2*Exp_FalsePos(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num(Node_1,Node_2))
  return weight
#%% 6. Topological Sort
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
    
    et=0
    et+=Exp_Test_Num(node_list[0],pop)
    efp=0
    efp+=Exp_FalsePos(node_list[0],pop)
    efn=0
    efn+=Exp_FalseNeg(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        et+=Exp_Test_Num(node_list[ele+1],node_list[ele])
        efp+=Exp_FalsePos(node_list[ele+1],node_list[ele])
        efn+=Exp_FalseNeg(node_list[ele+1],node_list[ele])
        
    print(node_list)
    return dist[pop],et,efp,efn
# 
#%% 7. Simulations
a=np.empty(simulations)
t=np.empty(simulations)
fp=np.empty(simulations)
fn=np.empty(simulations)
'''
CREATING THE NETWORK
USING PROPER WEIGHTS
'''
for sim in tqdm(range(0,simulations)):
  # get_risk_vector()
  g = Graph(pop+1) 

  for i in range(0,pop):
    for j in range(i+1,pop+1):
      g.addEdge(i,j,weight(i,j))


  # g.graph

  # # source = 0
  s = 0

  # # print ("Following are shortest distances from source %d " % s) 
  a[sim],t[sim],fp[sim],fn[sim] = g.shortestPath(s) 


len(a)
#%% 8. Results
print("OF:\n")
interval_a = stats.t.ppf(1.0 - (CI / 2.0),a.size-1) * (np.std(a)/ np.sqrt(a.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(a),4))+"+-"+str(round(interval_a,4)))
print("\nET:\n")
interval_t = stats.t.ppf(1.0 - (CI / 2.0),t.size-1) * (np.std(t)/ np.sqrt(t.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(t),4))+"+-"+str(round(interval_t,4)))
print("\nEFP:\n")
interval_fp = stats.t.ppf(1.0 - (CI / 2.0),fp.size-1) * (np.std(fp)/ np.sqrt(fp.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(fp),4))+"+-"+str(round(interval_fp,4)))
print("\nEFN:\n")
interval_fn = stats.t.ppf(1.0 - (CI / 2.0),fn.size-1) * (np.std(fn)/ np.sqrt(fn.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(fn),4))+"+-"+str(round(interval_fn,4)))
#%%
sum+=Exp_FalsePos(95,100)
print(sum)