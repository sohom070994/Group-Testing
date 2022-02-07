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

#%% 1.Import Libraries
import matplotlib.pyplot as plt
import random,math 
from tqdm import tqdm
import numpy as np 
import scipy.stats as stats
#%% 2. Establish constants


#Total population of subjects
pop = 100
#Risk vector initiate
# risk_vec = [0] * pop


#Defining weights for Objective Function
lambda_1 = 0.96 #EFN
lambda_2 = 0.02 #EFP

#testing constants
spec = 0.99
alpha= 0.0514
gamma= 0.2416

#confidence interval
CI = 0.05

#simulations
simulations=1
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
  if 0.00 <= randnum <= 16.34:
    return 0.09*3
  elif randnum <= 33.45:
    return 0.48*3
  elif randnum <= 48.89:
    return 1.07*3
  elif randnum <= 64.91:
    return 0.81*3
  elif randnum <= 80.78:
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


def Exp_FalseNeg(Node_1,Node_2): 
  # if Node_2==Node_1+1: 
  #   return (1-sens)*risk_vec[Node_1]
  # else:
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
    
def Exp_FalsePos(Node_1,Node_2):
  # if Node_2==Node_1+1: 
  #   return (1-spec)*(1-risk_vec[Node_1])
  # else:
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

def Exp_Test_Num(Node_1,Node_2): 
  # if Node_2==Node_1+1: 
  #   return 1
  # else:
    wt = 1
    for ele in range(Node_1,Node_2):
      wt = wt*(1-alpha_pr(Node_1,Node_2)*risk_vec[ele])
    return 1+(Node_2-Node_1)*(1-spec*(wt))
  

'''
TOTAL WEIGHT (OBJ FUNC VALUE)
FOR EVERY EDGE
'''
def weight(Node_1,Node_2):
  weight = (lambda_1*Exp_FalseNeg(Node_1,Node_2))+(lambda_2*Exp_FalsePos(Node_1,Node_2))+((1-lambda_1-lambda_2)*Exp_Test_Num(Node_1,Node_2))
  return weight
#%% 6. Shortest Path
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

    # Print the distance from source to Node_100
    # for i in range(self.V):
    # print ("Dist - %f, Node - %d" %(dist[100],100)) #if dist[i] != float ("Inf") else "Inf" ,
    print(node_list)
    return dist[100]
#%% 7. Run Simulations

a=np.empty(simulations)
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
  a[sim] = g.shortestPath(s) 


len(a)
#%% Results

#95% CI

interval = stats.t.ppf(1.0 - (CI / 2.0),a.size-1) * (np.std(a)/ np.sqrt(a.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(a),4))+"+-"+str(round(interval,4)))
#%%
of=0
of+=weight(node_list[0],pop)
for ele in range(0,len(node_list)-1):
    of+=weight(node_list[ele+1],node_list[ele])