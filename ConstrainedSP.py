# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%% 1. IMPORT LIBRARIES
import random
from tqdm import tqdm
import numpy as np 
import scipy.stats as stats
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
#%% 3.FUNCTIONS

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
                
    node_list=[]
    a=pathNode[B][pop]
    while (a!=0):
        node_list.append(a)
        a=pathNode[B][a]
    node_list.append(0) 


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
        
    return d[B][pop],et,efp,efn,len(node_list)

#%% 4. Shortest Path
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
    et+=Exp_Test_Num_D(node_list[0],pop)
    efp=0
    efp+=Exp_FalsePos_D(node_list[0],pop)
    efn=0
    efn+=Exp_FalseNeg_D(node_list[0],pop)
    for ele in range(0,len(node_list)-1):
        et+=Exp_Test_Num_D(node_list[ele+1],node_list[ele])
        efp+=Exp_FalsePos_D(node_list[ele+1],node_list[ele])
        efn+=Exp_FalseNeg_D(node_list[ele+1],node_list[ele])
        
    # print(part)
    return dist[pop],et,efp,efn
# 
#%% 5. ESTABLISH CONSTANT PARAMETERS
#population
pop = 100
risk_vec = [0] * pop

#Defining weights for Objective Function
lambda_1 = 1.0 #weight on EFN
lambda_2 = 0.0 #weight on EFP
# get_risk_vector()
#Parameters
spec = 0.99
alpha= 0.0514
gamma= 0.2416
const_sens = sens(1,1)

#Confidence Interval
CI = 0.05

# Simulations
simulations=200
#%% 6. Simulations
a=np.empty(simulations)
t=np.empty(simulations)
fp=np.empty(simulations)
fn=np.empty(simulations)
'''
CREATING THE NETWORK
USING PROPER WEIGHTS
'''



for sim in tqdm(range(0,simulations)):
  get_risk_vector()
  # risk_vec=[risk_vec[i]*.97/100 for i in range(len(risk_vec))]

 
  wt =[[None for i in range(pop+1)]for i in range(pop+1)]
  for i in range(0,pop):
      for j in range(i+1,pop+1):
          wt[i][j]=weight_D(i,j)
          
  g = Graph(pop+1)        
  for i in range(0,pop):
    for j in range(i+1,pop+1):
      g.addEdge(i,j,wt[i][j])


  # g.graph

  # # source = 0
  s = 0

  # # print ("Following are shortest distances from source %d " % s) 
  a[sim],t[sim],fp[sim],fn[sim] = g.shortestPath(s) 


# len(a)
#%% 7. Results
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

#%% 8. Bellman-Ford

of_b=np.empty(simulations)
t_b=np.empty(simulations)
fp_b=np.empty(simulations)
fn_b=np.empty(simulations)
part_b=np.empty(simulations)
for sim in tqdm(range(0,simulations)):
  get_risk_vector()
  
  wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
  for i in range(0,pop):
      for j in range(i+1,pop+1):
          wt_D[i][j]=weight_D(i,j)
          
  # B= max allowable edges         
  B=pop+1
  of_b[sim],t_b[sim],fp_b[sim],fn_b[sim],part_b[sim]=Bellman(B)
#%% 9. Bellman Ford Results
print("OF Bellman:\n")
interval_of_b = stats.t.ppf(1.0 - (CI / 2.0),of_b.size-1) * (np.std(of_b)/ np.sqrt(of_b.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(of_b),4))+"+-"+str(round(interval_of_b,4)))
print("\nET Bellman:\n")
interval_t_b = stats.t.ppf(1.0 - (CI / 2.0),t_b.size-1) * (np.std(t_b)/ np.sqrt(t_b.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(t_b),4))+"+-"+str(round(interval_t_b,4)))
print("\nEFP Bellman:\n")
interval_fp_b = stats.t.ppf(1.0 - (CI / 2.0),fp_b.size-1) * (np.std(fp_b)/ np.sqrt(fp_b.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(fp_b),4))+"+-"+str(round(interval_fp_b,4)))
print("\nEFN Bellman:\n")
interval_fn_b = stats.t.ppf(1.0 - (CI / 2.0),fn_b.size-1) * (np.std(fn_b)/ np.sqrt(fn_b.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(fn_b),4))+"+-"+str(round(interval_fn_b,4)))
print("\nPartitions Bellman:\n")
interval_part_b = stats.t.ppf(1.0 - (CI / 2.0),part_b.size-1) * (np.std(part_b)/ np.sqrt(part_b.size))
# ci = (a.mean() - interval, a.mean() + interval)
print(str(round(np.mean(part_b),4))+"+-"+str(round(interval_part_b,4)))