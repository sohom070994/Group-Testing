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
# import random
from tqdm import tqdm
import numpy as np 
import pandas as pd

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


def calc_measofnd(node_list):

    of=0
    of+=weight_ND(node_list[0],pop)
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
        of+=weight_ND(node_list[ele+1],node_list[ele])
        
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
risk_vec = [2.15/100] * pop

#COVID Parameters
spec = 1
alpha= 0.012595
gamma= 0.257138
const_sens = sens(1,1)

#%% 7. SIMULATIONS
results = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'OptOF',
                                 'OptGS']) 
for w1 in tqdm(list(range(0, 12, 1)),desc='lambda loop',position=0,leave=True):
   for w2 in list(range(0,12,1)):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10
           
           print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
           minOF=100000
           group_size=-1
           for gs in tqdm(range(1,pop+1),desc='sim loop',position=0,leave=True):
               node_list=list(range(0,pop,gs))
               node_list.reverse()
               OF=calc_measof(node_list)
               if (OF<minOF):
                   minOF=OF
                   group_size=gs
           results = results.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'OptOF':minOF,
                                     'OptGS':group_size},ignore_index=True)
           
           
   #%% 7. SIMULATIONS ND
resultsnd = pd.DataFrame(columns= ['w1',
                                 'w2',
                                 'OptOF',
                                 'OptGS']) 
for w1 in tqdm(list(range(0, 12, 1)),desc='lambda loop',position=0,leave=True):
   for w2 in list(range(0,12,1)):
       if(w1+w2<=10):
           lambda_1=w1/10
           lambda_2=w2/10
           
           print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
           minOF=100000
           group_size=-1
           for gs in tqdm(range(1,pop+1),desc='sim loop',position=0,leave=True):
               node_list=list(range(0,pop,gs))
               node_list.reverse()
               OF=calc_measofnd(node_list)
               if (OF<minOF):
                   minOF=OF
                   group_size=gs
           resultsnd = resultsnd.append({'w1':lambda_1,
                                     'w2':lambda_2,
                                     'OptOF':minOF,
                                     'OptGS':group_size},ignore_index=True)
#%%
lambda_1=0.3
lambda_2=0.3
minOF=1000000
group_size=-1
OF_array=np.empty(pop)
for gs in range(1,pop+1):
    print(gs)
    node_list=list(range(0,pop,gs))
    node_list.reverse()
    print(node_list)
    OF=calc_measof(node_list)
    OF_array[gs-1]=OF
    print(OF)
    if (OF<minOF):
        minOF=OF
        group_size=gs
OF_array.tofile('OptGS.csv',sep=',',format='%10.5f')   
#%% 8. OUTPUT TO CSV
results.to_csv('COVID-Opt_GS.csv')
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