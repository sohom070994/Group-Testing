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
import time
#%% 2. RISK VECTOR REALIZATION
'''
LIST OF FUNCTIONS 
FOR GENERATING RISK VECTOR
'''

def get_risk_vector():
   # '''
   # fill risk_vec
   # using coin flips
   # sort risk_vec
   # '''
   # for i in range(0,pop):
   #   risk_vec[i]=round(coin_flip()/100,4)

    elements = [0.0198,0.0138,0.046,0.0395,0.0126]
    probabilities = [0.025,0.3715,0.2685,0.0031,0.3319]
    a=np.random.choice(elements, pop, p=probabilities)
    return a

#%% 3. Debugging Block
# abc=get_risk_vector()
# print(abc)
# print(len(abc))
#%% 4. DILUTION CASE PERFORMANCE MEASURES

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
        
    return et,of,efp,efn

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
    return et,of,efp,efn

def ConfInt(a):
    return 1.96*np.std(a)/(np.size(a)**0.5)
#%% 6. ESTABLISH CONSTANT PARAMETERS
#population
pop = 100

#COVID Parameters
spec = 1
alpha= 0.012595
gamma= 0.257138
const_sens = sens(1,1)
multip=(1-const_sens)
#Confidence Interval
CI = 0.05

#Simulations
simulations=100

lambda_2=0


#%% 7. SIMULATIONS
# sys.stdout = open("FINALRESULTS.txt", "w")
start_time=time.time()
results = pd.DataFrame(columns= ['L1',
                                 'TESTS',
                                 'EFN_D',
                                 'EFN_D_HW',
                                 'EFN_ITS',
                                 'EFN_ITS_HW',
                                 'EFN_RITS',
                                 'EFN_RITS_HW']) 
           

           
efn_d=np.empty(simulations)
efn_RB=np.empty(simulations)
efn_high=np.empty(simulations)
numtests=np.empty(simulations)

for w1 in tqdm(list(range(0, 11, 2)),desc='lambda loop',position=0,leave=True):
       lambda_1=w1/10  
       print('W1(EFN) :'+str(lambda_1)+','+'W2(EFP) :'+str(lambda_2)+'\n')
       # if(w1!=9.6)&(w2==0.2):
       #     continue
       # if (w1==9.5)&(w2!=0.25):
       #      continue
       # # if not (((w1==10)&(w2==0))|((w1==0)&(w2==10))):
       # #      continue             
       for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
                       risk_vec=get_risk_vector()
                       uns_risk_vec=risk_vec
                       
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
                       
                        
                       # node_list=g_d.shortestPath(s)
                       #  # part_d[sim]=len(node_list)
                       # tests,*_, efn_d[sim]=calc_meas(node_list)    
                       tests,*_,efn_d[sim]=Bellman(wt_D,20)
                       numtests[sim]=tests
                       tests=int(round(tests))
                       
                       
                       pop_to_test=random.sample(list(uns_risk_vec),tests)
                       pop_to_test_FN=[ele*multip for ele in pop_to_test]
                       efn_RB[sim]=sum(pop_to_test_FN)+(2.23*(100-tests)/100)               
                       #Optimal Scenario
                       
                       
                       
                       pop_to_test_high=list(risk_vec)[-tests:]
                       pop_to_test_high_others=list(risk_vec)[0:pop-tests]
                       pop_to_test_high_FN=[ele*multip for ele in pop_to_test_high]
                       efn_high[sim]=sum(pop_to_test_high_FN)+sum(pop_to_test_high_others) 
                       
                       
                       # wt_D =[[None for i in range(pop+1)]for i in range(pop+1)]
                       # for i in range(0,pop):
                       #    for j in range(i+1,pop+1):
                       #          wt_D[i][j]=weight_D(i,j)            
                
                       # s = 0
                
                       # g_d = Graph(pop+1) 
              
                       # for i in range(0,pop):
                       #     for j in range(i+1,pop+1):
                       #         g_d.addEdge(i,j,wt_D[i][j])
                       
                        
                       # node_list=g_d.shortestPath(s)
                       #  # part_d[sim]=len(node_list)
                       # *_, efn_d[sim]=calc_meas(node_list)
        
        
       results = results.append({'L1':lambda_1,
                                 'TESTS':round(np.mean(numtests),4),
                                 'EFN_D':round(np.mean(efn_d),4),
                                 'EFN_D_HW':ConfInt(efn_d),
                                 'EFN_ITS':round(np.mean(efn_RB),4),
                                 'EFN_ITS_HW':ConfInt(efn_RB),
                                 'EFN_RITS':round(np.mean(efn_high),4),
                                 'EFN_RITS_HW':ConfInt(efn_high),},ignore_index=True)
# print("Time:")
# print(time.time() - start_time) 
#%% 8. OUTPUT TO CSV
results.to_csv('COVIDGTvsIND_B20June21.csv')          
#%% 7. VARYING N SIMULATIONS
lambda_1=0
lambda_2=0
simulations = 100
B=pop
results = pd.DataFrame(columns= ['OF_ND',
                                 'OF_ND_N10',
                                 'OF_ND_N17',
                                 'OF_ND_N18',
                                 'OF_ND_N19',
                                 'OF_ND_N20',
                                 'OF_ND_N21',
                                 'OF_ND_N22',
                                 'OF_ND_N23',
                                 'OF_ND_N24',
                                 'OF_ND_N25',
                                 'OF_ND_N26',
                                 'OF_ND_N27',
                                 'OF_ND_N30',
                                 'OF_ND_N40',    
                                 'OF_ND_N50',    
                                 'OF_ND_N60',   
                                 'OF_ND_N70',                                 
                                 'OF_ND_N80',
                                 'OF_ND_N90',
                                 'PART_ND',
                                 'PART_ND_N10',
                                 'PART_ND_N17',
                                 'PART_ND_N18',
                                 'PART_ND_N19',
                                 'PART_ND_N20',
                                 'PART_ND_N21',
                                 'PART_ND_N22',
                                 'PART_ND_N23',
                                 'PART_ND_N24',
                                 'PART_ND_N25',
                                 'PART_ND_N26',
                                 'PART_ND_N27',
                                 'PART_ND_N30',
                                 'PART_ND_N40',    
                                 'PART_ND_N50',    
                                 'PART_ND_N60',   
                                 'PART_ND_N70',                                 
                                 'PART_ND_N80',
                                 'PART_ND_N90',
                                 'OF_D',
                                 'OF_D_N10',
                                 'OF_D_N17',
                                 'OF_D_N18',
                                 'OF_D_N19',
                                 'OF_D_N20',
                                 'OF_D_N21',
                                 'OF_D_N22',
                                 'OF_D_N23',
                                 'OF_D_N24',
                                 'OF_D_N25',
                                 'OF_D_N26',
                                 'OF_D_N27',
                                 'OF_D_N30',
                                 'OF_D_N40',    
                                 'OF_D_N50',    
                                 'OF_D_N60',   
                                 'OF_D_N70',                                 
                                 'OF_D_N80',
                                 'OF_D_N90',
                                 'PART_D',
                                 'PART_D_N10',
                                 'PART_D_N17',
                                 'PART_D_N18',
                                 'PART_D_N19',
                                 'PART_D_N20',
                                 'PART_D_N21',
                                 'PART_D_N22',
                                 'PART_D_N23',
                                 'PART_D_N24',
                                 'PART_D_N25',
                                 'PART_D_N26',
                                 'PART_D_N27',
                                 'PART_D_N30',
                                 'PART_D_N40',    
                                 'PART_D_N50',    
                                 'PART_D_N60',   
                                 'PART_D_N70',
                                 'PART_D_N80',
                                 'PART_D_N90',
                                 'DIFF',
                                 'DIFF_N10',
                                 'DIFF_N17',
                                 'DIFF_N18',
                                 'DIFF_N19',
                                 'DIFF_N20',
                                 'DIFF_N21',
                                 'DIFF_N22',
                                 'DIFF_N23',
                                 'DIFF_N24',
                                 'DIFF_N25',
                                 'DIFF_N26',
                                 'DIFF_N27',
                                 'DIFF_N30',
                                 'DIFF_N40',    
                                 'DIFF_N50',    
                                 'DIFF_N60',   
                                 'DIFF_N70',
                                 'DIFF_N80',
                                 'DIFF_N90',
                                 'DIFF_N100'
                                 ])

of_nd=np.empty(simulations)
of_nd_n10=np.empty(simulations)
of_nd_n17=np.empty(simulations)
of_nd_n18=np.empty(simulations)
of_nd_n19=np.empty(simulations)
of_nd_n20=np.empty(simulations)
of_nd_n21=np.empty(simulations)
of_nd_n22=np.empty(simulations)
of_nd_n23=np.empty(simulations)
of_nd_n24=np.empty(simulations)
of_nd_n25=np.empty(simulations)
of_nd_n26=np.empty(simulations)
of_nd_n27=np.empty(simulations)
of_nd_n30=np.empty(simulations)
of_nd_n40=np.empty(simulations)
of_nd_n50=np.empty(simulations)
of_nd_n60=np.empty(simulations)
of_nd_n70=np.empty(simulations)
of_nd_n80=np.empty(simulations)
of_nd_n90=np.empty(simulations)

part_nd=np.empty(simulations)
part_nd_n10=np.empty(simulations)
part_nd_n17=np.empty(simulations)
part_nd_n18=np.empty(simulations)
part_nd_n19=np.empty(simulations)
part_nd_n20=np.empty(simulations)
part_nd_n21=np.empty(simulations)
part_nd_n22=np.empty(simulations)
part_nd_n23=np.empty(simulations)
part_nd_n24=np.empty(simulations)
part_nd_n25=np.empty(simulations)
part_nd_n26=np.empty(simulations)
part_nd_n27=np.empty(simulations)
part_nd_n30=np.empty(simulations)
part_nd_n40=np.empty(simulations)
part_nd_n50=np.empty(simulations)
part_nd_n60=np.empty(simulations)
part_nd_n70=np.empty(simulations)
part_nd_n80=np.empty(simulations)
part_nd_n90=np.empty(simulations)

of_d=np.empty(simulations)
of_d_n10=np.empty(simulations)
of_d_n17=np.empty(simulations)
of_d_n18=np.empty(simulations)
of_d_n19=np.empty(simulations)
of_d_n20=np.empty(simulations)
of_d_n21=np.empty(simulations)
of_d_n22=np.empty(simulations)
of_d_n23=np.empty(simulations)
of_d_n24=np.empty(simulations)
of_d_n25=np.empty(simulations)
of_d_n26=np.empty(simulations)
of_d_n27=np.empty(simulations)
of_d_n30=np.empty(simulations)
of_d_n40=np.empty(simulations)
of_d_n50=np.empty(simulations)
of_d_n60=np.empty(simulations)
of_d_n70=np.empty(simulations)
of_d_n80=np.empty(simulations)
of_d_n90=np.empty(simulations)

part_d=np.empty(simulations)
part_d_n10=np.empty(simulations)
part_d_n17=np.empty(simulations)
part_d_n18=np.empty(simulations)
part_d_n19=np.empty(simulations)
part_d_n20=np.empty(simulations)
part_d_n21=np.empty(simulations)
part_d_n22=np.empty(simulations)
part_d_n23=np.empty(simulations)
part_d_n24=np.empty(simulations)
part_d_n25=np.empty(simulations)
part_d_n26=np.empty(simulations)
part_d_n27=np.empty(simulations)
part_d_n30=np.empty(simulations)
part_d_n40=np.empty(simulations)
part_d_n50=np.empty(simulations)
part_d_n60=np.empty(simulations)
part_d_n70=np.empty(simulations)
part_d_n80=np.empty(simulations)
part_d_n90=np.empty(simulations)

diff=np.empty(simulations)
diff_n10=np.empty(simulations)
diff_n17=np.empty(simulations)
diff_n18=np.empty(simulations)
diff_n19=np.empty(simulations)
diff_n20=np.empty(simulations)
diff_n21=np.empty(simulations)
diff_n22=np.empty(simulations)
diff_n23=np.empty(simulations)
diff_n24=np.empty(simulations)
diff_n25=np.empty(simulations)
diff_n26=np.empty(simulations)
diff_n27=np.empty(simulations)
diff_n30=np.empty(simulations)
diff_n40=np.empty(simulations)
diff_n50=np.empty(simulations)
diff_n60=np.empty(simulations)
diff_n70=np.empty(simulations)
diff_n80=np.empty(simulations)
diff_n90=np.empty(simulations)


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
               
               
               


               of_nd_n10[sim],part_nd_n10[sim]=Bellmanmodof(10,wt_ND,B)
               of_nd_n17[sim],part_nd_n17[sim]=Bellmanmodof(17,wt_ND,B)
               of_nd_n18[sim],part_nd_n18[sim]=Bellmanmodof(18,wt_ND,B)
               of_nd_n19[sim],part_nd_n19[sim]=Bellmanmodof(19,wt_ND,B)
               of_nd_n20[sim],part_nd_n20[sim]=Bellmanmodof(20,wt_ND,B)
               of_nd_n21[sim],part_nd_n21[sim]=Bellmanmodof(21,wt_ND,B)
               of_nd_n22[sim],part_nd_n22[sim]=Bellmanmodof(22,wt_ND,B)
               of_nd_n23[sim],part_nd_n23[sim]=Bellmanmodof(23,wt_ND,B)
               of_nd_n24[sim],part_nd_n24[sim]=Bellmanmodof(24,wt_ND,B)
               of_nd_n25[sim],part_nd_n25[sim]=Bellmanmodof(25,wt_ND,B)
               of_nd_n26[sim],part_nd_n26[sim]=Bellmanmodof(26,wt_ND,B)
               of_nd_n27[sim],part_nd_n27[sim]=Bellmanmodof(27,wt_ND,B)
               of_nd_n30[sim],part_nd_n30[sim]=Bellmanmodof(30,wt_ND,B)
               of_nd_n40[sim],part_nd_n40[sim]=Bellmanmodof(40,wt_ND,B)           
               of_nd_n50[sim],part_nd_n50[sim]=Bellmanmodof(50,wt_ND,B) 
               of_nd_n60[sim],part_nd_n60[sim]=Bellmanmodof(60,wt_ND,B)
               of_nd_n70[sim],part_nd_n70[sim]=Bellmanmodof(70,wt_ND,B)
               of_nd_n80[sim],part_nd_n80[sim]=Bellmanmodof(80,wt_ND,B)
               of_nd_n90[sim],part_nd_n90[sim]=Bellmanmodof(90,wt_ND,B)
             
               of_d_n10[sim],part_d_n10[sim]=Bellmanmodof(10,wt_D,B)
               of_d_n17[sim],part_d_n17[sim]=Bellmanmodof(17,wt_D,B)
               of_d_n18[sim],part_d_n18[sim]=Bellmanmodof(18,wt_D,B)
               of_d_n19[sim],part_d_n19[sim]=Bellmanmodof(19,wt_D,B)
               of_d_n20[sim],part_d_n20[sim]=Bellmanmodof(20,wt_D,B)
               of_d_n21[sim],part_d_n21[sim]=Bellmanmodof(21,wt_D,B)
               of_d_n22[sim],part_d_n22[sim]=Bellmanmodof(22,wt_D,B)
               of_d_n23[sim],part_d_n23[sim]=Bellmanmodof(23,wt_D,B)
               of_d_n24[sim],part_d_n24[sim]=Bellmanmodof(24,wt_D,B)
               of_d_n25[sim],part_d_n25[sim]=Bellmanmodof(25,wt_D,B)
               of_d_n26[sim],part_d_n26[sim]=Bellmanmodof(26,wt_D,B)
               of_d_n27[sim],part_d_n27[sim]=Bellmanmodof(27,wt_D,B)
               of_d_n30[sim],part_d_n30[sim]=Bellmanmodof(30,wt_D,B)
               of_d_n40[sim],part_d_n40[sim]=Bellmanmodof(40,wt_D,B)           
               of_d_n50[sim],part_d_n50[sim]=Bellmanmodof(50,wt_D,B) 
               of_d_n60[sim],part_d_n60[sim]=Bellmanmodof(60,wt_D,B)
               of_d_n70[sim],part_d_n70[sim]=Bellmanmodof(70,wt_D,B)
               of_d_n80[sim],part_d_n80[sim]=Bellmanmodof(80,wt_D,B)
               of_d_n90[sim],part_d_n90[sim]=Bellmanmodof(90,wt_D,B)
                
               
               diff[sim]=(of_d[sim]-of_nd[sim])*100/of_nd[sim]
               diff_n10[sim]=(of_d_n10[sim]-of_nd_n10[sim])*100/of_nd_n10[sim]
               diff_n17[sim]=(of_d_n17[sim]-of_nd_n17[sim])*100/of_nd_n17[sim]
               diff_n18[sim]=(of_d_n18[sim]-of_nd_n18[sim])*100/of_nd_n18[sim]
               diff_n19[sim]=(of_d_n19[sim]-of_nd_n19[sim])*100/of_nd_n19[sim]
               diff_n20[sim]=(of_d_n20[sim]-of_nd_n20[sim])*100/of_nd_n20[sim]
               diff_n21[sim]=(of_d_n21[sim]-of_nd_n21[sim])*100/of_nd_n21[sim]
               diff_n22[sim]=(of_d_n22[sim]-of_nd_n22[sim])*100/of_nd_n22[sim]
               diff_n23[sim]=(of_d_n23[sim]-of_nd_n23[sim])*100/of_nd_n23[sim]
               diff_n24[sim]=(of_d_n24[sim]-of_nd_n24[sim])*100/of_nd_n24[sim]
               diff_n25[sim]=(of_d_n25[sim]-of_nd_n25[sim])*100/of_nd_n25[sim]
               diff_n26[sim]=(of_d_n26[sim]-of_nd_n26[sim])*100/of_nd_n26[sim]
               diff_n27[sim]=(of_d_n27[sim]-of_nd_n27[sim])*100/of_nd_n27[sim]
               diff_n30[sim]=(of_d_n30[sim]-of_nd_n30[sim])*100/of_nd_n30[sim]
               diff_n40[sim]=(of_d_n40[sim]-of_nd_n40[sim])*100/of_nd_n40[sim]
               diff_n50[sim]=(of_d_n50[sim]-of_nd_n50[sim])*100/of_nd_n50[sim]
               diff_n60[sim]=(of_d_n60[sim]-of_nd_n60[sim])*100/of_nd_n60[sim]
               diff_n70[sim]=(of_d_n70[sim]-of_nd_n70[sim])*100/of_nd_n70[sim]
               diff_n80[sim]=(of_d_n80[sim]-of_nd_n80[sim])*100/of_nd_n80[sim]
               diff_n90[sim]=(of_d_n90[sim]-of_nd_n90[sim])*100/of_nd_n90[sim]

               
 
              
results = results.append({           
                                 'OF_ND':round(np.mean(of_nd),4),
                                 'OF_ND_N10':round(np.mean(of_nd_n10),4),
                                 'OF_ND_N17':round(np.mean(of_nd_n17),4),
                                 'OF_ND_N18':round(np.mean(of_nd_n18),4),
                                 'OF_ND_N19':round(np.mean(of_nd_n19),4),
                                 'OF_ND_N20':round(np.mean(of_nd_n20),4),
                                 'OF_ND_N21':round(np.mean(of_nd_n21),4),
                                 'OF_ND_N22':round(np.mean(of_nd_n22),4),
                                 'OF_ND_N23':round(np.mean(of_nd_n23),4),
                                 'OF_ND_N24':round(np.mean(of_nd_n24),4),
                                 'OF_ND_N25':round(np.mean(of_nd_n25),4),
                                 'OF_ND_N26':round(np.mean(of_nd_n26),4),
                                 'OF_ND_N27':round(np.mean(of_nd_n27),4),
                                 'OF_ND_N30':round(np.mean(of_nd_n30),4),    
                                 'OF_ND_N40':round(np.mean(of_nd_n40),4),
                                 'OF_ND_N50':round(np.mean(of_nd_n50),4),
                                 'OF_ND_N60':round(np.mean(of_nd_n60),4),
                                 'OF_ND_N70':round(np.mean(of_nd_n70),4),
                                 'OF_ND_N80':round(np.mean(of_nd_n80),4), 
                                 'OF_ND_N90':round(np.mean(of_nd_n90),4),  
                        
                                 'PART_ND':round(np.mean(part_nd),4),
                                 'PART_ND_N10':round(np.mean(part_nd_n10),4),
                                 'PART_ND_N17':round(np.mean(part_nd_n17),4),
                                 'PART_ND_N18':round(np.mean(part_nd_n18),4),
                                 'PART_ND_N19':round(np.mean(part_nd_n19),4),
                                 'PART_ND_N20':round(np.mean(part_nd_n20),4),
                                 'PART_ND_N21':round(np.mean(part_nd_n21),4),
                                 'PART_ND_N22':round(np.mean(part_nd_n22),4),
                                 'PART_ND_N23':round(np.mean(part_nd_n23),4),
                                 'PART_ND_N24':round(np.mean(part_nd_n24),4),
                                 'PART_ND_N25':round(np.mean(part_nd_n25),4),
                                 'PART_ND_N26':round(np.mean(part_nd_n26),4),
                                 'PART_ND_N27':round(np.mean(part_nd_n27),4),
                                 'PART_ND_N30':round(np.mean(part_nd_n30),4),    
                                 'PART_ND_N40':round(np.mean(part_nd_n40),4),
                                 'PART_ND_N50':round(np.mean(part_nd_n50),4),
                                 'PART_ND_N60':round(np.mean(part_nd_n60),4),
                                 'PART_ND_N70':round(np.mean(part_nd_n70),4),
                                 'PART_ND_N80':round(np.mean(part_nd_n80),4), 
                                 'PART_ND_N90':round(np.mean(part_nd_n90),4),  
                                
                                 'OF_D':round(np.mean(of_d),4),
                                 'OF_D_N10':round(np.mean(of_d_n10),4),
                                 'OF_D_N17':round(np.mean(of_d_n17),4),
                                 'OF_D_N18':round(np.mean(of_d_n18),4),
                                 'OF_D_N19':round(np.mean(of_d_n19),4),
                                 'OF_D_N20':round(np.mean(of_d_n20),4),
                                 'OF_D_N21':round(np.mean(of_d_n21),4),
                                 'OF_D_N22':round(np.mean(of_d_n22),4),
                                 'OF_D_N23':round(np.mean(of_d_n23),4),
                                 'OF_D_N24':round(np.mean(of_d_n24),4),
                                 'OF_D_N25':round(np.mean(of_d_n25),4),
                                 'OF_D_N26':round(np.mean(of_d_n26),4),
                                 'OF_D_N27':round(np.mean(of_d_n27),4),
                                 'OF_D_N30':round(np.mean(of_d_n30),4),    
                                 'OF_D_N40':round(np.mean(of_d_n40),4),
                                 'OF_D_N50':round(np.mean(of_d_n50),4),
                                 'OF_D_N60':round(np.mean(of_d_n60),4),
                                 'OF_D_N70':round(np.mean(of_d_n70),4),
                                 'OF_D_N80':round(np.mean(of_d_n80),4), 
                                 'OF_D_N90':round(np.mean(of_d_n90),4),
                                
                                 'PART_D':round(np.mean(part_d),4),
                                 'PART_D_N10':round(np.mean(part_d_n10),4),
                                 'PART_D_N17':round(np.mean(part_d_n17),4),
                                 'PART_D_N18':round(np.mean(part_d_n18),4),
                                 'PART_D_N19':round(np.mean(part_d_n19),4),
                                 'PART_D_N20':round(np.mean(part_d_n20),4),
                                 'PART_D_N21':round(np.mean(part_d_n21),4),
                                 'PART_D_N22':round(np.mean(part_d_n22),4),
                                 'PART_D_N23':round(np.mean(part_d_n23),4),
                                 'PART_D_N24':round(np.mean(part_d_n24),4),
                                 'PART_D_N25':round(np.mean(part_d_n25),4),
                                 'PART_D_N26':round(np.mean(part_d_n26),4),
                                 'PART_D_N27':round(np.mean(part_d_n27),4),
                                 'PART_D_N30':round(np.mean(part_d_n30),4),    
                                 'PART_D_N40':round(np.mean(part_d_n40),4),
                                 'PART_D_N50':round(np.mean(part_d_n50),4),
                                 'PART_D_N60':round(np.mean(part_d_n60),4),
                                 'PART_D_N70':round(np.mean(part_d_n70),4),
                                 'PART_D_N80':round(np.mean(part_d_n80),4), 
                                 'PART_D_N90':round(np.mean(part_d_n90),4),   
                                 
                                 'DIFF':round(np.mean(diff),4),
                                 'DIFF_N10':round(np.mean(diff_n10),4),
                                 'DIFF_N17':round(np.mean(diff_n17),4),
                                 'DIFF_N18':round(np.mean(diff_n18),4),
                                 'DIFF_N19':round(np.mean(diff_n19),4),
                                 'DIFF_N20':round(np.mean(diff_n20),4),
                                 'DIFF_N21':round(np.mean(diff_n21),4),
                                 'DIFF_N22':round(np.mean(diff_n22),4),
                                 'DIFF_N23':round(np.mean(diff_n23),4),
                                 'DIFF_N24':round(np.mean(diff_n24),4),
                                 'DIFF_N25':round(np.mean(diff_n25),4),
                                 'DIFF_N26':round(np.mean(diff_n26),4),
                                 'DIFF_N27':round(np.mean(diff_n27),4),
                                 'DIFF_N30':round(np.mean(diff_n30),4),    
                                 'DIFF_N40':round(np.mean(diff_n40),4),
                                 'DIFF_N50':round(np.mean(diff_n50),4),
                                 'DIFF_N60':round(np.mean(diff_n60),4),
                                 'DIFF_N70':round(np.mean(diff_n70),4),
                                 'DIFF_N80':round(np.mean(diff_n80),4), 
                                 'DIFF_N90':round(np.mean(diff_n90),4)                                 
                                                                  
                                },ignore_index=True)

# #%% 8. VARYING N SIMULATIONS
# lambda_1=0
# lambda_2=0
# simulations = 100
# B=pop
# results = pd.DataFrame(columns= ['OF_ND',
#                                  'OF_ND_N10',
#                                  'OF_ND_N13',
#                                  'OF_ND_N17',
#                                  'OF_ND_N20',
#                                  'OF_ND_N23',
#                                  'OF_ND_N25',
#                                  'OF_ND_N30',
#                                  'OF_ND_N40',    
#                                  'OF_ND_N50',    
#                                  'OF_ND_N60',   
#                                  'OF_ND_N70',                                 
#                                  'OF_ND_N80',
#                                  'OF_ND_N90',
#                                  'PART_ND',
#                                  'PART_ND_N10',
#                                  'PART_ND_N13',
#                                  'PART_ND_N17',
#                                  'PART_ND_N20',
#                                  'PART_ND_N23',
#                                  'PART_ND_N25',
#                                  'PART_ND_N30',
#                                  'PART_ND_N40',    
#                                  'PART_ND_N50',    
#                                  'PART_ND_N60',   
#                                  'PART_ND_N70',                                 
#                                  'PART_ND_N80',
#                                  'PART_ND_N90',
#                                  'OF_D',
#                                  'OF_D_N10',
#                                  'OF_D_N13',
#                                  'OF_D_N17',
#                                  'OF_D_N20',
#                                  'OF_D_N23',
#                                  'OF_D_N25',
#                                  'OF_D_N30',
#                                  'OF_D_N40',    
#                                  'OF_D_N50',    
#                                  'OF_D_N60',   
#                                  'OF_D_N70',                                 
#                                  'OF_D_N80',
#                                  'OF_D_N90',
#                                  'PART_D',
#                                  'PART_D_N10',
#                                  'PART_D_N13',
#                                  'PART_D_N17',
#                                  'PART_D_N20',
#                                  'PART_D_N23',
#                                  'PART_D_N25',
#                                  'PART_D_N30',
#                                  'PART_D_N40',    
#                                  'PART_D_N50',    
#                                  'PART_D_N60',   
#                                  'PART_D_N70',
#                                  'PART_D_N80',
#                                  'PART_D_N90',
#                                  'DIFF',
#                                  'DIFF_N10',
#                                  'DIFF_N13',
#                                  'DIFF_N17',
#                                  'DIFF_N20',
#                                  'DIFF_N23',
#                                  'DIFF_N25',
#                                  'DIFF_N30',
#                                  'DIFF_N40',    
#                                  'DIFF_N50',    
#                                  'DIFF_N60',   
#                                  'DIFF_N70',
#                                  'DIFF_N80',
#                                  'DIFF_N90',
#                                  'DIFF_N100'
#                                  ])

# of_nd=np.empty(simulations)
# of_nd_n10=np.empty(simulations)
# of_nd_n13=np.empty(simulations)
# of_nd_n17=np.empty(simulations)
# of_nd_n20=np.empty(simulations)
# of_nd_n23=np.empty(simulations)
# of_nd_n25=np.empty(simulations)
# of_nd_n30=np.empty(simulations)
# of_nd_n40=np.empty(simulations)
# of_nd_n50=np.empty(simulations)
# of_nd_n60=np.empty(simulations)
# of_nd_n70=np.empty(simulations)
# of_nd_n80=np.empty(simulations)
# of_nd_n90=np.empty(simulations)

# part_nd=np.empty(simulations)
# part_nd_n10=np.empty(simulations)
# part_nd_n13=np.empty(simulations)
# part_nd_n17=np.empty(simulations)
# part_nd_n20=np.empty(simulations)
# part_nd_n23=np.empty(simulations)
# part_nd_n25=np.empty(simulations)
# part_nd_n30=np.empty(simulations)
# part_nd_n40=np.empty(simulations)
# part_nd_n50=np.empty(simulations)
# part_nd_n60=np.empty(simulations)
# part_nd_n70=np.empty(simulations)
# part_nd_n80=np.empty(simulations)
# part_nd_n90=np.empty(simulations)

# of_d=np.empty(simulations)
# of_d_n10=np.empty(simulations)
# of_d_n13=np.empty(simulations)
# of_d_n17=np.empty(simulations)
# of_d_n20=np.empty(simulations)
# of_d_n23=np.empty(simulations)
# of_d_n25=np.empty(simulations)
# of_d_n30=np.empty(simulations)
# of_d_n40=np.empty(simulations)
# of_d_n50=np.empty(simulations)
# of_d_n60=np.empty(simulations)
# of_d_n70=np.empty(simulations)
# of_d_n80=np.empty(simulations)
# of_d_n90=np.empty(simulations)

# part_d=np.empty(simulations)
# part_d_n10=np.empty(simulations)
# part_d_n13=np.empty(simulations)
# part_d_n17=np.empty(simulations)
# part_d_n20=np.empty(simulations)
# part_d_n23=np.empty(simulations)
# part_d_n25=np.empty(simulations)
# part_d_n30=np.empty(simulations)
# part_d_n40=np.empty(simulations)
# part_d_n50=np.empty(simulations)
# part_d_n60=np.empty(simulations)
# part_d_n70=np.empty(simulations)
# part_d_n80=np.empty(simulations)
# part_d_n90=np.empty(simulations)

# diff=np.empty(simulations)
# diff_n10=np.empty(simulations)
# diff_n13=np.empty(simulations)
# diff_n17=np.empty(simulations)
# diff_n20=np.empty(simulations)
# diff_n23=np.empty(simulations)
# diff_n25=np.empty(simulations)
# diff_n30=np.empty(simulations)
# diff_n40=np.empty(simulations)
# diff_n50=np.empty(simulations)
# diff_n60=np.empty(simulations)
# diff_n70=np.empty(simulations)
# diff_n80=np.empty(simulations)
# diff_n90=np.empty(simulations)


# for sim in tqdm(range(0,simulations),desc='sim loop',position=0,leave=True):
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
#                of_nd[sim]=calc_measof(node_list)
        
#                g_d = Graph(pop+1) 
      
#                for i in range(0,pop):
#                    for j in range(i+1,pop+1):
#                        g_d.addEdge(i,j,wt_D[i][j])
                
#                node_list=g_d.shortestPath(s) 
#                part_d[sim]=len(node_list)
#                of_d[sim]=calc_measof(node_list)
               
               
               


#                of_nd_n10[sim],part_nd_n10[sim]=Bellmanmodof(10,wt_ND,B)
#                of_nd_n13[sim],part_nd_n13[sim]=Bellmanmodof(13,wt_ND,B)
#                of_nd_n17[sim],part_nd_n17[sim]=Bellmanmodof(17,wt_ND,B)
#                of_nd_n20[sim],part_nd_n20[sim]=Bellmanmodof(20,wt_ND,B)
#                of_nd_n23[sim],part_nd_n23[sim]=Bellmanmodof(23,wt_ND,B)
#                of_nd_n25[sim],part_nd_n25[sim]=Bellmanmodof(25,wt_ND,B)
#                of_nd_n30[sim],part_nd_n30[sim]=Bellmanmodof(30,wt_ND,B)
#                of_nd_n40[sim],part_nd_n40[sim]=Bellmanmodof(40,wt_ND,B)           
#                of_nd_n50[sim],part_nd_n50[sim]=Bellmanmodof(50,wt_ND,B) 
#                of_nd_n60[sim],part_nd_n60[sim]=Bellmanmodof(60,wt_ND,B)
#                of_nd_n70[sim],part_nd_n70[sim]=Bellmanmodof(70,wt_ND,B)
#                of_nd_n80[sim],part_nd_n80[sim]=Bellmanmodof(80,wt_ND,B)
#                of_nd_n90[sim],part_nd_n90[sim]=Bellmanmodof(90,wt_ND,B)
             
#                of_d_n10[sim],part_d_n10[sim]=Bellmanmodof(10,wt_D,B)
#                of_d_n13[sim],part_d_n13[sim]=Bellmanmodof(13,wt_D,B)
#                of_d_n17[sim],part_d_n17[sim]=Bellmanmodof(17,wt_D,B)
#                of_d_n20[sim],part_d_n20[sim]=Bellmanmodof(20,wt_D,B)
#                of_d_n23[sim],part_d_n23[sim]=Bellmanmodof(23,wt_D,B)
#                of_d_n25[sim],part_d_n25[sim]=Bellmanmodof(25,wt_D,B)
#                of_d_n30[sim],part_d_n30[sim]=Bellmanmodof(30,wt_D,B)
#                of_d_n40[sim],part_d_n40[sim]=Bellmanmodof(40,wt_D,B)           
#                of_d_n50[sim],part_d_n50[sim]=Bellmanmodof(50,wt_D,B) 
#                of_d_n60[sim],part_d_n60[sim]=Bellmanmodof(60,wt_D,B)
#                of_d_n70[sim],part_d_n70[sim]=Bellmanmodof(70,wt_D,B)
#                of_d_n80[sim],part_d_n80[sim]=Bellmanmodof(80,wt_D,B)
#                of_d_n90[sim],part_d_n90[sim]=Bellmanmodof(90,wt_D,B)
                
               
#                diff=(of_d[sim]-of_nd[sim])*100/of_nd[sim]
#                diff_n10=(of_d_n10[sim]-of_nd_n10[sim])*100/of_nd_n10[sim]
#                diff_n13=(of_d_n13[sim]-of_nd_n13[sim])*100/of_nd_n13[sim]
#                diff_n17=(of_d_n17[sim]-of_nd_n17[sim])*100/of_nd_n17[sim]
#                diff_n20=(of_d_n20[sim]-of_nd_n20[sim])*100/of_nd_n20[sim]
#                diff_n23=(of_d_n23[sim]-of_nd_n23[sim])*100/of_nd_n23[sim]
#                diff_n25=(of_d_n25[sim]-of_nd_n25[sim])*100/of_nd_n25[sim]
#                diff_n30=(of_d_n30[sim]-of_nd_n30[sim])*100/of_nd_n30[sim]
#                diff_n40=(of_d_n40[sim]-of_nd_n40[sim])*100/of_nd_n40[sim]
#                diff_n50=(of_d_n50[sim]-of_nd_n50[sim])*100/of_nd_n50[sim]
#                diff_n60=(of_d_n60[sim]-of_nd_n60[sim])*100/of_nd_n60[sim]
#                diff_n70=(of_d_n70[sim]-of_nd_n70[sim])*100/of_nd_n70[sim]
#                diff_n80=(of_d_n80[sim]-of_nd_n80[sim])*100/of_nd_n80[sim]
#                diff_n90=(of_d_n90[sim]-of_nd_n90[sim])*100/of_nd_n90[sim]

               
 
              
# results = results.append({           
#                                  'OF_ND':round(np.mean(of_nd),4),
#                                  'OF_ND_N10':round(np.mean(of_nd_n10),4),
#                                  'OF_ND_N13':round(np.mean(of_nd_n13),4),
#                                  'OF_ND_N17':round(np.mean(of_nd_n17),4),
#                                  'OF_ND_N20':round(np.mean(of_nd_n20),4),
#                                  'OF_ND_N23':round(np.mean(of_nd_n23),4),
#                                  'OF_ND_N25':round(np.mean(of_nd_n25),4),
#                                  'OF_ND_N30':round(np.mean(of_nd_n30),4),    
#                                  'OF_ND_N40':round(np.mean(of_nd_n40),4),
#                                  'OF_ND_N50':round(np.mean(of_nd_n50),4),
#                                  'OF_ND_N60':round(np.mean(of_nd_n60),4),
#                                  'OF_ND_N70':round(np.mean(of_nd_n70),4),
#                                  'OF_ND_N80':round(np.mean(of_nd_n80),4), 
#                                  'OF_ND_N90':round(np.mean(of_nd_n90),4),  
                        
#                                  'PART_ND':round(np.mean(part_nd),4),
#                                  'PART_ND_N10':round(np.mean(part_nd_n10),4),
#                                  'PART_ND_N13':round(np.mean(part_nd_n13),4),
#                                  'PART_ND_N17':round(np.mean(part_nd_n17),4),
#                                  'PART_ND_N20':round(np.mean(part_nd_n20),4),
#                                  'PART_ND_N23':round(np.mean(part_nd_n23),4),
#                                  'PART_ND_N25':round(np.mean(part_nd_n25),4),
#                                  'PART_ND_N30':round(np.mean(part_nd_n30),4),    
#                                  'PART_ND_N40':round(np.mean(part_nd_n40),4),
#                                  'PART_ND_N50':round(np.mean(part_nd_n50),4),
#                                  'PART_ND_N60':round(np.mean(part_nd_n60),4),
#                                  'PART_ND_N70':round(np.mean(part_nd_n70),4),
#                                  'PART_ND_N80':round(np.mean(part_nd_n80),4), 
#                                  'PART_ND_N90':round(np.mean(part_nd_n90),4),  
                                
#                                  'OF_D':round(np.mean(of_d),4),
#                                  'OF_D_N10':round(np.mean(of_d_n10),4),
#                                  'OF_D_N13':round(np.mean(of_d_n13),4),
#                                  'OF_D_N17':round(np.mean(of_d_n17),4),
#                                  'OF_D_N20':round(np.mean(of_d_n20),4),
#                                  'OF_D_N23':round(np.mean(of_d_n23),4),
#                                  'OF_D_N25':round(np.mean(of_d_n25),4),
#                                  'OF_D_N30':round(np.mean(of_d_n30),4),    
#                                  'OF_D_N40':round(np.mean(of_d_n40),4),
#                                  'OF_D_N50':round(np.mean(of_d_n50),4),
#                                  'OF_D_N60':round(np.mean(of_d_n60),4),
#                                  'OF_D_N70':round(np.mean(of_d_n70),4),
#                                  'OF_D_N80':round(np.mean(of_d_n80),4), 
#                                  'OF_D_N90':round(np.mean(of_d_n90),4),
                                
#                                  'PART_D':round(np.mean(part_d),4),
#                                  'PART_D_N10':round(np.mean(part_d_n10),4),
#                                  'PART_D_N13':round(np.mean(part_d_n13),4),
#                                  'PART_D_N17':round(np.mean(part_d_n17),4),
#                                  'PART_D_N20':round(np.mean(part_d_n20),4),
#                                  'PART_D_N23':round(np.mean(part_d_n23),4),
#                                  'PART_D_N25':round(np.mean(part_d_n25),4),
#                                  'PART_D_N30':round(np.mean(part_d_n30),4),    
#                                  'PART_D_N40':round(np.mean(part_d_n40),4),
#                                  'PART_D_N50':round(np.mean(part_d_n50),4),
#                                  'PART_D_N60':round(np.mean(part_d_n60),4),
#                                  'PART_D_N70':round(np.mean(part_d_n70),4),
#                                  'PART_D_N80':round(np.mean(part_d_n80),4), 
#                                  'PART_D_N90':round(np.mean(part_d_n90),4),   
                                 
#                                  'DIFF':round(np.mean(diff),4),
#                                  'DIFF_N10':round(np.mean(diff_n10),4),
#                                  'DIFF_N13':round(np.mean(diff_n13),4),
#                                  'DIFF_N17':round(np.mean(diff_n17),4),
#                                  'DIFF_N20':round(np.mean(diff_n20),4),  
#                                  'DIFF_N23':round(np.mean(diff_n23),4),
#                                  'DIFF_N25':round(np.mean(diff_n25),4),
#                                  'DIFF_N30':round(np.mean(diff_n30),4),    
#                                  'DIFF_N40':round(np.mean(diff_n40),4),
#                                  'DIFF_N50':round(np.mean(diff_n50),4),
#                                  'DIFF_N60':round(np.mean(diff_n60),4),
#                                  'DIFF_N70':round(np.mean(diff_n70),4),
#                                  'DIFF_N80':round(np.mean(diff_n80),4), 
#                                  'DIFF_N90':round(np.mean(diff_n90),4)                                 
                                                                  
#                                 },ignore_index=True)


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
results.to_csv('NVar_001_Bpop_100extra22222.csv')
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