# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""

import numpy
n=5
p=[0.1*(i+1) for i in range(n)]
NSim=20000
alpha=0.6
gamma=1.2
Sp=0.95
FN=0.0
def Se(n,k):
    return 1-Sp*alpha**(float(k)/((n)**gamma))

def FNExact(n):

    prod=1.0

    summ=0.0

    for i in range(n):

        prod=prod*(1-p[i]*(1-(alpha**(1/(n**gamma)))))

        summ=summ+p[i]/(1-p[i]*(1-(alpha**(1/(n**gamma)))))

       

    return (1-Se(1,1))*sum(p)+Se(1,1)*Sp*(alpha**(1/(n**gamma)))*prod*summ

    """prod=1.0

    for i in range(n):

        prod=prod*(1-p[i]*(1-(alpha**(1/(n**gamma)))))

    summ=0.0

    for i in range(n):

        summ=summ+(1-Se(1,1)*(1-Sp*(alpha**(1/(n**gamma)))*prod/(1-p[i]*(1-(alpha**(1/(n**gamma)))))))*p[i]

    return summ"""

   

    

 

for s in range(NSim):

    true=[numpy.random.binomial(1,p[i]) for i in range(n)]

    totalK=sum(true)

    testOutcome=numpy.random.binomial(1,Se(n,totalK))

    if testOutcome==1:

        FN=FN+numpy.random.binomial(totalK,1-Se(1,1))

    else:

        FN=FN+totalK

       

 

print FN/NSim

 

print FNExact(n)

 


 

