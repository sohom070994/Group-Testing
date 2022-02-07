#%%
import cmath
import math

N=5

p=[0.1,0.2,0.3,0.4,0.5]

pi=math.pi

j=complex(0,1)

def Q(m):
    summ=0.0
    for n in range(1,N+1):
        prod=1.0
        for k in range(0,N):
            prod=prod*(p[k]*cmath.exp(j*2*pi*n/(N+1.0))+(1-p[k]))
        summ=summ+(1.0-cmath.exp(-j*2*pi*n*m/(N+1.0)))/(1.0-cmath.exp(-j*2*pi*n/(N+1.0)))*prod
            
    return 1.0-float(m)/(N+1.0)-summ/(N+1.0)

# print Q(5)

def P(m):
    summ=0.0
    for n in range(0,N+1):
        prod=1.0
        for k in range(0,N):
            prod=prod*(p[k]*cmath.exp(j*2*pi*n/(N+1.0))+(1-p[k]))
        summ=summ+cmath.exp(-j*2*pi*n*m/(N+1.0))*prod            
    return (summ/(N+1.0)).real

print(P(5))

# %%
import math,cmath
pi=math.pi
compuv=complex(0,1)
risk_vec=[0.1,0.2,0.3,0.4,0.5]
def test_prod(Node_1,Node_2):
    prod=1.0
    for k in range(Node_1,Node_2):
        prod=prod*(risk_vec[k]*cmath.exp(compuv*2*pi*2/(Node_2-Node_1+1.0))+(1-risk_vec[k]))
    return prod

# %%
prod_array=[]
prod_array.append((risk_vec[0]*cmath.exp(compuv*2*pi*2/(2.0))+(1-risk_vec[0])))
for i in range(1,len(risk_vec)):
    prod_array.append(risk_vec[i]*cmath.exp(compuv*2*pi*2/(i+1+1.0))+(1-risk_vec[i])*prod_array[i-1])
# %%
