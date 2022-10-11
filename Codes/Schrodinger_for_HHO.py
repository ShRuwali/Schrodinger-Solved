##This code calculates the eigenvalues of the full and the half harmonic oscillator depending on where the potential
##has been kept. The orthogonal wavefunction that is being used over here is the sine functions.

from sympy import Symbol
import numpy as np
x=Symbol('x')
from math import sqrt, pi, sin
m=30                ## This many eigenvalues will be created. Also a 30 by 30 matrix
n=600              ## This many roots will be found
I=2                ## Initial value of integration
F=3                ## Final value of integration
a=1                ## Width of the potential
sh=4 ## Shifting the potential so that it originates from the initial value of the integration. If placed between the integration will yield eigenvalues of full harmonic oscillator
w=100              ##The value of omega for the potential

def Poly(x):                                        ##This gives the nth order of the Polynomial in Polynomial form
    P0=1
    P1=x
    for i in range(1,n):
       P2=(1/(i+1))*((2*i+1)*x*P1-i*P0)
       P0=P1
       P1=P2
    return P2

def PolyDer(x):                  #This gives the first derivative of the nth Polynomial in Polynomial form
    P0=1
    P1=x
    for i in range(1,n):
        P2=(1/(i+1))*((2*i+1)*x*P1-i*P0)
        P0=P1
        PL=P1
        P1=P2
    P3=(1/(n+1))*((2*n+1)*x*P2-n*PL)
    PND=((n+1)/(1-x **2))*(x*P2-P3)
    return PND

def Roots():                                 #We now find the roots of the Legendre Polynomial
    Root=[]
    Rp = -1
    while Rp <= 1:
        Lp=Rp
        Rp=Rp+0.00001
        if Poly(Lp)*Poly(Rp)<0:
            LP=Lp
            RP=Rp
            Mp=(RP+LP)/2
            while abs(LP-RP)>0.00001:
                if Poly(Mp)*Poly(LP)>0:
                    LP=Mp
                else:
                    RP=Mp
                Mp=(LP+RP)/2
            Root.append(Mp)
    return Root

rootsinlist=Roots()                                 #Roots in list form
print('roots',rootsinlist)

def XV():                                          #Scaling so that integration is between desired values
    x=[]
    for i in range(0,n):
        z=((F-I)/2)*rootsinlist[i]+(F+I)/2
        x.append(z)
    return x
X=XV()
print('wav val',X)

def Derivativ():                                     #Derivatives to find the weights
    Deri=[]
    for i in range (0,n):
        z = PolyDer(rootsinlist[i])
        Deri.append(z)
    return Deri

derivativesinlist=Derivativ()                        #Derivatives in list form

def weights():
    weight=[]
    for i in range (0,n):
        z=2/((1-rootsinlist[i]**2)*(derivativesinlist[i]**2))
        Z=z
        weight.append(Z)
    return weight

weightsinlist=weights()                               #Weights in list form
print('weights',weightsinlist)

h=0.00001                                #The step in calculating derivative

def RX0():                              #Gives (root-h) used for derivative
    x0=[]
    for i in range(0,n):
        xL=X[i]-h
        x0.append(xL)
    return x0
X0=RX0()
print('root-h',X0)

def RX2():                              #Gives (root+h) used for derivative
    x0=[]
    for i in range(0,n):
        xL=X[i]+h
        x0.append(xL)
    return x0
X2=RX2()
print('root+h',X2)

def WavL():
    D=[]
    for j in range(0,n):
        for i in range(1,m+1):
            V=sqrt(2/a)*sin((i*pi*X0[j])/a)
            D.append(V)
    return D

WavLL=WavL()
print('WF in left',WavLL)

def WavR():
    D=[]
    for j in range(0,n):
        for i in range(1,m+1):
            V=sqrt(2/a)*sin((i*pi*X2[j])/a)
            D.append(V)
    return D

WavRL = WavR()
print('WF in right', WavRL)

def derwav():
    Der=[]
    for i in range(0,len(WavRL)):
        D=(WavRL[i]-WavLL[i])/(2*h)
        Der.append(D)
    return Der

DWIL=derwav()
print('Derivative',DWIL)

def Amn():                         #We find the value of Amn
    Add=[]
    for i in range(0,m):
        for j in range (0,m):
            h=i
            k=j
            C=0
            Su=[]
            while C<=(n-1):
                S=DWIL[h]*DWIL[k]*weightsinlist[C]*((F-I)/2)
                Su.append(S)
                h=h+m
                k=k+m
                C=C+1
            A=(sum(Su))/2
            Add.append(A)
    return(Add)

AmnL=Amn()
print('Amn',AmnL)

#####################Calculation of Bmn##############################

def WFL():
    D=[]
    for j in range(0,n):
        for i in range(1,m+1):
            V=sqrt(2/a)*sin((i*pi*X[j])/a)
            D.append(V)
    return D

WFIL=WFL()
print('WF at point',WFIL)

def Bmn():
    Add=[]
    for i in range(0,m):
        for j in range (0,m):
            h=i
            k=j
            C=0
            Su=[]
            while C<=(n-1):
                S=WFIL[h]*WFIL[k]*     (0.5*(w**2)*(X[C]-sh/2)**2)*      weightsinlist[C]*((F-I)/2)   ####This is where the potential is used####
                Su.append(S)
                h=h+m
                k=k+m
                C=C+1
            A=sum(Su)
            Add.append(A)
    return(Add)

BmnL=Bmn()
print('The values of Bmn are:::',BmnL)

def Fmn():
    y=[]
    for i in range(0,len(AmnL)):
        S=AmnL[i]+BmnL[i]
        y.append(S)
    return y

FmnL=Fmn()
print('The value of Fmn',FmnL)

def Mat():
    M=np.zeros((m,m))
    k=0
    for i in range(0,m):
        for j in range(0,m):
            M[i][j]=FmnL[k]
            k=k+1
    return M

print('The eigen matrix',Mat())

def eig():
    A=Mat()
    e=np.linalg.eigvals(A)
    return e

val=eig()/w

val.sort()

print('the eigen values are:::',val)

