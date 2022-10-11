##This code calculates the eigenvalues of the full harmonic oscillator using Hermite polynomials as orthogonal functions.
##As the number of eigenvalues are increased, the roots of the polynomial should also be increased.
from sympy import Symbol
import sympy as s
import numpy as np
x=Symbol('x')
from math import sqrt, factorial, pi

m=6      ##This many eigenvalues will be found
n=200    ##This many roots of the Legendre polynomial will be found, larger the better
a=-6     ##initial value of integration
b=6      ##Final value of integration

def Poly(x):               ##This gives the nth order of the Polynomial in Polynomial form. used to find the derivative of the Polynomial
    P0=1                   ##first order Legendre polynomial
    P1=x                   ##second order Legendre Polynomial
    for i in range(1,n):
       P2=(1/(i+1))*((2*i+1)*x*P1-i*P0)
       P0=P1
       P1=P2
    return P2

def PolyDer(x):       #This gives the first derivative of the nth Polynomial in Polynomial form. Used to find the weights
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

def Roots():                                 #We now find the roots of the Legendre Polynomial using Bisection method
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
        z=((b-a)/2)*rootsinlist[i]+(b+a)/2
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

def weights():    #weights of the Legendre Polynomial. Used in the integration
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

def WavL(): ##the left point of the wavefunction, to calculate the derivative of the wavefunction
    GNP=[]
    for j in range (0,n):
        H0=1
        G0=(1/pi**(1/4))*H0*s.exp(-(X0[j])**2/2)
        H1=2*X0[j]
        G1=(1/(sqrt(2)*pi**(1/4)))*H1*s.exp(-(X0[j])**2/2)
        GN=[G0,G1]
        for i in range(1,m-1):
            H2=2*X0[j]*H1-2*i*H0
            G2=(1/(sqrt(2**(i+1)*factorial(i+1))*pi**(1/4)))*H2*s.exp(-(X0[j])**2/2)
            GN.append(G2)
            H0=H1
            H1=H2
        GNP.append(GN)
    return GNP

WavLL=WavL()
print('left point WF in concatenated',WavLL)

def WavSL():          #This gives wavefunctions in a single list
    WE=[]
    for i in range(0,n):
        j=0
        while j<=m-1:
            z=WavLL[i][j]
            WE.append(z)
            j=j+1
    return WE

WavSLL=WavSL()
print('left point WF single list', WavSLL)

def RX2():                              #Gives (root+h) used for derivative
    x0=[]
    for i in range(0,n):
        xL=X[i]+h
        x0.append(xL)
    return x0
X2=RX2()
print('root+h',X2)

def WavR(): #the right point of the wavefunction, used to calculate the derivative of the wavefunction
    GNP=[]
    for j in range (0,n):
        H0=1
        G0=(1/pi**(1/4))*H0*s.exp(-(X2[j])**2/2)
        H1=2*X2[j]
        G1=(1/(sqrt(2)*pi**(1/4)))*H1*s.exp(-(X2[j])**2/2)
        GN=[G0,G1]
        for i in range(1,m-1):
            H2=2*X2[j]*H1-2*i*H0
            G2=(1/(sqrt(2**(i+1)*factorial(i+1))*pi**(1/4)))*H2*s.exp(-(X2[j])**2/2)
            GN.append(G2)
            H0=H1
            H1=H2
        GNP.append(GN)
    return GNP

WavRL=WavR()
print('right point WF in concatenated',WavRL)

def WavSR():          #This gives wavefunctions in a single list
    WE=[]
    for i in range(0,n):
        j=0
        while j<=m-1:
            z=WavRL[i][j]
            WE.append(z)
            j=j+1
    return WE

WavSLR=WavSR()
print('right point WF single list', WavSLR)

def derwav():  #the derivative of the wavefunction
    Der=[]
    for i in range(0,len(WavSLR)):
        D=(WavSLR[i]-WavSLL[i])/(2*h)
        Der.append(D)
    return Der

DWIL=derwav()
print('The derivatives are:::',DWIL)

def Amn():                         #We find the value of Amn
    Add=[]
    for i in range(0,m):
        for j in range (0,m):
            h=i
            k=j
            C=0
            Su=[]
            while C<=(n-1):
                S=DWIL[h]*DWIL[k]*weightsinlist[C]*((b-a)/2)
                Su.append(S)
                h=h+m
                k=k+m
                C=C+1
            A=(sum(Su))/2
            Add.append(A)
    return(Add)

AmnL=Amn()
print('The values of Amn are:::',AmnL)

######################Calculation of Bmn###############################################

def hermitepoly():
    GNP=[]
    for j in range (0,n):
        H0=1
        G0=(1/pi**(1/4))*H0*s.exp(-(X[j])**2/2)
        H1=2*X[j]
        G1=(1/(sqrt(2)*pi**(1/4)))*H1*s.exp(-(X[j])**2/2)
        GN=[G0,G1]
        for i in range(1,m-1):
            H2=2*X[j]*H1-2*i*H0
            G2=(1/(sqrt(2**(i+1)*factorial(i+1))*pi**(1/4)))*H2*s.exp(-(X[j])**2/2)
            GN.append(G2)
            H0=H1
            H1=H2
        GNP.append(GN)
    return GNP

WIL=hermitepoly()
print('The wavefunctions in concatenated list',WIL)

def WF():          #This gives wavefunctions in a single list
    WE=[]
    for i in range(0,n):
        j=0
        while j<=m-1:
            z=WIL[i][j]
            WE.append(z)
            j=j+1
    return WE

WFIL=WF()
print('The wavefunctions in a list', WFIL)

def Bmn():
    Add=[]
    for i in range(0,m):
        for j in range (0,m):
            h=i
            k=j
            C=0
            Su=[]
            while C<=(n-1):
                S=WFIL[h]*WFIL[k]*(X[C]*X[C]*0.5)*weightsinlist[C]*((b-a)/2)   ####This is where the potential is used####
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

print(Mat())

def eig():
    A=Mat()
    e=np.linalg.eigvals(A)
    return e

val=eig()

val.sort()

print('the eigen values are:::',val)

