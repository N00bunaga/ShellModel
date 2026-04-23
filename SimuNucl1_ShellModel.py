import numpy as np
from scipy.special import factorial

N = 4# Número de cuerpos en el shell model

def Cl_Go_CAL(j1,j2,J,m1,m2,M):
    #READ ME
    #Clebsh-Gordan calculadora con j1>=j2 si se quiere el caso j1<=j2 (coeficiente <j2,j1,m2,m1|J,M>) hacer C*FSigno
    if j1<j2:
        return print("Seguir el READ ME")
    
    kmax = min(j1+j2-J,j1-m1,j2+m2)
    kmin = max(0,j2-J-m1,j1+m2-J)
    k = np.linspace(kmin,kmax,1)
    C1 = np.sqrt((2*J+1)*factorial(J+j1-j2,exact=False)*factorial(J-j1+j2,exact=False)*factorial(j1+j2-J,exact=False)/factorial(j1+j2+J+1,exact=False))
    C2 = np.sqrt(factorial(J+M,exact=False)*factorial(J-M,exact=False)*factorial(j1+m1,exact=False)*factorial(j1-m1,exact=False)*factorial(j2+m2,exact=False)*factorial(j2-m2,exact=False))
    C3 = np.sum(pow(-1,k)*pow(factorial(k,exact=False)*factorial(j1+j2-J-k,exact=False)*factorial(j1-m1-k,exact=False)*factorial(j2+m2-k,exact=False)*factorial(J-j2+m1+k,exact=False)*factorial(J-j1-m2+k,exact=False),-1))
    C = C1*C2*C3
    FSigno = pow(-1,J-j1-j2)
    if m1+m2 != M:
        return 0

def KDelta(i,j):
    return np.equal(i,j).astype(int)

import numpy as np
from sympy.physics.wigner import wigner_6j

def KDelta(i,j):
    return np.equal(i,j).astype(int)

def CFP_3Cuerpos(j,J):
    J12 = np.arange(0, 2*j+1, 2)
    J12P = 4 #np.random.choice(J12)

    # 6j para cada J12
    w6j_vals = np.array([
        float(wigner_6j(j, j, J12P, J, j, J12_i))
        for J12_i in J12
    ])

    CFP_as = (
        KDelta(J12P, J12)
        + 2*np.sqrt(2*J12P + 1)*np.sqrt(2*J12 + 1)*w6j_vals
    )

    # 6j diagonal
    w6j_diag = float(wigner_6j(j, j, J12P, J, j, J12P))

    # normalización analítica
    Norm = (3*(1 + 2*(2*J12P + 1)*w6j_diag))**(-0.5)

    CFP = Norm * CFP_as

    return CFP, J12, J12P
j = 5/2
J = 3/2
print(CFP_3Cuerpos(j,J))