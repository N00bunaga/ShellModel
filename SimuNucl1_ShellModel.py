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


def CFP_3Cuerpos(j,J): #j el m. angular para particulas identicas, J el resultado
    from sympy.physics.wigner import wigner_6j
    J = np.atleast_1d(J)  # asegura array
    J12 = np.arange(0,2*j,2) #Estados intermedios o progenitor (jj)
    J12P = np.random.choice(J12) #Generador escogido aleatoriamente
    CFP_as = np.zeros(shape=(len(J12),len(J)), dtype=float) #Genero memoria
    #Calculo de los 6j (como matriz o elemento de)
    w6j = np.zeros((len(J12), len(J)), dtype=float)
    for a, J12_i in enumerate(J12):
        for b, J_i in enumerate(J):
            w6j[a, b] = float(wigner_6j(j, j, J12P, J_i, j, J12_i))
    delta = KDelta(J12P, J12)[:, None]
    factor = 2*np.sqrt(2*J12P + 1) * np.sqrt(2*J12)[:, None]
    ###
    CFP_as = delta + factor*w6j

    #Normalización metodo 1 (Analitico)
    CFP_nas1 = np.zeros_like(CFP_as)
    CFP_nas2 = np.zeros_like(CFP_as)
    for b, J_i in enumerate(J):
        #Metodo 1 (Analítico (depende de J))
        w6j_diag = float(wigner_6j(j, j, J12P, J_i, j, J12P))
        denom = 3*(1 + 2*(2*J12P + 1)*w6j_diag)
        if not np.isclose(denom, 0):
            Norm1 = denom**(-0.5)
            CFP_nas1[:, b] = Norm1 * CFP_as[:, b]
        else:
            CFP_nas1[:, b] = np.nan  # estado no físico

        #Metodo 2 (Numerico)
        norm2 = np.sqrt(np.sum(CFP_as[:, b]**2))
        if norm2 > 1e-12:
            CFP_nas2[:, b] = CFP_as[:, b] / norm2
        else:
            CFP_nas2[:, b] = np.nan
    print("Norma analítica:", np.sum(CFP_nas1[:, b]**2))
    print("Norma numérica:", np.sum(CFP_nas2[:, b]**2))
    return [CFP_nas1,"*",CFP_nas2], J12, J12P
E = np.array([3/2,5/2,9/2]) #Los momentos angulares de prueba en el libro
print(CFP_3Cuerpos(5/2,E))
#print(CFP_3Cuerpos(5/2,3/2))