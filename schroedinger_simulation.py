import numpy as np 
import matplotlib.pyplot as plt 

def phif(j):
    phi=0
    if j==0:
        phi=0
    elif j==N:
        phi=0
    else:
        phi=np.exp(j*kot*1j)*np.exp(-(8*(4*j-N)**2)/N**2)

    return phi

valores=[10,100,200]

probabilidades=[]
for n in valores:
    N=100
    nciclos=10
    lamda=0.7
    h=0.1
    s=0
    ko=2*np.pi*nciclos/(N)

    #Calculo s,k tilda
    kot=ko
    st=0.25*1/(kot**2)

    #Calculo vt y phi
    phi=np.zeros((N+1,n+1),dtype=complex)
    Vt=np.zeros(N)

    Vt[range(int(2*N/5),int(3*N/5))]=lamda*kot*kot

    for j in range(N):
        phi[j][0]=phif(j)

    AO=-2+(2j)/st-Vt
    #Calculo alpha
    alfa=np.zeros(N,dtype=complex)
    for j in range(N-1,0,-1):
        alfa[j-1]=-1/(AO[j]+1*alfa[j])

    #Iteracion
    norma_cuadrado = np.zeros(n)
    norma = np.zeros(n)
    prob = np.zeros(N)
    beta=np.zeros((N,n),dtype=complex)
    chi=np.zeros((N,n),dtype=complex)
    b=np.zeros((N,n),dtype=complex)
    tiempo=[]
    for i in range(n):
        for j in range(N):
            b[j][i]=(4*1j/st)*phi[j][i]

        for j in range(N-1, 0, -1):
            beta[j-1][i] = (1/(AO[j]+1*alfa[j]))*(b[j][i] - 1*beta[j][i])

        for j in range(1, N):
            chi[j][i] = alfa[j-1]*chi[j-1][i]+beta[j-1][i]

        for j in range(N):
            norma_cuadrado[i] += np.abs(phi[j][i])**2
        
        norma[i]=np.sqrt(norma_cuadrado[i])
        for j in range(N):
            phi[j][i+1] = chi[j][i] - phi[j][i]
            
        tiempo.append(i)
    print(norma[0])
        
    espacio = []
    for j in range(N):
        prob[j] = np.abs(phi[j][n]**2) / norma[9]**2
        espacio.append(j)
        
    probabilidades.append(prob)
    
    #Graficamos la probabilidad 
    plt.plot(espacio, probabilidades[-1], label=f'n={n}', )
    

plt.xlabel(f'Espacio')
plt.ylabel('Probabilidad')
plt.plot(espacio,Vt,label=f'Potencial\nlambda={lamda}')
plt.title('Probabilidad en función del espacio')
plt.legend()
plt.grid(True)
plt.show()
