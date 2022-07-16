import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interpol

#6. Espectro de Kanai/Tajimi
   
def Kanai_Tajimi(Ap,tipo,duraçao,dt):
    """ 
    Monta o espectro de Kanai/Tajimi e gera um sinal temporal.
        
    :type Ap: int
    :param Ap: Aceleração de pico
    
    :type tipo: string
    :param tipo: define o tipo de solo que a estrutura está ascente
        
    :type duraçao: int
    :param duraçao: tempo de duração do sismo
        
    :type dt: int
    :param dt: Passo temporal adotado para a montagem do sinal
        
    """


        
    ## Ap : Peek Ground Acceleration
    ## tipo : tipo de solo 
    ## duraçao : tempo de duração do sinal
    ## dt : Passo temporal
    
    g = 9.806
    Ap *= g
    pg = 3
    tf = int(duraçao/dt)

    if tipo == 'rocha':
        wg = 8* np.pi
        zg = 0.6
    elif tipo == 'solo_rígido':
        wg = 5* np.pi
        zg = 0.6
    elif tipo == 'solo_mole':
        wg = 2.4* np.pi
        zg = 0.85
        
    f = np.linspace(0,25,tf)
    df = f[1]-f[0]
    w = 2*np.pi*f
    S0 = (Ap**2)/((pg**2)*(np.pi*wg*((1/(2*zg))+2*zg)))
    Sg = S0*((1+4*(zg**2)*(w/wg)**2)/(((1-(w/wg)**2)**2)+4*(zg**2)*(w/wg)**2))


    P = np.random.sample(tf)*2*np.pi
        
   
    ag= np.zeros(tf)
    t = np.linspace(0,duraçao,tf)
    S = np.zeros(tf)



    for i in range(tf):

        S =np.sqrt(2*Sg*df)*np.cos(w*t[i]+ P)
        ag[i] = sum(S)

   

    ag*= Ap/np.max(abs(ag))  ## Normalização das acelerações
    plt.figure(1, figsize=(8,4)) 
    plt.plot(f,ag,'b')
    plt.xlabel('frequência(Hz)'); plt.ylabel('Aceleração(m/s²)')
    plt.xlim(0,20); plt.ylim(-max(ag)*2,max(ag)*2);plt.title(' Aceleração sísmica')
    plt.grid(True)
#---------------------------------------------------------------------------------
    return t,ag


# 15. Duhamel

def Duhamel(wn,a,t,zn =0.05,u0 =0,v0 =0):
    """
    Resolve a equação de equilíbrio pelo método de Duhamel.
    Necessita da utlizacão de Modal_Analysis para ser executado para
    modelos com mais de um grau de liberdade.
        
    """
        
       
    
    tf = int(len(t))
    n  = len(wn)
    dt = t[1] - t[0]
        
    U  = np.zeros((n,tf))
    
    for i in range(n):
        
        
        wd = wn[i]*np.sqrt(1- zn**2)
        e = np.exp(zn*wn[i]*t)
        s = np.sin(wd*t)
        c = np.cos(wd*t)
        
        A = dt*np.cumsum(e*c*a)
        B = dt*np.cumsum(e*s*a)
        
        U[i,:] = (u0*c+(v0+u0*zn*wn[i])*s/wd)/e
        
        U[i,:] +=(A*s-B*c)/e/wd
        
    return U

def Pseudo(ag,classe):
    """
    Define um espectro de pseudoacelerações conforme a NBR15421:2006

    """
    index   = np.array([0,1,2,3,4])
    columns = np.array([0,0.10,0.15])
    Ca      = np.array([[0.8,0.8,0.8],
                        [0.8,1,1],
                        [1.2,1.2,1.2],
                        [1.6,1.6,1.5],
                        [2.5,2.5,2.1]])
    Cv      = np.array([[0.8,0.8,0.8],
                        [1,1,1],
                        [1.7,1.7,1.7],
                        [2.4,2.4,2.2],
                        [3.5,3.5,3.4]])
    f1 = interpol.interp2d(columns,index,Ca,kind = 'linear')
    f2 = interpol.interp2d(columns,index,Cv,kind = 'linear')

    ca = f1(ag,classe)
    cv = f2(ag,classe)
    
    T  = np.linspace(0,10,1000)
    Sa =np.zeros(T.size)

    ags0 = ca*ag*9.81
    ags1 = cv*ag*9.81
    for i in range(T.size):
        if T[i] <=  0.08*cv/ca:
            Sa[i] = ags0*(18.75*T[i]*ca/cv+1)
        elif T[i]<= 0.4*cv/ca:
            Sa[i] = 2.5*ags0
        else:
            Sa[i] = ags1/T[i] 

    return Sa,T
