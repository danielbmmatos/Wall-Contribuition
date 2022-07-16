# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:25:16 2019

@author: Daniel Matos
"""

import numpy as np
import scipy.linalg as sc
import matplotlib.pyplot as plt
import pandas as pd
from pandas import ExcelWriter
from numba import jit
import os
#-----------------------------------------------------------------------------------
#1. Inicializando a classe
#-----------------------------------------------------------------------------------
""" 
Esse programa tem como finalidade a criação de uma classe para a disciplina
Dinâmica de sistemas estruturais/UFRGS

"""

class Vibmec():
    """ 
    A classe Vibmec será utilizada para seguir o paradigma da orientação
    a objetos
    
    """
    
    def __init__(self,Nome =0,Nn=0,E=0,Kr = 0,Mr = 0):
        
        """
        Inicializa a classe a partir de um arquivo .xlsx
    
        :type Nome: String
        :param Nome: Planilha para importação dos dados
        
        :type Nn: int
        :param Nn: Número de nós do problema.
        
        """
        if hasattr(Kr,'shape') == False:
#-------------------------------------------------------------------------    
#1.1. Importar planilha do Excel  
#-------------------------------------------------------------------------  
            Arquivo = pd.read_excel(Nome)
#-------------------------------------------------------------------------     
#1.2. Montar vetor de coordenadas ( X e Y de cada nó) 
#------------------------------------------------------------------------- 
            cx = list(Arquivo['Cx'])[0:Nn]
            cy = list(Arquivo['Cy'])[0:Nn]
        
#-------------------------------------------------------------------------     
#1.3. Montar matriz identidade( Nó inicial e final de cada barra)
#-------------------------------------------------------------------------  
            Id1     = list(Arquivo['barra (nó 1)'].dropna())
            Id2     = list(Arquivo['barra (nó 2)'].dropna())
            Nb      = len(Id1)
            ID      = np.zeros((2,Nb))
            ID[0,:] = Id1
            ID[1,:] = Id2
#-------------------------------------------------------------------------     
#1.4. Alocar as propriedades de cada barra
#-------------------------------------------------------------------------  
            A   = list(Arquivo['Area(m2)'])
            I   = list(Arquivo['Inércia(m4)'])
            RHO = list(Arquivo['Densidade'])
          
#------------------------------------------------------------------------- 
#1.5. Matriz identidade em relação aos graus de liberdade
#-------------------------------------------------------------------------  
            IDG = np.zeros((6,Nb))
    
            for i in range(3):
        
                IDG[i,:]   = ID[0,:]*3-2+i
                IDG[i+3,:] = ID[1,:]*3-2+i
#-------------------------------------------------------------------------        
#1.6. Comprimento de cada barra e cossenos diretores
#-------------------------------------------------------------------------  
            Lx   = np.zeros(Nb)
            Ly   = np.zeros(Nb)
            cosx = np.zeros(Nb)
            cosy = np.zeros(Nb)
            L    = np.zeros(Nb)
        
            for n in range (Nb):
    
                k1      = int(ID[0,n] -1)  # 'k' são os buscadores da matriz Id
                k2      = int(ID[1,n] -1)
                Lx[n]   = cx[k2] - cx[k1]
                Ly[n]   = cy[k2] - cy[k1]
                L[n]    = np.sqrt(Lx[n]**2 + Ly[n]**2)
                cosx[n] = Lx[n]/L[n]
                cosy[n] = Ly[n]/L[n]
    
#-------------------------------------------------------------------------        
#1.7. Montagem das matrizes de massa e rigidez
#-------------------------------------------------------------------------  

            K = np.zeros((Nn*3,Nn*3))
            M = np.zeros((Nn*3,Nn*3))

            for i in range (Nb):
    
    #1.7.1 Matriz de rigidez local da barra
                Ke =np.array([[E[i]*A[i]/L[i], 0, 0, -E[i]*A[i]/L[i],0 ,0 ],
                  [0, 12*E[i]*I[i]/(L[i]**3), 6*E[i]*I[i]/(L[i]**2), 0,
                   -12*E[i]*I[i]/(L[i]**3),6*E[i]*I[i]/(L[i]**2)],
                  [0,6*E[i]*I[i]/(L[i]**2), 4*E[i]*I[i]/L[i], 0, 
                   -6*E[i]*I[i]/(L[i]**2), 2*E[i]*I[i]/L[i] ],
                  [-E[i]*A[i]/L[i], 0, 0, E[i]*A[i]/L[i],0 ,0 ],
                  [0, -12*E[i]*I[i]/(L[i]**3), -6*E[i]*I[i]/(L[i]**2),
                   0,12*E[i]*I[i]/(L[i]**3),-6*E[i]*I[i]/(L[i]**2)],
                  [0,6*E[i]*I[i]/(L[i]**2), 2*E[i]*I[i]/L[i], 0,
                   -6*E[i]*I[i]/(L[i]**2), 4*E[i]*I[i]/L[i] ]])

    #1.7.2 Matriz de massa local da barra
                Me = ((RHO[i]*A[i]*L[i])/420)*np.array([[140, 0, 0, 70, 0, 0],
         [0, 156, 22*L[i], 0, 54, -13*L[i]],
         [0, 22*L[i], 4*(L[i]**2), 0, 13*L[i], -3*(L[i]**2)],
         [70, 0, 0, 140, 0, 0],
         [0, 54, 13*L[i], 0, 156, -22*L[i]],
         [0, -13*L[i], -3*(L[i]**2), 0, -22*L[i], 4*(L[i]**2)]])
   
    #1.7.3 Matriz de rotação 
                R = np.array([[cosx[i], cosy[i], 0, 0 ,0 ,0],
                  [-cosy[i], cosx[i],0, 0, 0, 0],
                  [0,0,1,0,0,0],                     
                  [0,0,0,cosx[i], cosy[i], 0],
                  [0, 0, 0,-cosy[i], cosx[i],0],
                  [0,0,0,0,0,1]])
        
    #1.7.4 Rotação das matrizes para coordendas globais
                     
                KT = np.dot(np.dot(R.T, Ke),R)             
                MT = np.dot(np.dot(R.T, Me),R)
    
    #1.7.5 Matrizes temporárias
                k_temp1 = np.zeros((Nn*3,Nn*3))    
                m_temp1 = np.zeros((Nn*3,Nn*3))
        
   #1.7.6 Alocação das matrizes temporárias na matriz global
                m = int(IDG[0,i]-1)
                n = int(IDG[2,i])
                o = int(IDG[3,i]-1)
                p = int(IDG[5,i])
    
                k_temp1[m:n,m:n] = KT[0:3,0:3]
                k_temp1[o:p,m:n] = KT[3:6,0:3]
                k_temp1[m:n,o:p] = KT[0:3,3:6]
                k_temp1[o:p,o:p] = KT[3:6,3:6]
    
                K += k_temp1 
    
                m_temp1[m:n,m:n] = MT[0:3,0:3]
                m_temp1[o:p,m:n] = MT[3:6,0:3]
                m_temp1[m:n,o:p] = MT[0:3,3:6]
                m_temp1[o:p,o:p] = MT[3:6,3:6]
    
                M += m_temp1 
#-------------------------------------------------------------------------
# Atributos da classe.
            self.Stiff  = K
            self.Mass   = M
            self.xcoord = cx
            self.ycoord = cy
            self.ID     = ID
            self.L      = L
            
        else:
            self.Stiff  = Kr
            self.Mass   = Mr
#-------------------------------------------------------------------------
#2. Desenhar a estrutura
#-------------------------------------------------------------------------
    def DRAW(self):
        """
        Desenha a estrutura com base nas coordenadas disponibilizadas na
        planilha
        
        """
        cx  = self.xcoord
        cy  = self.ycoord
        ID  = self.ID
        CX  = np.zeros(2)
        CY  = np.zeros(2)
        fig = plt.figure(111,figsize=(8,10))
        ax  = fig.add_subplot(111)
        plt.xlabel('x (m)');plt.ylabel('y (m)')
        #plt.grid(True)
        Nb  = len(ID[0,:])
        for n in range (Nb):
            k1      = int(ID[0,n] -1) 
            k2      = int(ID[1,n] -1)
            CX[0]   = cx[k1] 
            CX[1]   = cx[k2]
            CY[0]   = cy[k1] 
            CY[1]   = cy[k2] 
            ax.plot(CX,CY,'k',linewidth=0.6)
            if CX[0]==CX[1]:
                ax.text((CX[0]+CX[1])/2,(CY[0]+CY[1])/2,n, horizontalalignment='left',color = 'red')
            else:
                ax.text((CX[0]+CX[1])/2,(CY[0]+CY[1])/2,n, horizontalalignment='center',color = 'red')
                
            ax.plot(cx[k1],cy[k1],'ko')
            ax.text(cx[k1]+0.5,cy[k1]+0.5,k1, horizontalalignment='center',color = 'black')
            ax.plot(cx[k2],cy[k2],'ko')
            ax.text(cx[k2]+0.5,cy[k2]+0.5,k2, horizontalalignment='center',color = 'black')
            fig.canvas.draw()

#--------------------------------------------------------------
#3. Eliminação dos graus de liberdade restritos
#--------------------------------------------------------------
  
    def Restr(self,Nr): 
        """ 
        Restringe os graus de liberdade a partir de uma lista de
        graus restritos
        
        :type Nr: list
        :param Nr: Lista com os graus de liberdade restritos
        
        """
        
        K = self.Stiff
        M = self.Mass

        Kr_1 = np.delete(K,Nr,0)
        Kr  = np.delete(Kr_1,Nr,1)
        Mr_1 = np.delete(M,Nr,0)
        Mr  = np.delete(Mr_1,Nr,1)

        df = pd.DataFrame(Kr)
        writer = ExcelWriter('Matriz de rigidez.xlsx')
        df.to_excel(writer,'Sheet1', index=False) 
        writer.save()
                               ##Exportação das matrizes de massa e rigidez
        df1 = pd.DataFrame(Mr)
        writer = ExcelWriter('Matriz de massa.xlsx')
        df1.to_excel(writer,'Sheet1', index=False)
        writer.save()
    
#----------------------------------------------------------------------
        self.RestK = Kr
        self.RestM = Mr
#------------------------------------------------------------------------
#4. Calculo das frequencias e modos de vibração
#-----------------------------------------------------------------------------------
    def Eig(self,N,nn,Wall = False,Verbose = True):
        """ 
        Resolve um problema de autovalores e autovetores
        
        :type N: int
        :param N: número de nós desejados
        
        """
        
    ## Kr: Matriz de rigidez restringida
    ## Mr : Matriz de rigidez restringida
    ## N : Número de modos desejados
        
        if Wall == True:
            Kr = self.W
        else:
            Kr = self.RestK
            
        Mr  = self.RestM
        cx  = self.xcoord
        cy  = self.ycoord
        ID  = self.ID
        
        w21,Phi1 = sc.eig(Kr,Mr)

        iw = w21.argsort()
        w21 = w21[iw]          ## Garantindo a ordem dos autovalores e autovetores
        Phi1 = Phi1[:,iw]

        wr = np.real(w21)
        wk = np.sqrt(wr)
        fk = np.real(wk/(2*np.pi))
        wk = 2*np.pi*fk
        
        if Verbose == True:
            for k in range(N):
                print(k+1, "ª frequencia natural = {0:3.2f}Hz".format(fk[k]),"\n")
    
        #fig = plt.figure(1, figsize=(12,8))
        #plt.xlabel('x (m)');plt.ylabel('y (m)')
        #plt.grid(True)
        #plt.axis('off')
    
        #Nb  = len(ID[0,:])
        #ph  =np.zeros(nn)
        #pv  =np.zeros(nn)
    
    
        #CX  = np.zeros(2)
        #CY  = np.zeros(2)
        #for i in range (N):
            #ax  = fig.add_subplot(1,N,i+1)
            #plt.grid(True)
            #ph[4:]+=Phi1[::3,i]
            #pv[4:]+=Phi1[1::3,i]
            #for n in range (Nb):
                #k1      = int(ID[0,n] -1) 
                #k2      = int(ID[1,n] -1)
                #CX[0]   = cx[k1] 
                #CX[1]   = cx[k2]
                #CY[0]   = cy[k1] 
                #CY[1]   = cy[k2] 
                #plt.axis('off')
                #ax.plot(CX,CY,'g:')
                
            
                #CX[0]   = cx[k1] + 10*ph[k1] 
                #CX[1]   = cx[k2] + 10*ph[k2]
                #CY[0]   = cy[k1] + 10*pv[k1]
                #CY[1]   = cy[k2] + 10*pv[k2]
                #plt.axis('off')
                #ax.plot(CX,CY,'b')
                #fig.canvas.draw()
                
#------------------------------------------------------------------------------
        self.Frequency  = fk
        self.RFrequency = wk
        self.Modes      = Phi1
#------------------------------------------------------------------------------
#5. Matriz de amortecimento
#-----------------------------------------------------------------------------------
    
    
    def Rayleigh(self,z1,z2): 
        """
        Define a matriz de amortecimento proporcional as matrizes de 
        massa e rigidez
        
        """
    ## z1 : Fator de amortecimento do primeiro modo 
    ## z2 : Fator de amortecimento do segundo modo
        Kr   = self.RestK
        Mr   = self.RestM
        wk   = self.RFrequency
        zeta =np.zeros(2)
        zeta[0] = z1
        zeta[1] = z2
        a0 = -2*(wk[0]*wk[1])/(wk[0]**2 - wk[1]**2)
        a1 = a0*( zeta[0]*wk[1] - zeta[1]*wk[0])
        a2 = a0*(-zeta[0]/wk[1] + zeta[1]/wk[0])      
    

        Cr  = a1*Mr + a2*Kr              



        df2 = pd.DataFrame(Cr)
        writer = ExcelWriter('Matriz de Amortecimento.xlsx')
        df2.to_excel(writer,'Sheet1', index=False)
        writer.save()
    
#-------------------------------------------------------------------------------
        self.Damping = Cr
#-----------------------------------------------------------------------------------
#6. Espectro de Kanai/Tajimi
#-----------------------------------------------------------------------------------
   
    def SeismicAceleration(self,Ap,tipo,duraçao,dt,w =0,sismoreal = None):
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


        if sismoreal == None:
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
            
            if w!=0:
                 wg = w*2*np.pi
                 zg = 0.3
            
            
            f = np.linspace(0,25,tf)
            df = f[1]-f[0]
            w = 2*np.pi*f
            S0 = (Ap**2)/((pg**2)*(np.pi*wg*((1/(2*zg))+2*zg)))
            Sg = S0*((1+4*(zg**2)*(w/wg)**2)/(((1-(w/wg)**2)**2)+4*(zg**2)*(w/wg)**2))


            plt.figure(2, figsize=(8,4)) 
            plt.plot(f,Sg,'b')
            plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(m²/s³)');
            plt.xlim(0,20); plt.ylim(0,max(Sg)*2);#plt.title(' Espectro de aceleração')
            plt.grid(True)
            
            plt.savefig('Resources/KanaiTajimi.svg')
            np.random.seed(50)
            P = np.random.sample(tf)*2*np.pi
        
   
            ag= np.zeros(tf)
            t = np.linspace(0,duraçao,tf)
            S = np.zeros(tf)



            for i in range(tf):

                S =np.sqrt(2*Sg*df)*np.cos(w*t[i]+ P)
                ag[i] = sum(S)
                
        

            ag*= Ap/np.max(abs(ag))  ## Normalização das acelerações
            
        else:
            
            df = pd.read_fwf(sismoreal,header = None).dropna().values
            ag = []
            print(np.amax(df))
            for i in range(df[:,0].size):
                agl = df[i,:]
                ag = np.hstack((ag,agl))
             
            duraçao = dt* ag.size
            
            t = np.linspace(0,duraçao,ag.size)
            
#---------------------------------------------------------------------------------
        self.time      = t
        self.signal    = ag
        self.Td        = duraçao
        self.Ap        = Ap
#------------------------------------------------------------------------------
#7.Função de envoltória
#-----------------------------------------------------------------------------------
    def Envolve(self):
        """
        Utiliza uma função de envoltória para tornar o sinal mais semelhante
        a um carregamento sismico
        
        """
        ag      = self.signal
        duraçao = self.Td
        t       = self.time
        Ap      = self.Ap
        tf      = len(t)
        env     = np.ones(tf)
        env1    = np.arange(int(0.05*tf))/(0.05*tf)
        env2    = (1.148698355**t[0:int(0.8*tf)])/64/4
        env[0:int(0.05*tf)] = env1
        env[int(0.2*tf):tf] = env2[::-1]

        plt.figure(3,figsize=(8,4))

        plt.plot(t,ag,'b')
        plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (m/s²)');
        plt.xlim(0,duraçao); plt.ylim(-Ap,Ap);plt.title(' Aceleração do solo')
        plt.grid(True)
        plt.savefig('Resources/Sinal_artificial.svg')
        plt.figure(4,figsize=(8,4))
        #plt.plot(t,ag,'c')
        plt.plot(t,env,'r')
        plt.xlabel('Tempo (s)'); plt.ylabel('Envoltória (t)');
        plt.xlim(0,duraçao); #plt.ylim(-Ap,Ap);#plt.title(' Função de envoltória')
        plt.savefig('envoltoria.svg')
        plt.grid(True)

        age = ag*env
        age*= Ap/np.max(abs(age))
        
        plt.figure(5,figsize=(8,4))
        plt.plot(t,age,'k')
        plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (m/s²)');
        plt.xlim(0,0.6*duraçao); plt.ylim(-Ap,Ap);
        #plt.title('Aceleração do solo parametrizada')
        plt.grid(True)
        plt.savefig('Resources/Sinal_envoltoria.svg')
    
        return age

#------------------------------------------------------------------------------
#8. Montagem do vetor de forças
#-----------------------------------------------------------------------------------

    def Sismo(self,age):
        """ 
        Monta um vetor de forças com base nas acelerações geradas pelo espectro
        de Kanai e Tajimi
        
        :type age: list
        :param age: Sinal temporal utilizado para gerar o vetor de forças
        
        """
        Mr = self.RestM
        t  = self.time
    
        n = int (len(Mr[0,:]))
        B = np.zeros((n,1))
        B[::3,0] = np.ones(int(n/3))
        ag = np.zeros((1,len(t)))
        ag[0,:] = age
        F1 = np.dot(B.T,Mr)
        F = np.dot(F1.T,ag) 
    
        #plt.figure(6,figsize=(8,4))   
        #plt.plot(t,F[n-3,:],'b')
        #plt.xlabel('Tempo (s)'); plt.ylabel('Força N')
        #plt.xlim(0,max(t)); plt.ylim(-max(F[n-3])*1.2,max(F[n-3])*1.2);
        #plt.title('Força no último pavimento')
        #plt.grid(True)
        self.Direct = B[:,0]
        return F
#------------------------------------------------------------------------------
# 9. Newmark
#-----------------------------------------------------------------------------------
 # Resolve a equação de equilíbrio por meio do método de Newmark.

    def Newmark(self,F,u0,v0,t,Wall = False):
        """
        Resolve a equação de equilíbrio dinâmico de forma explícita pelo método de 
        Newmark
        
        :type F: list
        :param F: Vetor de forças externas
        
        :type u0: list
        :param u0: Vetor de deslocamentos iniciais
        
        :type v0: list
        :param v0: Vetor de velocidades iniciais
        
        """
        if Wall == True:
            Kr = self.W
        else:
            Kr = self.RestK
            
        Mr = self.RestM
        Cr = self.Damping
        #t  = self.time
#------------------------------------------------------------------------- 
#9.1. Criação dos arrays necessários
#-------------------------------------------------------------------------       
        tf  = int(len(t))
        n   = len(F[:,0])
        A   = np.zeros((n,tf))
        v   = np.zeros((n,tf))
        d   = np.zeros((n,tf))
        dt  = t[1]-t[0]
    
        d[0:n,0] = u0
        v[0:n,0] = v0
    
#------------------------------------------------------------------------- 
#9.2. Determinação das constantes do método de Newmark
#-------------------------------------------------------------------------
        delta = 0.5
        alfa  = 0.25
        a0    = 1/(alfa*(dt**2))
        a1    = 1/(alfa*dt)
        a2    = (1/(2*alfa))-1
        a3    = delta/(dt*alfa)
        a4    = delta/alfa - 1
        a5    = (dt/2)*(delta/alfa - 2)

        A[:,0] = np.dot(np.linalg.inv(Mr),(F[:,0]-np.dot(Cr,
         v[:,0])-np.dot(Kr,d[:,0])))
    
        d4 = a0*Mr + a3*Cr + Kr
        D  = np.linalg.inv(d4)
#------------------------------------------------------------------------- 
#9.3. Resolução da equação de equilíbrio dinâmico
#-------------------------------------------------------------------------
        for i in range(tf-1):
            d1       = np.dot(Mr,(a0*d[:,i]+ a1*v[:,i] + a2*A[:,i]))
            d2       = np.dot(Cr,(a3*d[:,i]+ a4*v[:,i] + a5*A[:,i]))
            d3       = F[:,i+1]+ d1 + d2
            d[:,i+1] = np.dot(D,d3)
            v[:,i+1] = a3*(d[:,i+1] - d[:,i]) - a4*v[:,i] - a5*A[:,i]
            A[:,i+1] = a0*(d[:,i+1] - d[:,i]) - a1*v[:,i] - a2*A[:,i]
    
        return d,v,A
#------------------------------------------------------------------------------
#10. Velocidade média
#-----------------------------------------------------------------------------------
    def Vmedia(z,S1,S3,V0):
        """
        Cria um perfil de velocidades média conforme a altura
        
        :type z: int
        :param z: Altura
        
        :type S1: int
        :param S1: Fator topográfico conforme a NBR6123
        
        :type S3: int
        :param S3: Fator estatístico conforme a NBR6123
        
        :type V0: int
        :param V0: Velocidade básica do vento conforme a NBR6123
        """
        b = .71               # Constante de Karman
        p = .23               # Rugosidade Categoria IV
        Fr = .69              # Altura de referência 
        S2 = b*Fr*(z/10)**p
        Vm = S1*S2*S3*V0
        
        return Vm
#------------------------------------------------------------------------------
#11. Espectro de Davemport
#-----------------------------------------------------------------------------------
  
    def Davenport(self,z,S1,S3,V0,duraçao,dt):
        """
        Usa o Espectro de Davemport para montar um sinal temporal da velocidade 
        flutuante.
        
        :type z: int
        :param z: Altura
        
        :type S1: int
        :param S1: Fator topográfico conforme a NBR6123
        
        :type S3: int
        :param S3: Fator estatístico conforme a NBR6123
        
        :type V0: int
        :param V0: Velocidade básica do vento conforme a NBR6123
        
        :type duraçao:int
        :param duraçao: tempo de duração do sismo
        
        :type dt: int
        :param dt: Passo temporal adotado para a montagem do sinal
        
        """
    
        f1 = 0.0001            # Frequência inicial
        f2 = 3                 # Frequência final
        tf = int(duraçao/dt)
        t = np.linspace(0,duraçao,tf)
        f = np.linspace(f1,f2,tf)           
        z0 = 1                         # Comprimento de Rugosidade (Categoria 4)
        b = .71                        # Constante de Karman
        p = .23                        # Rugosidade Categoria IV
        Fr = .69                       # Altura de referência 
        S2 = b*Fr*(z/10)**p
        V10 = S1*S2*S3*V0
        n = f*1200/V10               
        u = 0.4*V10/np.log(10/z0)
        Sw = 4*n**2*u**2/(f*(1+n**2)**(4/3))
        plt.plot(f,Sw,'b')
        plt.xlim((f1, f2))     
        plt.xscale('log')
        plt.title('Espectro de Davenport')
        plt.xlabel('f (Hz)')
        plt.ylabel('Densidade espectral(m²/s³)')
        plt.grid(True)

        df = f[1]-f[0]
        w  = 2*np.pi*f

        
        P = np.random.sample(tf)*2*np.pi
       
        ag = np.zeros(tf)
        for i in range(tf):
            ag[i] =np.dot(np.sqrt(2*Sw*df),np.cos(w*t[i] + P))
            
        self.time      = t
        self.signal    = ag
        self.frequency = f
        self.Td        = duraçao
#------------------------------------------------------------------------------
# 12. Wind force
#-----------------------------------------------------------------------------------

    def windforce(V,ca,A):
        """
        Retorna o vetor de forças devido ao vento incidindo lateralmente no edifício
        
        :type V: 
        :param V: Vetor de velocidade do vento
        
        :type ca: int
        :param ca: Coeficiente de arrasto do vento
        
        :type A: int
        :param A: Área de influência
        
        """
        F = np.zeros((420,len(V[0,:])))
        q = 0.613*V**2 
        Fp = ca*A*q
        F[::21,:] = Fp
        return F
#------------------------------------------------------------------------------  
#13 . Diferenças finitas centrais
#-----------------------------------------------------------------------------------
    def Finite_diff(self,F,u0,v0):
        """
        Resolve a equação de equilíbrio dinâmico de forma explícita pelo método das
        diferenças finitas centrais
        
        :type F: numpy.ndarray
        :param F: Vetor de forças externas
        
        :type u0: numpy.ndarray
        :param u0: Vetor de deslocamentos iniciais
        
        :type v0: numpy.ndarray
        :param v0: Vetor de velocidades iniciais
        
        """
        Kr = self.RestK
        Mr = self.RestM
        Cr = self.Damping
        t  = self.time
        tf = len(t)

        n      = len(F[:,0])
        d      = np.zeros((n,tf+1))
        d[:,1] = u0
    
        dt     = t[1]-t[0]
    
        a0     =np.dot(np.linalg.inv(Mr),(F[:,0]-np.dot(Cr,v0)-np.dot(Kr,u0)))
        d[:,0] = (dt**2)/2*a0 - dt*v0 + u0
    
    
        C1     = np.linalg.inv(1/(dt**2)*Mr + 1/(2*dt)*Cr)
        C2     = Kr - 2/(dt**2)*Mr 
        C3     = 1/(dt**2)*Mr - 1/(2*dt)*Cr
    
    
        for i in range(1,tf):
     
            d[:,i+1] = np.dot(C1,(F[:,i-1] - np.dot(C2,d[:,i]) - np.dot(C3,d[:,i-1])))
                   
                   
        return d
#------------------------------------------------------------------------------  
# 14. Análise modal
#----------------------------------------------------------------------------------- 
    def Modal_Analysis(self,Fr,u0,v0,n,Wall = False):
        """
         Desacopla as equações de equilíbrio no número de modos definido.
        
        :type Fr: numpy.ndarray
        :param Fr: Vetor de forças externas
        
        :type u0: numpy.ndarray
        :param u0: Vetor de deslocamentos iniciais
        
        :type v0: numpy.ndarray
        :param v0: Vetor de velocidades iniciais
        
        :type n: int
        :param n: Número de modos desejados
        
        """
        if Wall == True:
            Kr = self.W
        else:
            Kr = self.RestK
           
        Kr = self.RestK
        Mr = self.RestM
        Cr = self.Damping
        Phi1 = self.Modes
       
    
        Km = np.diagonal(np.dot(np.dot(Phi1[:,0:n].T,Kr),Phi1[:,0:n]))
        Mm = np.diagonal(np.dot(np.dot(Phi1[:,0:n].T,Mr),Phi1[:,0:n]))
        Cm = np.diagonal(np.dot(np.dot(Phi1[:,0:n].T,Cr),Phi1[:,0:n]))
        Fm = np.dot(Phi1[:,0:n].T,Fr)
        

      
    
        V0 = (np.dot(np.linalg.inv(Phi1),v0)) 
        D0 = (np.dot(np.linalg.inv(Phi1),u0))
#----------------------------------
    
        self.ModalW = np.sqrt(Km/Mm)
        self.ModalM = Mm
        self.Modalz = Cm/(2*Mm*self.ModalW)
        self.ModalF = Fm
        self.V0     = V0
        self.U0     = D0
#------------------------------------------------------------------------------  
# 15. Duhamel
#-----------------------------------------------------------------------------------
    def Duhamel(self):
        """
        Resolve a equação de equilíbrio pelo método de Duhamel.
        Necessita da utlizacão de Modal_Analysis para ser executado para
        modelos com mais de um grau de liberdade.
        
        """
        Mm = self.ModalM
        Wn = self.ModalW
        Zn = self.Modalz
        v0 = self.V0
        u0 = self.U0
        F  = self.ModalF
        t  = self.time
    
        tf = int(len(t))
        n  = len(Wn)
        dt = t[1] - t[0]
        U0 = np.zeros((n,tf))
        U  = np.zeros((n,tf))
    
        for i in range(n):
        
            wn = Wn[i]
            zn = Zn[i]
        
            wd = wn*np.sqrt(1- zn**2)
            e = np.exp(zn*wn*t)
            s = np.sin(wd*t)
            c = np.cos(wd*t)
        
            A = dt*np.cumsum(e*c*F[i,:])
            B = dt*np.cumsum(e*s*F[i,:])
        
            U[i,:] = (u0[i]*c+(v0[i]+u0[i]*zn*wn)*s/wd)/e
        
            U[i,:] +=(A*s-B*c)/e/wd/Mm[i]
        
        return U

#-----------------------------------------------------------------------------------
# 16. Wall
#-----------------------------------------------------------------------------------
    def Wall(self,Nn,andar,Nr,tap):
        
        #"""
        #Determina e inclui na matriz de rigidez da estrutura as barras referentes a 
        #alvenaria conforme suas propriedades e o andar que deve ser aplicado.
        
        #"""
        
#-------------------------------------------------------------------------    
#1.1. Importar planilha do Excel  
#-------------------------------------------------------------------------  
            Arquivo = pd.read_excel('Resources/Wall3.xlsx')
#-------------------------------------------------------------------------     
#1.2. Montar vetor de coordenadas ( X e Y de cada nó) 
#------------------------------------------------------------------------- 
            cx = list(Arquivo['Cx'])[0:Nn]
            cy = list(Arquivo['Cy'])[0:Nn]
#-------------------------------------------------------------------------     
#1.3. Montar matriz identidade( Nó inicial e final de cada barra)
#-------------------------------------------------------------------------  
            Id1     = list(Arquivo['barra (nó 1)'].dropna())
            Id2     = list(Arquivo['barra (nó 2)'].dropna())
            Nb      = len(Id1)
            ID      = np.zeros((2,Nb))
            ID[0,:] = np.array(Id1) + np.ones(Nb)*4*andar
            ID[1,:] = np.array(Id2) + np.ones(Nb)*4*andar
#-------------------------------------------------------------------------     
#1.4. Alocar as propriedades de cada barra
#-------------------------------------------------------------------------  
            #A   = list(Arquivo['Area(m2)'])
            I   = list(Arquivo['Inércia(m4)'])
            RHO = list(Arquivo['Densidade'])
          
#------------------------------------------------------------------------- 
#1.5. Matriz identidade em relação aos graus de liberdade
#-------------------------------------------------------------------------  
            IDG = np.zeros((6,Nb))
    
            for i in range(3):
        
                IDG[i,:]   = ID[0,:]*3-2+i
                IDG[i+3,:] = ID[1,:]*3-2+i
#-------------------------------------------------------------------------        
#1.6. Comprimento de cada barra e cossenos diretores
#-------------------------------------------------------------------------  
            E    = np.ones(Nb)*1.2*10**9
            Ep   = 28*10**9
            Ev   = Ep
            Ea   = 1.2*10**9
            Iv   = 6.75*10**-4
            Ip   = Iv
            Lx   = np.zeros(Nb)
            Ly   = np.zeros(Nb)
            ah   = np.zeros(Nb)
            al   = np.zeros(Nb)
            W    = np.zeros(Nb)
            cosx = np.zeros(Nb)
            cosy = np.zeros(Nb)
            L    = np.zeros(Nb)
        
            for n in range (Nb):
    
                k1      = int(np.round(ID[0,n] -1,0))  # 'k' são os buscadores da matriz Id
                k2      = int(np.round(ID[1,n] -1,0))
         
                Lx[n]   = cx[k2] - cx[k1]
                Ly[n]   = cy[k2] - cy[k1]
                L[n]    = np.sqrt(Lx[n]**2 + Ly[n]**2)
                
                theta   = Ly[n]/Lx[n]
                ah[n]   = np.pi/2*((4*Ep*Ip*Ly[n])/(Ea*tap*np.sin(2*theta)))**(1/4)
                al[n]   = np.pi*((4*Ep*Iv*Lx[n])/(Ea*tap*np.sin(2*theta)))**(1/4)
                W[n]    = (ah[n]**2 + al[n]**2)**0.5
                                 
                cosx[n] = Lx[n]/L[n]
                cosy[n] = Ly[n]/L[n]
                                 
            A   = W*0.15
    
#-------------------------------------------------------------------------        
#1.7. Montagem das matrizes de massa e rigidez
#-------------------------------------------------------------------------  

            K = np.zeros((Nn*3,Nn*3))

            for i in range (Nb):
    
    #1.7.1 Matriz de rigidez local da barra
                Ke =np.array([[E[i]*A[i]/L[i], 0, 0, -E[i]*A[i]/L[i],0 ,0 ],
                  [0, 12*E[i]*I[i]/(L[i]**3), 6*E[i]*I[i]/(L[i]**2), 0,
                   -12*E[i]*I[i]/(L[i]**3),6*E[i]*I[i]/(L[i]**2)],
                  [0,6*E[i]*I[i]/(L[i]**2), 4*E[i]*I[i]/L[i], 0, 
                   -6*E[i]*I[i]/(L[i]**2), 2*E[i]*I[i]/L[i] ],
                  [-E[i]*A[i]/L[i], 0, 0, E[i]*A[i]/L[i],0 ,0 ],
                  [0, -12*E[i]*I[i]/(L[i]**3), -6*E[i]*I[i]/(L[i]**2),
                   0,12*E[i]*I[i]/(L[i]**3),-6*E[i]*I[i]/(L[i]**2)],
                  [0,6*E[i]*I[i]/(L[i]**2), 2*E[i]*I[i]/L[i], 0,
                   -6*E[i]*I[i]/(L[i]**2), 4*E[i]*I[i]/L[i] ]])


    #1.7.3 Matriz de rotação 
                R = np.array([[cosx[i], cosy[i], 0, 0 ,0 ,0],
                  [-cosy[i], cosx[i],0, 0, 0, 0],
                  [0,0,1,0,0,0],                     
                  [0,0,0,cosx[i], cosy[i], 0],
                  [0, 0, 0,-cosy[i], cosx[i],0],
                  [0,0,0,0,0,1]])
        
    #1.7.4 Rotação das matrizes para coordendas globais
                     
                KT = np.dot(np.dot(R.T, Ke),R)             
    
    #1.7.5 Matrizes temporárias
                k_temp1 = np.zeros((Nn*3,Nn*3))    
        
   #1.7.6 Alocação das matrizes temporárias na matriz global
                m = int(IDG[0,i]-1)
                n = int(IDG[2,i])
                o = int(IDG[3,i]-1)
                p = int(IDG[5,i])
    
                k_temp1[m:n,m:n] = KT[0:3,0:3]
                k_temp1[o:p,m:n] = KT[3:6,0:3]
                k_temp1[m:n,o:p] = KT[0:3,3:6]
                k_temp1[o:p,o:p] = KT[3:6,3:6]
    
                K += k_temp1 
    
            Kr_1 = np.delete(K,Nr,0)
            Kr  = np.delete(Kr_1,Nr,1)

#-------------------------------------------------------------------------
            self.W = Kr + self.RestK
    
#-------------------------------------------------------------------------
#17. Story drift
#-------------------------------------------------------------------------
    def Storydrift(d):
        hd    = d[::12,:]
        sd    = np.zeros(hd.shape)
        sdmax = np.zeros(len(hd[:,0])+1)
        
        sd[0,:] = hd[0,:]
        for i in range(1,len(hd[:,0])):
            
            sd[i,:] = hd[i,:] - hd[i-1,:]
        
        for i in  range(1,len(hd[:,0])+1):
            
            sdmax[i] = np.amax(sd[i-1,:])
        
        return sdmax
#-------------------------------------------------------------------------
#17. Story drift2
#-------------------------------------------------------------------------
    def Storydrift2(umax):
        sdmax = np.zeros(umax.shape)
        sdmax[0] = umax[0]
        sdmax[1:]= umax[1:]-umax[0:umax.size-1]
        
        return sdmax
#-------------------------------------------------------------------------
#17. Freqdomain
#-------------------------------------------------------------------------            
    def FreqDomain(self,Sa,T,Verbose = True):
        """ Calcula do vetor Força máxima para determinado espectro de Pseu
        doaceleração
        """
        w = self.RFrequency
        Mt = 8000
        Wi = 0
        i = 0
        
        while Wi<= 0.9:
                Ln = -self.Modes[:,i].T@self.RestM@self.Direct
                Mn = self.Modes[:,i].T@self.RestM@self.Modes[:,i]
                Wi+= (Ln**2/Mn)/Mt
                i+=1
        if Verbose ==True:
            print('Número de modos necessários =',i+1,'Wi = ',Wi*100,'%')
        F = np.zeros([len(self.Modes[:,0]),i+1])

        for n in range (i+1): 
            Ln     = -self.Modes[:,n]@self.RestM@self.Direct
            Mn     = self.Modes[:,n].T@self.RestM@self.Modes[:,n]
            Sai    = np.interp(2*np.pi/w[n],T,Sa)
            F[:,n] = self.RestM@(Ln/Mn*Sai*self.Modes[:,n])

        return F
#-------------------------------------------------------------------------    
#19. SQRSS
#------------------------------------------------------------------------- 
    def SQRSS(u):
        umax = np.sqrt(np.sum(u**2),1)

        return umax
#------------------------------------------------------------------------- 
# 20. CQC
#------------------------------------------------------------------------- 
    def CQC(self,u,zt):
        umax = np.zeros(u.size)
        n = len(u[0,:])
        w = self.RFrequency[0:n]
        for i in range (n):
            for j in range (n):
                bij = w[i]/w[j]
                pij = (8*zt**2*(1+bij)*bij**(1.5))/((1-bij**2)**2+4*zt**2*bij
                *(1+bij**2)**2)
                umax+= pij*u[:,i]*u[:,j]

        return umax            
#-------------------------------------------------------------------------
#2. Newmark Não-linear
#-------------------------------------------------------------------------
def NLNewmarkM(k,m,c,p,u0,v0,t,e,Wall,Port,E,model):
    
# 0.0 iniciando os parâmetros

    nt  = int(len(t))
    n   = int(len(p[:,0]))
    a   = np.zeros((n,nt))
    v   = np.zeros((n,nt))
    u   = np.zeros((n,nt))
    fs  = np.zeros((n,nt))
    kt  = np.zeros((nt,n,n))
    kht = np.zeros((nt,n,n))
    ph  = np.zeros((n,nt))
    R   = np.zeros((n,nt))
    
#1.0 Determinação dos parâmetros iniciais:
    b = 0.25
    y = 0.50

#1.1 passo temporal

    dt = t[1] - t[0]
    
#1.2 variáveis de estado    
    kt[0] = k
    fs[:,0] = np.dot(kt[0],u[:,0]) 
    
#1.3 cálculo da aceleração inicial:

    a[:,0] = np.linalg.inv(m)@(p[:,0] - c@v[:,0] - fs[:,0])
            
# 1.4 Determinação das consntantes:
            
    a1 = 1/(b*dt**2)*m + y/(b*dt)*c
            
    a2 = 1/(b*dt)*m + (y/b -1)*c
            
    a3 = (1/(2*b)-1)*m + dt*(y/(2*b)-1)*c
    
#2.0 início da iteração temporal
    for i in range(nt-1):
            print('Passo temporal:',i)
# 2.1 iniciação dos valores
            
            u[:,i+1]  = u[:,i]
            fs[:,i+1] = fs[:,i]
            kt[i+1]   = kt[i]

# 2.2 Cálculo da força equivalente 

            ph[:,i+1] = p[:,i+1] + a1@u[:,i] + a2@v[:,i] + a3@a[:,i]
            
# 3.0 início da iteração de Newton-Raphson
            j = 0

            while True:

                R[:,i+1] = ph[:,i+1] - fs[:,i+1] - a1@u[:,i+1]
                
                j+=1
                
                if j>50:
                    
                    print('erro: não converge na iteração',i)
               
                    
                if np.linalg.norm(R[:,i+1]) <= e:
                    
                    print('Numero de iterações:',j)
                    break
                    
                else:
                    
                    kht[i+1] = kt[i+1] + a1
                
                    du = np.linalg.inv(kht[i+1])@R[:,i+1]
                    
                    u[:,i+1] += du
                    
                    kt[i+1] = model(u[:,i+1],Wall,Port,E)
                    
                    fs[:,i+1] = np.dot(kt[i+1],u[:,i+1])

# 4.0 Determinação da velocidade e da aceleração

            v[:,i+1] = y/(b*dt)*(u[:,i+1] - u[:,i]) + (1-y/b)*v[:,i] + dt*(1-y/(2*b))*a[:,i]
            
            a[:,i+1] = 1/(b*dt**2)*(u[:,i+1] - u[:,i]) - 1/(b*dt)*v[:,i] - (1/(2*b) - 1)*a[:,i]
            
    return u,v,a,kt,fs,R

