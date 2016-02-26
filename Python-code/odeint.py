from pylab import *
from scipy.integrate import odeint
import numpy as np
import time
timestr = time.strftime("%m%d%Y")
    
def solution(y20):  #solve for odes
    error=1e-4
    y10=y30=y40=1e-12
    Sigma=-1.1*pow(P,2.0/3.0)
    yinit = array([y10,y20,y30,y40]) # initial values
    S1,S2,Sinc=0.0,20.0,5.0
    while S2<=20.0:
        S = np.arange(S1,S2,0.001) # space interval 
        S1=S2
        S2=S2+Sinc
        y = odeint(deriv,yinit,S) 
        yinit=y[-1,:]       
        region=list(range(0,len(y)-1))
        for i in region:  
            if (y[i,0]<=pi+error and y[i,0]>=pi-error):
                data1=[y20,y[i,0],y[i,1],y[i,2],y[i,3],S[i],Sigma]
                print(data1)
                return data1
            
def deriv(y,S): # Coupled odes

    C0=0.0
    Sigma=-1.1*pow(P,2.0/3.0)
      
    dtheta=y[1] # equation 3.5a
    dU=-y[1]*cos(y[0])/y[3]+cos(y[0])*sin(y[0])/(y[3]*y[3])+\
        y[2]*sin(y[0])/y[3]+P*y[3]*cos(y[0])/2 # equation 3.5b
    dgamma=(y[1]-C0)*(y[1]-C0)/2-sin(y[0])*sin(y[0])/(2*y[3]*y[3])+\
            P*y[3]*sin(y[0])+Sigma # equation 3.5c
    dX=cos(y[0]) # equation 3.5d
    return [dtheta,dU,dgamma,dX]
 
def ChangeP(VarP):
    
    global P,scale
    P=VarP # Model kinetic parameters
    scale=pow(P,-1.0/3.0)
    y20range=np.arange(-2.0/scale,2.0/scale,0.001)
    total = []
    for j in y20range:
        mm=solution(j)
        if mm!=None:
            total.append(mm) 
    final=list(total)
    return final

def batchsave():
    Prange=np.arange(1,2,1)
    for kk in Prange:
        start_time = time.time()
        with open('%s_P=%-.2f.txt'%(timestr,kk), "wb") as outfile:
            final=ChangeP(kk)
            np.savetxt(outfile, final, fmt='%-7.8f')
        print("--- %s seconds ---" % (time.time() - start_time)) 
    
batchsave()