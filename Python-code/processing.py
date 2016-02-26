from scipy.integrate import odeint
import numpy as np
import math
import pylab as plt
import time
timestr = time.strftime("%m%d%Y")

 # Coupled odes
def deriv(y,S):
    
    C0=0.0
    Sigma=-1.1*pow(P,2.0/3.0) 
     
    dtheta=y[1] # equation 3.5a
    dU=-y[1]*math.cos(y[0])/y[3]+math.cos(y[0])*math.sin(y[0])/(y[3]*y[3])+\
        y[2]*math.sin(y[0])/y[3]+P*y[3]*math.cos(y[0])/2 # equation 3.5b
    dgamma=(y[1]-C0)*(y[1]-C0)/2-math.sin(y[0])*math.sin(y[0])/(2*y[3]*y[3])+\
            P*y[3]*math.sin(y[0])+Sigma # equation 3.5c
    dX=math.cos(y[0]) # equation 3.5d
    return [dtheta,dU,dgamma,dX]     

# Change P values and read data
Prange=np.arange(1,2,1)
for kk in Prange:
    total = np.loadtxt('02242016_P=%-.2f.txt'%kk)
# define paramaters
    P=kk;xi0=0;x0=0.0001;ganma0=0;deltas=0.001;
    scale=pow(P,-1.0/3.0);scale1=pow(P,1.0/3.0);Sigma=-1.1*pow(P,2.0/3.0);    
# Plot fig 16
    fig=plt.figure()
    plt.rc('font', size=6)
    plt.subplot(2,1,1)
    plt.hold(True)
    u0x=[];x0x=[];
    u0x=total[:,0]
    x0x=total[:,4]
    plt.plot(np.array(u0x)*scale,np.array(x0x)*scale1,'.r')     
# find minimal values
    data=np.array(total) 
    [row1,col]=data.shape
    solutions=[]
    for i in np.arange(1,row1-1,1):
        if True:
#        if data[i,6]==Sigma:
            if data[i,4]*scale1<=0.15 and (data[i,4]<=data[i-1,4] and data[i,4]<=data[i+1,4]):
                plt.hold(True)
                u0y=total[i,0]
                x0y=total[i,4]
                plt.plot(np.array(u0y)*scale,np.array(x0y)*scale1,'.b')
                printa=[i,data[i,:]]
                solutions.append(printa)
#                print(printa)
                 
#######################################
####### start of bulk analysis  #######
#######################################                
    solength=len(solutions)
    for l in np.arange(0,solength,1):  
        sl=solutions[l][1][5]
#        print ('sl=%f'%sl)
        u0=solutions[l][1][0]
#        print ('u0=%f'%u0)    
        yinit=[1e-12,u0,1e-12,1e-12]#    
        S = np.arange(0.,sl,deltas) # space interval 
        y = odeint(deriv,yinit,S) 
        a1=[row[0] for row in y];b1=[row[1] for row in y];
        c1=[row[2] for row in y];d1=[row[3] for row in y];
        s1=S; 
        length2=len(s1)
        x=np.zeros(length2)
        z=np.zeros(length2)
        c2=np.zeros(length2)
        for n in np.arange(0,length2-1,1):
            x[n+1]=x[n]+deltas*math.cos(a1[n]);
            z[n+1]=z[n]-deltas*math.sin(a1[n]); 
            c2[n+1]=c2[n]+deltas*(math.sin(a1[n])/d1[n]);    
        plt.subplot(4,solength,3*solength+1+l)          
        plt.gca().set_aspect('equal', adjustable='box')
        plt.plot(x,z,-x,z)
        plt.title('sl=%.3f,u0=%.3f'%(sl,u0))
    plt.savefig('02242016_P=%-.2f.png'%kk,dpi =300)
    plt.close(fig)    
#######################################    
####### end of bulk analysis  #########
#######################################
    
                   
#########################################
####### start of detail analysis  #######
#########################################
#    solength=len(solutions)
#    for l in np.arange(0,solength,1):  
#        sl=solutions[l][1][5]
##        print ('sl=%f'%sl)
#        u0=solutions[l][1][0]
##        print ('u0=%f'%u0)    
#        yinit=[1e-12,u0,1e-12,1e-12]#    
#        S = np.arange(0.,sl,deltas) # space interval 
#        y = odeint(deriv,yinit,S) 
#        a1=[row[0] for row in y];b1=[row[1] for row in y];
#        c1=[row[2] for row in y];d1=[row[3] for row in y];
#        s1=S; 
#        length2=len(s1)
#        x=np.zeros(length2);z=np.zeros(length2);c2=np.zeros(length2)
#        for n in np.arange(0,length2-1,1):
#            x[n+1]=x[n]+deltas*math.cos(a1[n]);
#            z[n+1]=z[n]-deltas*math.sin(a1[n]); 
#            c2[n+1]=(math.sin(a1[n])/d1[n]);  
#        plt.figure()
#        plt.subplot(2,3,1)    
#        plt.plot(x,z,-x,z)
#        plt.title('sl=%f,u0=%f'%(sl,u0))
#        plt.gca().set_aspect('equal', adjustable='box')
#        plt.subplot(2,3,2)
#        plt.plot(s1,x)
#        plt.xlabel('S');plt.ylabel('X')       
#        plt.subplot(2,3,3)
#        plt.plot(s1,z)
#        plt.xlabel('S');plt.ylabel('z')
#        plt.subplot(2,3,4)
#        plt.plot(s1,b1)
#        plt.xlabel('S'); plt.ylabel('c1')
#        plt.subplot(2,3,5)
#        plt.plot(s1,c2)
#        plt.xlabel('S');plt.ylabel('c2')
#        plt.subplot(2,3,6)
#        plt.plot(s1,a1)
#        plt.xlabel('S');plt.ylabel('PSI')
#        plt.savefig('02172016_P=1.00.png'%kk,dpi =300)
#        plt.close() 
#######################################        
####### End of detail analysis  #######
####################################### 