import numpy as np
from numba import jit
import time
from scipy import linalg


@jit(nopython=True)
def f1(x1, x2, z):
    # print(x1)
    N = int(x1.shape[0])
    f = np.zeros(N)
    for i in range(N):
        if x1[i] < 0:
            f[i] = (x1[i])**3.0 - 3.0 * (x1[i]**2)
        elif x1[i] >= 0:
            f[i] = x1[i] * (x2[i] - 0.6 * (z[i] - 4.0)**2)
    # print(f)
    return f


@jit(nopython=True)
def f2(x2):
    N = int(x2.shape[0])
    f = np.zeros(N)
    for i in range(N):
        if x2[i] < -0.25:
            f[i] = 0
        elif x2[i] >= -0.25:
            f[i] = 6 * (x2[i] + 0.25)
    return f


@jit(nopython=True)
def fz(z):
    N = z.shape[0]
    f = np.zeros(N)
    for i in range(N):
        if z[i] > 0:
            f[i] = z[i]
        elif z[i] <= 0:
            f[i] = z[i]+0.1 * z[i]**7
    # print(f)
    return f


@jit(nopython=True)
def History(x1, k1 , noise_x1, a0 , Tau, W, t0,h, t1):
    N = x1.shape[1]
    hist = np.zeros(N)
    
    for i in range(N):
        for j in range(N):
            m = Tau[i, j]
            tt1=t1-m
            if(tt1<0):
                tt1=t0+tt1
            a = x1[tt1, j] + a0*h*k1[tt1, j] + \
                a0*np.sqrt(h)*noise_x1[tt1, j]
            b = x1[t1, i] + a0*h*k1[t1, i] + a0*np.sqrt(h)*noise_x1[t1, i]
            hist[i] = hist[i] + W[i, j] * (a - b)

    return hist




def EZs(patient_name):
    patients = {
        "P1": np.array([60, 63]),
        "P2": np.array([47, 59, 80]),
        "P5": np.array([82, 86, 115, 125, 134, 153, 155]),
        "P3": np.array([134, 153, 155]),
        "P4": np.array([34, 50, 52, 53, 72, 73, 74]),
        "random": np.array([0])
    }
    return patients[patient_name]


def x0_w_init(patient_name,EZ_ind):
    patients = {
        "P1": np.array([[-2.33, -2.33],[0.0,0.0]]),
        "P2": np.array([[-2.33, -2.33,-2.33],[0.0,0.0,0.0]]),
        "P5": np.array([[-2.20, -2.20, -2.20, -2.20, -2.33, -2.33,-2.33],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]),
        "P3": np.array([[-2.20, -2.20,-2.20],[0.0,0.0,0.0]]),
        "P4": np.array([[-2.33, -2.33, -2.33, -2.33, -2.33, -2.33, -2.33],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]),
        "random": np.array([[-2.33],[0.0]])
    }    

    return patients[patient_name][:,EZ_ind]    



def ddff(x1,x2,z):

    if x1<0:
        f1=x1**3 -3*x1**2
        df1dx1= 3 * x1**2 - 6 * x1
        df1dx2=0
        df1dz=0
    else:
        f1=x1 * (x2 - 0.6 * (z-4)**2)
        df1dx1= x2 - 0.6 * (z-4)**2
        df1dx2=x1
        df1dz= -x1*1.2*(z-4)
 
    if x2>=-0.25:
        df2dx2=6
        f2=6*(x2+0.25)
    else:
        f2=0
        df2dx2=0
    
    return f1,f2,df1dx1,df1dx2,df2dx2,df1dz


def Jacobian6D(X,W,model):
    M=int(X.size)
    N=int(M/6)

    I1=model["I1"]
    I2=model["I2"]
    r=model["r"]
    gamma=model["gamma"]
    tau2=model["tau2"]
    x0=model["x0"]

    x1=X[ 0*N:1*N ]
    y1=X[ 1*N:2*N ]
    z= X[ 2*N:3*N ]
    x2=X[ 3*N:4*N ]
    y2=X[ 4*N:5*N ]
    g= X[ 5*N:6*N ]

    dx1=np.zeros(N)
    dy1=np.zeros(N)
    dz=np.zeros(N)
    dx2=np.zeros(N)
    dy2=np.zeros(N)
    dg=np.zeros(N)

    II=np.eye(N)

    x1x1=np.zeros([N,N]); y1x1=np.zeros([N,N]); zx1=-(r) * W; x2x1=np.zeros([N,N])
    x1y1=II;         y1y1=-II;        zy1=np.zeros([N,N]);    x2y1=np.zeros([N,N])
    x1z=np.zeros([N,N]);  y1z=np.zeros([N,N]);  zz=-(r)*II;   x2z=-0.3*II
    x1x2=np.zeros([N,N]); y1x2=np.zeros([N,N]); zx2=np.zeros([N,N]);    x2x2=np.zeros([N,N])
    x1y2=np.zeros([N,N]); y1y2=np.zeros([N,N]); zy2=np.zeros([N,N]);    x2y2=-II
    x1g=np.zeros([N,N]);  y1g=np.zeros([N,N]);  zg=np.zeros([N,N]);     x2g=0.002*II

    y2x1=np.zeros([N,N]); gx1=II
    y2y1=np.zeros([N,N]); gy1=np.zeros([N,N])
    y2z=np.zeros([N,N]);  gz=np.zeros([N,N])
    y2x2=np.zeros([N,N]); gx2=np.zeros([N,N])
    y2y2=(-1/tau2)*II; gy2=np.zeros([N,N])
    y2g=np.zeros([N,N]);  gg=-gamma*II

    for i in range(N):
   
        f1,f2,df1dx1,df1dx2,df2dx2,df1dz = ddff(x1[i],x2[i],z[i])
        
        x1z[i,i]=-1-df1dz
        zx1[i,i]=(r) * ( 4 + np.sum(W[:,i]) )
        x1x1[i,i]=-df1dx1
        x1x2[i,i]=df1dx2
        y2x2[i,i]=(1/tau2) * df2dx2
        y1x1[i,i]=-10*x1[i]
        x2x2[i,i]= 1-3*x2[i]**2
    
        
        
        dx1[i]=y1[i]-f1-z[i]+I1
        dy1[i]=1 - 5*x1[i]**2 -y1[i]
        dz[i]=(r) * ( 4*(x1[i]-x0[i]) - z[i] - np.dot(x1,W[:,i]) + np.sum( W[:,i] ) * x1[i] )
        dx2[i]=-y2[i] + x2[i] - x2[i]**3 +I2 +0.002*g[i] - 0.3*(z[i]-3.5)
        dy2[i]=(1/tau2)*(-y2[i]+f2)
        dg[i]= x1[i] - gamma*g[i]  
    

    J=np.concatenate((np.concatenate((x1x1 , x1y1 , x1z , x1x2 , x1y2 , x1g ),axis=1),\
    np.concatenate((y1x1 , y1y1 , y1z , y1x2 , y1y2 , y1g ),axis=1) ,\
    np.concatenate((zx1  , zy1  , zz  , zx2  , zy2  , zg),axis=1) ,\
    np.concatenate((x2x1 , x2y1 , x2z , x2x2 , x2y2 , x2g),axis=1) ,\
    np.concatenate((y2x1 , y2y1 , y2z , y2x2 , y2y2 , y2g),axis=1) ,\
    np.concatenate((gx1  , gy1  , gz  , gx2  , gy2  , gg),axis=1) ),axis=0)


    Xdot=np.concatenate((dx1,dy1,dz,dx2,dy2,dg))
    return J




def Xprime(X,W,model):
    M=int(X.size)
    N=int(M/6)

    I1=model["I1"]
    I2=model["I2"]
    r=model["r"]
    gamma=model["gamma"]
    tau2=model["tau2"]
    x0=model["x0"]

    x1=X[ 0*N:1*N ]
    y1=X[ 1*N:2*N ]
    z= X[ 2*N:3*N ]
    x2=X[ 3*N:4*N ]
    y2=X[ 4*N:5*N ]
    g= X[ 5*N:6*N ]

    dx1=np.zeros(N)
    dy1=np.zeros(N)
    dz=np.zeros(N)
    dx2=np.zeros(N)
    dy2=np.zeros(N)
    dg=np.zeros(N)

    II=np.eye(N)



    for i in range(N):
   
        f1,f2,df1dx1,df1dx2,df2dx2,df1dz = ddff(x1[i],x2[i],z[i])
    
        dx1[i]=y1[i]-f1-z[i]+I1
        dy1[i]=1 - 5*x1[i]**2 -y1[i]
        dz[i]=(r) * ( 4*(x1[i]-x0[i]) - z[i] - np.dot(x1,W[:,i]) + np.sum( W[:,i] ) * x1[i] )
        dx2[i]=-y2[i] + x2[i] - x2[i]**3 +I2 +0.002*g[i] - 0.3*(z[i]-3.5)
        dy2[i]=(1/tau2)*(-y2[i]+f2)
        dg[i]= x1[i] - gamma*g[i]  
    


    Xdot=np.concatenate((dx1,dy1,dz,dx2,dy2,dg))
    return Xdot










# def Random_connectivity(N,N0,p):
#     out_dir="./files/"
#     a=1
#     b=5
#     Nez=int(N/N0)
#     seed=int(1000*time.time())
#     seed = np.mod(seed,2**31)+1
#     np.random.seed(seed)
#     W=np.zeros([N,N])
#     Tau=np.zeros([N,N])
#     A=np.zeros([N,N])
#     for i in range(N-1):
#         if (i==0):
#             nn=Nez
#         else:
#             nn=1
#         for k in range(nn):
#             for j in range(i+1,N):
#                 r=np.random.random(1)
#                 if(r<p):
#                     A[i,j]=1
#                     A[j,i]=1
#                     W[i,j]=W[i,j]+np.random.beta(a,b,1)
#                     W[j,i]=W[i,j]
#                     Tau[i,j]=np.random.random(1)*200
#                     Tau[j,i]=Tau[i,j]
    
#     W=W*(N0/N)
#     return W , Tau
    # fw=open(out_dir+'weights.txt','w')
    # fTau=open(out_dir+'tract_lengths.txt','w')
    # fA=open(out_dir+'adj.txt','w')
    # np.savetxt(fw,W,fmt='%1.6f')
    # np.savetxt(fTau,Tau,fmt='%1.6f')
    # np.savetxt(fA,A,fmt='%1.6f')
    # zipfile.ZipFile(out_dir+'connectivity_random_N_'+str(N)+'_p_'+str(int(100*p) )+'.zip',mode='w').write(out_dir+'weights.txt','weights.txt')
    # zipfile.ZipFile(out_dir+'connectivity_random_N_'+str(N)+'_p_'+str(int(100*p))+'.zip',mode='a').write(out_dir+'tract_lengths.txt','tract_lengths.txt')
    # zipfile.ZipFile(out_dir+'connectivity_random_N_'+str(N)+'_p_'+str(int(100*p))+'.zip',mode='a').write(out_dir+'adj.txt','adj.txt')
