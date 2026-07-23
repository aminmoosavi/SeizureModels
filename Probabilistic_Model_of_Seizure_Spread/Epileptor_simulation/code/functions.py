import numpy as np
from numba import jit
import time


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
        "IL": np.array([60, 63]),
        "fb": np.array([47, 59, 80]),
        "id039": np.array([82, 86, 115, 125, 134, 153, 155]),
        "id001": np.array([134, 153, 155]),
        "id017": np.array([34, 50, 52, 53, 72, 73, 74]),
        "random": np.array([0])
    }
    return patients[patient_name]


def x0_w_init(patient_name,EZ_ind):
    patients = {
        "IL": np.array([[-2.311, -2.311],[0.0,0.0]]),
        "fb": np.array([[-2.311, -2.311,-2.311],[0.0,0.0,0.0]]),
        "id039": np.array([[-2.311, -2.311, -2.311, -2.311, -2.311, -2.311,-2.311],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]),
        "id001": np.array([[-2.311, -2.311,-2.311],[0.0,0.0,0.0]]),
        "id017": np.array([[-2.311, -2.311, -2.311, -2.311, -2.311, -2.311, -2.311],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]),
        "random": np.array([[-2.311],[0.0]])
    }    

    return patients[patient_name][:,EZ_ind]    




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
