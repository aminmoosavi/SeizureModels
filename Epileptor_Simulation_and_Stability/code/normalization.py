import numpy as np


def normalize(W,method):
    N = W.shape[0]
    for i in range(N):
        W[i,i]=0.0
    clipp(W, 0.95)

    if method=="global":
        W=W/np.max(W)

    elif method=="logarithmic":
        W=np.log(W+1)
        W=W/np.max(W)


    elif method=="square_root":
        W=np.sqrt(W)
        W = W / np.max(W)

    elif method=="local":
        W=local_norm(W)

    return W



def clipp(W,p):
    N=W.shape[0]
    M=int((N*(N-1))/2)
    U=np.zeros(M)
    k=0
    for i in range(N-1):
        for j in range(i+1,N):
            U[k]=W[i,j]
            k=k+1
    
    UU=np.sort(U)
    Mp=int(M*p)
    Up=UU[Mp]
    # print(Up)
    for i in range(N-1):
        for j in range(i+1 , N):
            
            if W[i,j]>Up:
                W[i,j]=Up
                W[j,i]=Up
    return W
    




def local_norm(W):
    N=W.shape[0]
    for i in range(N):
        W[i,:]=W[i,:]/np.sum(W[i,:])
    return W






# def clippp(W,p):
#     N=W.shape[0]
#     M=int((N*(N-1))/2)
#     U=np.zeros(M)
#     k=0
#     for i in range(N-1):
#         for j in range(i+1,N):
#             U[k]=W[i,j]
#             k=k+1
    
    
#     Up=np.percentile(U,int(100*p))
#     print(Up)
#     for i in range(N-1):
#         for j in range(i+1 , N):
#             # if i==j:
#             #     W[i,j]=0

#             if W[i,j]>Up:
#                 W[i,j]=Up
#                 W[j,i]=Up
#     return W
