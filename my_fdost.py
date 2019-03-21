#!/usr/bin/env python
#
# 1-D fast discrete orthonormal stockwell transform code 
#   for convenience , signal length should be pow of 2
# parameter:
#    data——numpy.ndarray , the origion 1-D signal
#
# Date:2019/3/20
#
# Copyright (c) by GaoSong


import numpy as np
import math
from scipy import signal


def fdost(data,origion=False):
    # return  DOST  coefficient , first positive frequency after negative(same as fft)
    if not origion: 
        N=len(data)
        fdata=np.fft.fft(data)
        p_limit=math.floor(math.log2(N))
        res=[]
        ns=0
        for i in range(p_limit):
            mu,beta,tau=vbt(i,N)
            r_matrix=np.array([])
            for j in range(beta):
                r_matrix=np.hstack((r_matrix,np.array((-1)**j)))
            
            v_matrix=r_matrix[:,np.newaxis].dot(np.fft.ifft(fdata[ns:ns+beta]).reshape(1,beta))*np.sqrt(beta)
            ns+=beta
            res.append(v_matrix)
        for i in range(-p_limit,0):
            mu,beta,tau=vbt(i,N)
            r_matrix=np.array([])
            for j in range(beta):
                r_matrix=np.hstack((r_matrix,np.array((-1)**j)))
            
            v_matrix=r_matrix[:,np.newaxis].dot(np.fft.ifft(fdata[ns:ns+beta]).reshape(1,beta))*np.sqrt(beta)
            ns+=beta
            res.append(v_matrix)
        return res

    # return total redundant array in positive frequency
    elif origion:
        N=len(data)
        fdata=np.fft.fft(data)
        p_limit=math.floor(math.log2(N))
        ns=0
        for i in range(p_limit):
            mu,beta,tau=vbt(i,N)
            r_matrix=np.array([])
            for j in range(beta):
                r_matrix=np.hstack((r_matrix,np.array((-1)**j)))
            v_matrix=r_matrix[:,np.newaxis].dot(np.fft.ifft(fdata[ns:ns+beta]).reshape(1,beta))*np.sqrt(beta)
            ns+=beta
            if i == 0 :
                fdost_res=np.tile(v_matrix,N)
            elif i == 1 :
                fdost_res=np.vstack((fdost_res,np.tile(v_matrix,N)))
            else :
                for m in range(beta):
                    tmp=np.array([])
                    for n in range(beta):
                        tmp=np.hstack((tmp,np.tile(np.asarray(v_matrix[m][n]),int(N/beta))))
                        
                    fdost_res=np.vstack((fdost_res,tmp))
                
        return fdost_res

def vbt(p,N):
    # compute orthonormal basis vector parrameter mu,beta,tau
    if p == 0:
        mu=0;beta=1;tau=0
    elif p == 1:
        mu=1;beta=1;tau=0
    elif p > 1:
        mu=2**(p-1)+2**(p-2)
        beta=2**(p-1)
        tau=list(range(beta-1))

    elif p == -1:
        mu=-1;beta=1;tau=0
    elif p == -math.floor(math.log2(N)):
        mu=-2**(-p-1);beta=1;tau=0
    else:
        mu=-2**(-p-1)-2**(-p-2)+1
        beta=2**(-p-1)
        tau=list(range(beta-1))
    return mu,beta,tau


def idost(fdost_res):
    N=2**(len(fdost_res)/2)
    p_limit=math.floor(math.log2(N))
    ori=np.array([])
    for i in range(p_limit):
        mu,beta,tau=vbt(i,N)
        r_matrix=np.array([])
        for j in range(beta):
            r_matrix=np.hstack((r_matrix,np.array((-1)**j)))
        iff=np.linalg.pinv(r_matrix[:,np.newaxis]).dot(fdost_res[i]/np.sqrt(beta)).flatten()
        ff=np.fft.fft(iff)
        ori=np.hstack((ori,ff))
    for i in range(-p_limit,0):
        mu,beta,tau=vbt(i,N)
        r_matrix=np.array([])
        for j in range(beta):
            r_matrix=np.hstack((r_matrix,np.array((-1)**j)))
        iff=np.linalg.pinv(r_matrix[:,np.newaxis]).dot(fdost_res[i+len(fdost_res)]/np.sqrt(beta)).flatten()
        ff=np.fft.fft(iff)
        ori=np.hstack((ori,ff))
    ori=np.fft.ifft(ori).real
    return ori




# if __name__=='__main__':

#     t=np.linspace(0,10,4096)
#     w=signal.chirp(t,f0=12.5,f1=2.5,t1=10,method='linear')

#     dostres=fdost(w,origion=True) 
#     dostres=dostres[0:256]
#     extent=(t[0],t[-1],0,25.6)

#     import matplotlib.pyplot as plt
    
#     fig,ax=plt.subplots(2,1)
#     ax[0].plot(t,w)
#     ax[1].imshow(np.abs(dostres),origin='lower',extent=extent)
#     ax[1].axis('tight')
#     ax[1].set(xlabel='time',ylabel='frequency')
#     plt.show()