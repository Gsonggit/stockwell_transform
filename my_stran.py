#!/usr/bin/env python
#
# 1-D redundant stockwell transform code modified by GaoSong from origion author Stockwell matlab code(st.m also included in module)
#
# parameter:
#    data——numpy.ndarray , the origion 1-D signal
#    minfreq——int  ,  the minimum frequence in result
#    maxfreq——int  ,  the maximum frequence in result
#    samprate——int  ,  the time interval of origion signal
#    freqsamprate——int   ,  the frequence interval in result 
#    remove_edge——boolean   ,  remove a least-squares fit parabola and put a 5% hanning taper on origion signal
#    analytic_signal——boolean   ,   turn the real_value signal into analytic signal and then start s-transform
#    factor——int   ,   the width factor of localizing gaussian window
# 
#
#Date:2019/1/1
#
#Copyright (c) by GaoSong

import numpy as np
from obspy.core import Stream,Trace
from scipy import signal 

def st(data,minfreq,maxfreq,samprate=1,freqsamprate=1,remove_edge=False,analytic_signal=False,factor=1):
    if data.shape[0] <= 1 or len(data.shape)>1 :
        raise TypeError('input data invalid ,please check!') 	
    orig=data
    st_res=np.zeros((int((maxfreq-minfreq)/freqsamprate)+1,len(data)),dtype='c8')	
    #data.reshape(-1,1)
    if remove_edge:
        print('remove_edge selected;  Remove trend with polynomial fit and taper!')
        tmp=Trace(data=orig)
        tmp.detrend('polynomial')
        tmp.taper(0.05)
        orig=tmp.data
    if analytic_signal:
        print('analytic_signal selected;  Calculating analytic signal!')
        orig=signal.hilbert(orig)
    vec=np.hstack((np.fft.fft(orig),np.fft.fft(orig)))
    
    if minfreq == 0:
       st_res[0]=np.mean(orig)*np.ones(len(data))
    else:
       st_res[0]=np.fft.ifft(vec[minfreq+1:minfreq+len(data)+1]*g_window(len(data),minfreq,factor))
       
    for i in range(freqsamprate,(maxfreq-minfreq)+1,freqsamprate):
        st_res[int(i/freqsamprate)]=np.fft.ifft(vec[minfreq+i+1:minfreq+i+len(data)+1]*g_window(len(data),minfreq+i,factor))
    return st_res


def g_window(length,freq,factor):
    #origin gausswindow function,i don't understand
    vector=np.zeros((2,length))
    for i in range(length):
        vector[0,i]=i
        vector[1,i]=-length+i
    vector=np.power(vector,2)
    vector=vector*(-factor*2*np.pi**2/freq**2)
    return sum(np.exp(vector))
	
	
# if __name__=='__main__':
    # t=np.linspace(0,10,5001)
    # w=signal.chirp(t,f0=12.5,f1=2.5,t1=10,method='linear')
    # fmin=0
    # fmax=250
    # stres=st(w,fmin,fmax)
    # dt=t[1]-t[0]
    # fmin=fmin/(len(w)*dt)
    # fmax=fmax/(len(w)*dt)
    # extent=(t[0],t[-1],fmin,fmax)
    
    # import matplotlib.pyplot as plt
    
    # fig,ax=plt.subplots(2,1,sharex=True)
    # ax[0].plot(t,w)
    # ax[1].imshow(np.abs(stres),origin='lower',extent=extent)
    # ax[1].axis('tight')
    # ax[1].set(xlabel='time',ylabel='frequency')
    # plt.show()