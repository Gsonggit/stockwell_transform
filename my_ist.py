#!/usr/bin/env python
#
# 1-D inverse stockwell transform code modified by GaoSong from origion ist.m code

import numpy as np

def ist(st_matrix):
    stsp=np.sum(st_matrix,axis=1)
    if st_matrix.shape[1] % 2 != 0:
        negsp=stsp[2:].T[::-1]
    else:
        negsp=stsp[2:-1].T[::-1]

    fullstsp=np.hstack((np.conjugate(stsp.T),negsp))
    ts=np.fft.ifft(fullstsp).real
    return ts