"""Equivalent moment SDP; its dual variables propose the SOS certificate."""
import numpy as np
from scipy import sparse
import cvxpy as cp


def problem(free,rhs,projected,row_scale):
    z=cp.Variable(len(rhs));y=cp.multiply(row_scale,z)
    objective=np.zeros(free.shape[1]);objective[0]=1
    equation=free.T@y==objective;cones=[]
    for P,ii,jj in projected:
        r=int(max(ii.max(),jj.max()))+1
        cols=np.arange(len(ii));off=ii!=jj
        S=sparse.csc_matrix((np.concatenate((np.where(off,.5,1.),np.full(int(off.sum()),.5))),
            (np.concatenate((ii*r+jj,jj[off]*r+ii[off])),np.concatenate((cols,cols[off])))),
            shape=(r*r,len(ii)))
        C=cp.reshape(S@(P.T@y),(r,r),order='C')
        cones.append(C>>0)
    return cp.Problem(cp.Minimize(rhs@y),[equation,*cones]),equation,y,cones
