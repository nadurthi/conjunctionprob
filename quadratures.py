import numpy as np


def GenerateIndex(ND,numbasis):
    """
    % (n,n*np.ones(1,n))
    % This function computes all permutations of 1-D basis functions.
    %
    """
    numbasis=[int(c) for c in numbasis]
    print numbasis[0]
    index = np.array([np.arange(0,numbasis[0])]).reshape(-1,1)  #%short for canonical_0 - first dimension's nodes: this will be loooped through the dimensions
    for ct in np.arange(1,ND)  :
        repel = index
        repsize = len(index[:,0]); # %REPetition SIZE
        repwith = np.zeros((repsize,1)); # %REPeat WITH this structure: initialization
        
        for rs in np.arange(1,numbasis[ct]):
            repwith = np.vstack((repwith, np.ones((repsize,1))*rs))   #%update REPeating structure
        index =  np.hstack(( np.tile(repel,(numbasis[ct],1)), repwith))   #%update canon0

    return index


def Generate_conjaxis(n):
    index=GenerateIndex(n,n*np.ones(n))
    roww,coll=index.shape
    # print roww,coll
    A=np.identity(n)

    dr=None;
    for i in range(roww):
        y=index[i,:]
        if y[y>1].shape[0]==0:
            if dr is None:
                dr=index[i,:]
            else:
                dr=np.vstack((dr,index[i,:] ) )
    
    mo=-1*np.ones(n)
    X=np.zeros((2**n,n))
    for i in range(2**n):
        rr=np.power(mo,dr[i,:])
        sig=0
        for j in range(n):
            sig=sig+rr[j]*A[:,j]

        X[i,:]=sig
    return X