from __future__ import division
import matplotlib.pyplot as plt
from sklearn import tree
import scipy as sc
import numpy as np
from scipy.stats import multivariate_normal
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.cluster import KMeans

def coloriter():
    for m in ['o','s','*','+']:
        for c in ['b','r','g','k','y']:
            yield c+m

def F(rth):
    x=rth[0]*np.cos(rth[1])
    y=rth[0]*np.sin(rth[1])
    return np.array([x,y])

def Normalize01(X,w=None,probs=None):
    Y=np.zeros(X.shape)
    M,P=GetMeanCov(X,w=w)
    invsqrtP=np.linalg.inv( sc.linalg.sqrtm(P) )
    M=M.reshape(-1,1)
    for i in range(X.shape[0]):
        Y[i,:]=np.dot(invsqrtP,X[i,:].reshape(-1,1)-M).reshape(1,-1)[0]
    
    if probs is not None:
        probs=probs/np.abs( np.linalg.det(invsqrtP) )
        return Y,probs
    else:
        return Y

def Frev(X):
    r=np.linalg.norm(X)
    th=np.arctan2(X[1],X[0])
    return np.array([r,th])


def FJac(rth):
    th=rth[1]
    r=rth[0]
    return np.array( [[np.cos(th),-r*np.sin(th)],[np.sin(th),r*np.cos(th)]] )


def GetMeanCov(X,w=None):
    M=np.average(X,axis=0,weights=w)
#     P=np.cov(X,rowvar=False,aweights=w)
    PP=0
    if w is None:
        w=np.ones(X.shape[0])/X.shape[0]
        
    for i in range(X.shape[0]):
        PP=PP+w[i]*np.dot((X[i,:]-M).reshape(-1,1),(X[i,:]-M).reshape(1,-1) )
    
    return (M,PP)

def UnscentedTransform_pts(mu,P):
    mu=mu.copy().reshape(1,-1)[0]
    n=P.shape[0];
    X=np.zeros((2*n+1,n))
    w=np.zeros(2*n+1)
    if n<4:
        k=3-n;
    else:
        k=1;

    X[0,:]=mu;
    w[0]=k/(n+k);
    A=sc.linalg.sqrtm((n+k)*P);
    for i in range(0,n):
        X[i+1,:]=(mu+A[:,i])
        X[i+n+1,:]=(mu-A[:,i])
        w[i+1]=1/(2*(n+k));
        w[i+n+1]=1/(2*(n+k));
    
    return (X,w)

def fitclusters(Z,w=None):
    M,P=GetMeanCov(Z,w=w)
    eigs,v=np.linalg.eig(P)
    print eigs
    maxeig=np.sqrt( max(eigs) )
    
    print "main clsuter"
    print M
    print maxeig
    
    savekmeans={}
    for Ncls in [2,3,4]:
        maxeig
        kmeans = KMeans(n_clusters=Ncls, random_state=10).fit(Z)
        Zids=kmeans.labels_
        # get distance between clusters
        d=[]
        for i in range(Ncls):
            for j in range(i+1,Ncls):
                d.append( np.linalg.norm(kmeans.cluster_centers_[i]-kmeans.cluster_centers_[j]) )
        
        savekmeans[Ncls]=kmeans
        print "d = ",d
        if min(d)< maxeig:
            break
            
    print "#cluster = ",Ncls-1
    
    N=Ncls-1
    kmeans=savekmeans[N]
    # split the points into clusters

    return (N,kmeans.labels_,kmeans.cluster_centers_)

class GaussMixModel(object):
    def __init__(self,ws,Ms,Ps):
        self.weights=list(ws)
        self.means=np.array(Ms).copy()
        self.covs=np.array(Ps).copy()
        self.N=len(ws)
    
    def purgecomponent(self,i):
#         print "purging i = ",i," with mean ",self.means[i]," and cov ",self.covs[i]
        del self.weights[i]
        self.means=np.delete(self.means,i,axis=0)
        self.covs=np.delete(self.covs,i,axis=0)
        self.N=self.N-1
        
    def AppendComponents(self,ws,Ms,Ps):
        self.weights=list(self.weights)+list(ws)
        self.means=np.vstack(( self.means,np.array(Ms).copy() ))
        self.covs=np.vstack(( self.covs,np.array(Ps).copy() ))
        self.N=len(self.weights)
    
    def evalaute_probability(self,X):
        p=np.zeros(X.shape[0])
        for i in range(self.N):
            p=p+self.weights[i]*multivariate_normal.pdf(X, self.means[i], self.covs[i])
        return p
    def get_component_probabilities(self,X):
        p=[]
        for i in range(self.N):
            p.append( multivariate_normal.pdf(X, self.means[i], self.covs[i]) )
        return np.array(p)
    
    def Prune_byweight(self):
        stdw=np.std(self.weights)
        mnw=np.mean(self.weights)
        
        while(1):
            flg=0
            for i in range(self.N):
                if np.abs(self.weights[i]-mnw)>3*stdw or self.weights[i]<1e-25:
                    self.purgecomponent(i)
                    flg=1
                    break
            if flg==0:
                break
    
    def Prune_byCov(self):
        eigs=np.abs(np.array( [max(np.linalg.eig(P)[0]) for P in self.covs]) )
        stdeig=np.std(eigs)
        mnw=np.mean(eigs)
        
        while(1):
            flg=0
            for i in range(self.N):
                if np.abs(eigs[i]-mnw)>3*stdeig or np.linalg.matrix_rank(self.covs[i])<self.covs[i].shape[0]:
                    self.purgecomponent(i)
                    eigs=np.delete(eigs,i)
                    flg=1
                    break
            if flg==0:
                break
        
    def getweight_rangeid(self,w):
        W={}
        for i in range(0,self.N):
            if i==0:
                W[0]=(0,self.weights[0])
            else:
                W[i]=( W[i-1][1],W[i-1][1]+self.weights[i])
            if w>=W[i][0] and w<W[i][1]:
                return i
            
#         D={}
#         for key,value in W.items():
#             D[value]=key
        
        
    def randomsamples(self,N,scale=1):
        
        X=np.zeros((N,self.covs[0].shape[0]))
        for i in range(N):
            w=np.random.random(1)
            ind=self.getweight_rangeid(w)
#             print w,ind
            X[i,:]=np.random.multivariate_normal(self.means[ind], scale**2*self.covs[ind], 1)
        return X
        
        
    def optimize_weights_original(self,X,M0,P0):
        M0=M0.reshape(1,-1)[0]
        p0=multivariate_normal.pdf(X, M0, P0)
        p=self.get_component_probabilities(X)
        A=p.transpose()
        b=p0
        weights=np.linalg.lstsq(A, b)[0]
        weights=weights/np.sum(weights)
        self.weights=weights
        
        return weights
    
    def optimize_weights_prob(self,X,pT):
        p0=pT
        p=self.get_component_probabilities(X)
        A=p.transpose()
        b=p0
        results=sc.optimize.lsq_linear(A, b,bounds=(np.zeros(self.N), np.ones(self.N)))
        weights=results.x
        weights=weights/np.sum(weights)
        self.weights=weights
        
        return weights
    
    def plotcomp_ellipsoids(self,k,dims=[0,1],comp=None,ax=None):
        if comp is None:
            rang=range(self.N)
        else:
            rang=[comp]
        for i in rang:
            MU=self.means[i].reshape(-1,1)[0:2]
            SIG=self.covs[i][0:2,0:2]
            sqrtmSIG=sc.linalg.sqrtm(k*SIG)
            th=np.linspace(0,2*np.pi,100)
            X=np.zeros((100,2))
            for j in range(100):
                X[j,:]=( np.dot(sqrtmSIG,np.array([[np.cos(th[j])],[np.sin(th[j])]]))+MU ).transpose()[0]
            if ax is None:
                plt.plot(X[:,0],X[:,1],'r')
            else:
                ax.plot(X[:,0],X[:,1],'r')
                
        
def transformGMM(GMM0,Func):
    GMM1=GaussMixModel(GMM0.weights,GMM0.means,GMM0.covs)
    for i in range(GMM0.N):
        X,w=UnscentedTransform_pts(GMM0.means[i],GMM0.covs[i])
        Y=np.zeros(X.shape)
        for j in range(X.shape[0]):
            Y[j,:]=Func(X[j,:])
        mu,P=GetMeanCov(Y,w=w)
        GMM1.means[i]=mu
        GMM1.covs[i]=P
    return GMM1

def getleafpaths(estimator):
    n_nodes = estimator.tree_.node_count
    children_left = estimator.tree_.children_left
    children_right = estimator.tree_.children_right
    feature = estimator.tree_.feature
    threshold = estimator.tree_.threshold


    # The tree structure can be traversed to compute various properties such
    # as the depth of each node and whether or not it is a leaf.
    node_depth = np.zeros(shape=n_nodes, dtype=np.int64)
    is_leaves = np.zeros(shape=n_nodes, dtype=bool)
    stack = [(0, -1)]  # seed is the root node id and its parent depth
    while len(stack) > 0:
        node_id, parent_depth = stack.pop()
        node_depth[node_id] = parent_depth + 1

        # If we have a test node
        if (children_left[node_id] != children_right[node_id]):
            stack.append((children_left[node_id], parent_depth + 1))
            stack.append((children_right[node_id], parent_depth + 1))
        else:
            is_leaves[node_id] = True
    
    nodepath={}
    leafnodes=[]
#     print n_nodes
    for i in range(n_nodes):
        nodepath[i]=[]
        if is_leaves[i]:
            leafnodes.append(i)
            
    for i in range(n_nodes):
            
        if is_leaves[i]:
#             print("%snode=%s leaf node." % (node_depth[i] * "\t", i))
            pass
        else:
            L=nodepath[i]+[(feature[i],'lte',threshold[i])]
            R=nodepath[i]+[(feature[i],'gt',threshold[i])]
            
            nodepath[ children_left[i]  ]=L
            nodepath[ children_right[i] ]=R
            
    del nodepath[0]
    return (nodepath,leafnodes)

def getleafids(X,nodepaths,leafnodes):
#     MajorIndex=np.array(range(X.shape[0]))
    Xleafids=np.zeros(X.shape[0])
    for lfn in leafnodes:
        path=nodepaths[lfn]
        ind=np.ones(X.shape[0], dtype=bool)
        for pp in path:
            if pp[1]=='lte':
                ind=ind & (X[:,pp[0]]<=pp[2])
            else:
                ind=ind & (X[:,pp[0]]>pp[2])
        Xleafids[ind]=lfn
    return Xleafids

def getclassifiedregions(pdiff,thres,X,rndT=0,min_samples_leaf=100,max_depth=2):
    err=pdiff.copy()
    indhigh=pdiff>=thres
    indlow=pdiff<thres
    err[indhigh]=1
    err[indlow]=0

    clf = tree.DecisionTreeClassifier(max_features=1,random_state=rndT,criterion='entropy',class_weight='balanced',min_samples_leaf=min_samples_leaf,max_depth=max_depth)
    clf = clf.fit(X, err)
    Errpred=clf.predict(X)
#     print np.sum( np.abs(Errpred-err) )
    
    nodepaths,leafnodes=getleafpaths(clf)
    Xleafids=getleafids(X,nodepaths,leafnodes)
    
    Xclass=err
    Xclass_errors={}
    for leafid in leafnodes:
        Y=Xclass[Xleafids==leafid]
        Xclass_errors[leafid]=np.average(Y)
        
    return (clf,nodepaths,leafnodes,Xleafids,Xclass,Xclass_errors)


def getpartitionMeanCovs(X,Xpartids):
    newmeans=[]
    newcovs=[]
    leafids=np.unique(Xpartids)
    for lfn in leafids:
        XX=X[Xpartids==lfn,:]
        mm,cc=GetMeanCov(XX)
        if np.linalg.matrix_rank(cc)<cc.shape[0]:
            eig,_=np.linalg.eig(cc)
            eig=np.abs(eig)
            meaneig=sorted(eig)[int(len(eig)/2)]
            cc=np.diag(np.ones(cc.shape[0])*meaneig)
            
        newmeans.append( mm )
        newcovs.append( 1.0*cc )
    return newmeans,newcovs



def printtree(estimator):
# estimator=clf
    n_nodes = estimator.tree_.node_count
    children_left = estimator.tree_.children_left
    children_right = estimator.tree_.children_right
    feature = estimator.tree_.feature
    threshold = estimator.tree_.threshold


    # The tree structure can be traversed to compute various properties such
    # as the depth of each node and whether or not it is a leaf.
    node_depth = np.zeros(shape=n_nodes, dtype=np.int64)
    is_leaves = np.zeros(shape=n_nodes, dtype=bool)
    stack = [(0, -1)]  # seed is the root node id and its parent depth
    while len(stack) > 0:
        node_id, parent_depth = stack.pop()
        node_depth[node_id] = parent_depth + 1

        # If we have a test node
        if (children_left[node_id] != children_right[node_id]):
            stack.append((children_left[node_id], parent_depth + 1))
            stack.append((children_right[node_id], parent_depth + 1))
        else:
            is_leaves[node_id] = True

    print("The binary tree structure has %s nodes and has "
          "the following tree structure:"
          % n_nodes)
    for i in range(n_nodes):
        if is_leaves[i]:
            print("%snode=%s leaf node." % (node_depth[i] * "\t", i))
        else:
            print("%snode=%s test node: go to node %s if X[:, %s] <= %s else to "
                  "node %s."
                  % (node_depth[i] * "\t",
                     i,
                     children_left[i],
                     feature[i],
                     threshold[i],
                     children_right[i],
                     ))
    print()
    