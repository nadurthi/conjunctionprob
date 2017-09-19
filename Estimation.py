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
import copy
import cvxpy as cxpy
import time

print "loaded Estimation module on : ",pd.datetime.today()


def getrawmoms(Nmoms,mu=0,sig=1):
    X=np.zeros(Nmoms)
    for i in range(0,Nmoms):
        C=np.abs( sc.special.hermitenorm(i,True) )
        C=[int(c) for c in C]
        x=np.repeat(mu/sig, len(C))
        p=range(len(C)-1,-1,-1)
        X[i]=sig**(len(C)-1)*np.sum( np.multiply(C,np.power(x,p)) )
    return X


def coloriter():
    for m in ['o','s','*','+']:
        for c in ['b','r','g','k','y']:
            yield c+m


def plotellipse(MU,SIG):
    MU=MU.reshape(-1,1)
    th=np.linspace(0,2*np.pi,100)
    X=np.zeros((100,2))
    for i in range(100):
        X[i,:]=( np.dot(sc.linalg.sqrtm(SIG),np.array([[np.cos(th[i])],[np.sin(th[i])]]))+MU ).transpose()[0]
    plt.plot(X[:,0],X[:,1],'r')


class Normalize0I(object):
    def __init__(self,X,scale=1,w=None):
        
        M,P=GetMeanCov(X,w=w)
        self.mean=M
        self.cov=scale*P
        
    def normalize(self,X,probs=None):
        P=self.cov
        M=self.mean
        invsqrtP=np.linalg.inv( sc.linalg.sqrtm(P) )
        M=M.reshape(-1,1)
        
        Y=np.zeros(X.shape)
        
        for i in range(X.shape[0]):
            Y[i,:]=np.dot(invsqrtP,X[i,:].reshape(-1,1)-M).reshape(1,-1)[0]

        if probs is not None:
            probs=probs/np.abs( np.linalg.det(invsqrtP) )
            return Y,probs
        else:
            return Y
    
    def revert(self,X,probs=None):
        P=self.cov
        M=self.mean
        sqrtP=sc.linalg.sqrtm(P) 
        M=M.reshape(-1,1)
        Y=np.zeros(X.shape)
        for i in range(X.shape[0]):
            Y[i,:]=(np.dot(sqrtP,X[i,:].reshape(-1,1))+M).reshape(1,-1)[0]

        if probs is not None:
            probs=probs/np.abs( np.linalg.det(sqrtP) )
            return Y,probs
        else:
            return Y

class Normalize_11(object):
    def __init__(self,X,w=None):
        
        M,P=GetMeanCov(X,w=w)
        self.mean=M
        Y=X-np.tile(M.reshape(1,-1)[0],(X.shape[0],1))
        Ymn=np.amax(np.abs(Y),axis=0)
        P=np.diag(np.power(Ymn,2))
#         Ymx=np.amax(Y,axis=0)
        
        self.cov=P
    
    def normalize(self,X,probs=None):
        P=self.cov
        M=self.mean
        invsqrtP=np.linalg.inv( sc.linalg.sqrtm(P) )
        M=M.reshape(-1,1)
        
        Y=np.zeros(X.shape)
        
        for i in range(X.shape[0]):
            Y[i,:]=np.dot(invsqrtP,X[i,:].reshape(-1,1)-M).reshape(1,-1)[0]

        if probs is not None:
            probs=probs/np.abs( np.linalg.det(invsqrtP) )
            return Y,probs
        else:
            return Y
    
    def revert(self,X,probs=None):
        P=self.cov
        M=self.mean
        sqrtP=sc.linalg.sqrtm(P) 
        M=M.reshape(-1,1)
        Y=np.zeros(X.shape)
        for i in range(X.shape[0]):
            Y[i,:]=(np.dot(sqrtP,X[i,:].reshape(-1,1))+M).reshape(1,-1)[0]

        if probs is not None:
            probs=probs/np.abs( np.linalg.det(sqrtP) )
            return Y,probs
        else:
            return Y


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

def fitclusters(Z,NN=2,w=None):
    M,P=GetMeanCov(Z,w=w)
    eigs,v=np.linalg.eig(P)
    print eigs
    maxeig=np.sqrt( max(eigs) )
    mineig=np.sqrt( min(eigs) )
    
#     print "main clsuter"
#     print M
#     print maxeig
    
    savekmeans={}
    flg=0
    for Ncls in range(1,NN+1):
        maxeig
        kmeans = KMeans(n_clusters=Ncls, random_state=10).fit(Z)
        Zids=kmeans.labels_
        # get distance between clusters
        d=[]
        for i in range(Ncls):
            for j in range(i+1,Ncls):
                d.append( np.linalg.norm(kmeans.cluster_centers_[i]-kmeans.cluster_centers_[j]) )
        
        savekmeans[Ncls]=kmeans
        print "cluster split ------------==========---------============-------"
        print d,maxeig,mineig
        if len(d)!=0:
            if min(d)< maxeig:
                N=Ncls-1
                flg=1
                break
    
    if flg==0:
        N=Ncls
    kmeans=savekmeans[N]
    # split the points into clusters

    return (N,kmeans.labels_,kmeans.cluster_centers_)

def weightcost(x,A,b):
    return np.linalg.norm( np.dot(A,x.reshape(-1,1))-b )**2




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


class GaussMixModel(object):
    def __init__(self,ws,Ms,Ps,**kwargs):
        self.weights=copy.deepcopy( list(ws) )
        self.means=np.array(Ms).copy()
        self.covs=np.array(Ps).copy()
        self.N=len(ws)
        self.illcond=set([])
        
        self.appendedstate={}
        
        for key,value in kwargs.items():
            setattr(self,key,value)
    
    def setfixedpoints(self,X,px):
        self.fixedpoints=X.copy()
        self.fixedprobs=px.copy()

    


    def MakeCopy(self):
        GMM=GaussMixModel(self.weights,self.means,self.covs,appendedstate=self.appendedstate)
        GMM.fixedpoints=self.fixedpoints.copy()
        GMM.fixedprobs=self.fixedprobs.copy()

        return GMM


    def scalecovs(self,c):
        for i in range(self.N):
            self.covs[i]=self.covs[i]*c
    
    def getmarginalized_moms(self,Nmoms):
        dim=len(self.means[0])
        M=np.zeros((dim,Nmoms))
        for d in range(dim):
            for i in range(self.N):
                M[d,:]=M[d,:]+self.weights[i]*getrawmoms(Nmoms,mu=self.means[i][d],sig=self.covs[i][d,d]) 
        return M

    def Gaussiantransform(self,M,P,inplace=False,fixedpts_also=True):
        b=M.reshape(-1,1)
        A=sc.linalg.sqrtm(P)
        return self.lineartransform(A,b,inplace=inplace,fixedpts_also=fixedpts_also)

    def lineartransform(self,A,b,inplace=False,fixedpts_also=True):
        """
        linear trasnform the GMM
        """
        if fixedpts_also:
            Npoints=self.fixedpoints.shape[0]
            newpoints=np.zeros(self.fixedpoints.shape)
            newprobs=np.zeros(Npoints)

            for i in range(Npoints):
                newpoints[i,:]=(np.dot( A,self.fixedpoints[i,:].reshape(-1,1))+b.reshape(-1,1)).reshape(1,-1)[0]
                newprobs[i]=self.fixedprobs[i]/np.abs(np.linalg.det(A))

        if inplace:
            for i in range(self.N):
                self.means[i]=(np.dot(A,self.means[i].reshape(-1,1))+b.reshape(-1,1)).reshape(1,-1)[0]
                self.covs[i]=np.dot(np.dot(A,self.covs[i]),A.transpose())

            if fixedpts_also:            
                self.fixedpoints=newpoints.copy()
                self.fixedprobs=newprobs.copy()


        else:

            GMM=self.MakeCopy()

            for i in range(self.N):
                GMM.means[i]=(np.dot(A,self.means[i].reshape(-1,1))+b.reshape(-1,1)).reshape(1,-1)[0]
                GMM.covs[i]=np.dot(np.dot(A,self.covs[i]),A.transpose())
            
            if fixedpts_also:
                GMM.fixedpoints=newpoints.copy()
                GMM.fixedprobs=newprobs.copy()
                
            return GMM

    def functiontransform(self,Func,Fjac=None,method='UT',inplace=False,fixedpts_also=True):
        """
        linear trasnform the GMM
        """
        if fixedpts_also and Fjac is not None:
            Npoints=self.fixedpoints.shape[0]
            newpoints=np.zeros(self.fixedpoints.shape)
            newprobs=np.zeros(Npoints)

            for i in range(Npoints):
                A=Fjac(self.fixedpoints[i,:])
                newpoints[i,:]=Func(self.fixedpoints[i,:])
                newprobs[i]=self.fixedprobs[i]/np.abs(np.linalg.det(A))

        if inplace:
            for i in range(self.N):
                X,w=UnscentedTransform_pts(self.means[i],self.covs[i])
                Y=np.zeros(X.shape)
                for j in range(X.shape[0]):
                    Y[j,:]=Func(X[j,:])
                mu,P=GetMeanCov(Y,w=w)
                self.means[i]=mu
                self.covs[i]=P

                if fixedpts_also and Fjac is not None:            
                    self.fixedpoints=newpoints.copy()
                    self.fixedprobs=newprobs.copy()

        else:

            GMM1=GaussMixModel(self.weights,self.means,self.covs)
            for i in range(GMM1.N):
                X,w=UnscentedTransform_pts(GMM1.means[i],GMM1.covs[i])
                Y=np.zeros(X.shape)
                for j in range(X.shape[0]):
                    Y[j,:]=Func(X[j,:])
                mu,P=GetMeanCov(Y,w=w)
                GMM1.means[i]=mu
                GMM1.covs[i]=P

            if fixedpts_also and Fjac is not None:
                GMM1.fixedpoints=newpoints.copy()
                GMM1.fixedprobs=newprobs.copy()

            return GMM1

    def purgecomponent(self,i):
#         print "purging i = ",i," with mean ",self.means[i]," and cov ",self.covs[i]
        self.weights=np.delete(self.weights,i)
        self.means=np.delete(self.means,i,axis=0)
        self.covs=np.delete(self.covs,i,axis=0)
        self.N=self.N-1
        
    def AppendComponents(self,ws,Ms,Ps):
        if len(ws)==0 or len(Ms)==0 or len(Ps)==0:
            return
        self.appendedstate={'ws':ws,'Ms':Ms,'Ps':Ps}
        
        self.weights=list(self.weights)+list(ws)
        self.means=np.vstack(( self.means,np.array(Ms).copy() ))
        self.covs=np.vstack(( self.covs,np.array(Ps).copy() ))
        self.N=len(self.weights)
    
    def ReplaceComponents(self,ws,Ms,Ps):
        
        self.weights=copy.deepcopy( list(ws) )
        self.means=np.array(Ms).copy()
        self.covs=np.array(Ps).copy()
        self.N=len(ws)
        self.illcond=set([])
        
        self.appendedstate={}


    def evalaute_probability(self,X):

        pp=np.zeros(X.shape[0])
        p=self.get_component_probabilities(X)
        
        for i in range(self.N):
            pp=pp+self.weights[i]*p[i]

        return pp
    
    def evalaute_probability_at_fixedpoints(self):
        return self.evalaute_probability(self.fixedpoints)
    
    def get_component_probabilities(self,X,USE0I=True):
        p=[]
        for i in range(self.N):
            if USE0I:
                pp=mvnpdf(X,self.means[i],self.covs[i])
            else:
                try:
                    pp=multivariate_normal.pdf(X, self.means[i], self.covs[i],allow_singular=False)

                except:
                    print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
                    print "                SINGULAR COMPONENT ", i,"           "
                    print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

                    pp=mvnpdf(X,self.means[i],self.covs[i])
            
            p.append( pp)
            
        return np.array(p)
    
    def Prune_byweight(self):
        stdw=np.std(self.weights)
        mnw=np.mean(self.weights)
        
        while(1):
            flg=0
            for i in range(self.N):
                if self.weights[i]<1e-25:
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
    
    def Prune_LowRank(self):
        eigs=np.abs( np.array([np.linalg.eig(P)[0] for P in self.covs]) )
        eigmean=np.mean(eigs,axis=0)
        eigstd=np.std(eigs,axis=0)
        
        while(1):
            flg=0
            for i in range(self.N):
                if np.linalg.matrix_rank(self.covs[i])<self.covs[i].shape[0]:
                    self.purgecomponent(i)
                    eigs=np.delete(eigs,i,axis=0)
                    
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
    
    def normalize_weights(self):
        weights=np.ones(self.N)
        weights=weights/np.sum(weights)
        self.weights=weights
        
        return weights
    
    def optimize_weights_prob(self,X,pT):
        p0=pT
        p=self.get_component_probabilities(X)
        A=p.transpose()
        b=p0


        results=sc.optimize.lsq_linear(A,b,bounds= (np.zeros(self.N), np.ones(self.N)) )
        weights=results.x
        weights=weights/np.sum(weights)
        self.weights=weights
        
        return weights
    
    def optimize_weights_prob_nonlinear(self,X,pT):
        p0=pT
        p=self.get_component_probabilities(X)
        A=p.transpose()
        b=p0
        cons = ({'type': 'eq',
                'fun' : lambda x: np.sum(x)-1,
                }) # 'jac' : lambda x: np.ones(len(x)) (np.zeros(self.N), np.ones(self.N))
        results=sc.optimize.minimize(weightcost, np.ones(self.N)/self.N,args=(A,b),constraints=cons,bounds=[(0,1) for i in range(self.N)] )
        weights=results.x
        weights=weights/np.sum(weights)
        self.weights=weights
        
        return weights
    
    def optimize_weights_cvx(self,X,pT):

        
        p0=pT
        p=self.get_component_probabilities(X)
        A=p.transpose()
        b=p0
        
        # Construct the problem.
        x = cxpy.Variable(self.N)
        objective = cxpy.Minimize(cxpy.norm(A*x - b,1))
        constraints = [0 <= x, x <= 1,cxpy.sum_entries(x)==1]
        prob = cxpy.Problem(objective, constraints)

        # The optimal objective is returned by prob.solve().
        result = prob.solve(solver='CVXOPT')
        if np.isfinite(result):
            # The optimal value for x is stored in x.value.
            self.weights=x.value
            return x.value
        else:
            return self.optimize_weights_prob(X,pT)
        # The optimal Lagrange multiplier for a constraint
        # is stored in constraint.dual_value.
#         print constraints[0].dual_value

    def plotcomp_ellipsoids(self,k,dims=[0,1],comp=None,ax=None,c='r'):
        if comp is None:
            rang=range(self.N)
        else:
            rang=[comp]
        for i in rang:
            MU=self.means[i].reshape(-1,1)[[dims[0],dims[1] ]]
            SIG=self.covs[i][[[ dims[0] ],[dims[1]]],[dims[0],dims[1]]]
            sqrtmSIG=sc.linalg.sqrtm(k*SIG)
            th=np.linspace(0,2*np.pi,100)
            X=np.zeros((100,2))
            for j in range(100):
                X[j,:]=( np.dot(sqrtmSIG,np.array([[np.cos(th[j])],[np.sin(th[j])]]))+MU ).transpose()[0]
            if ax is None:
                plt.plot(X[:,0],X[:,1],c)
            else:
                ax.plot(X[:,0],X[:,1],c)
                

def MakeGMMcopy(GMM):
    return GaussMixModel(GMM.weights,GMM.means,GMM.covs)




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
    
    Xclass=pdiff
    Xclass_errors={}
    for leafid in leafnodes:
        Y=Xclass[Xleafids==leafid]
        Xclass_errors[leafid]=np.average(Y)
        
    return (clf,nodepaths,leafnodes,Xleafids,Xclass,Xclass_errors)


def getpartitionMeanCovs(X,Xpartids,probs=None):
    newmeans=[]
    newcovs=[]
    leafids=np.unique(Xpartids)
    for lfn in leafids:
        XX=X[Xpartids==lfn,:]
        if len(XX)==0:
            continue
            
        if probs is not None:
            w=probs[Xpartids==lfn].copy()
            w=w/np.sum(w)
        else:
            w=None
            
        mm,cc=GetMeanCov(XX,w=w)
        A=[]
        meanprob=np.max(probs)
        if meanprob==0:
            continue
            
        for c in np.hstack((np.linspace(1e-5,2,100))):
#             pp=multivariate_normal.pdf(mm, mm, c*cc)
            pp=1/np.sqrt(np.linalg.det(2*np.pi*c*cc))
            A.append( (c,np.linalg.norm(pp-meanprob),pp,cc,1/np.sqrt(np.linalg.det(2*np.pi*cc))) )
        
#         for c in np.hstack((np.linspace(1e-15,5,1000),np.linspace(5,50,100))):
#             cq=np.identity(X.shape[1])
#             pq=1/np.sqrt(np.linalg.det(2*np.pi*c*cq))
#             A.append( (c,np.linalg.norm(pq-meanprob),pq,cq) )
        
        
        try:
            S=sorted(A,key=lambda x:x[1])
            c=S[0][0]
        except:
            print A
            print mm,cc
        c=S[0][0]
        pp=S[0][2]
        cc=S[0][3]
        
        
        print "meanprob = ",meanprob," c = ",c,pp ," max eig cc = ",max(np.linalg.eig(cc)[0]) ," cc0 = ",S[0][4]
        
        if c>20:
            continue


        cc=c*cc
   

        if np.linalg.matrix_rank(cc)<cc.shape[0]:
            continue
            
            eig,_=np.linalg.eig(cc)
            eig=np.abs(eig)
            meaneig=sorted(eig)[int(len(eig)/2)]
            cc=np.diag(np.ones(cc.shape[0])*meaneig)
            
        newmeans.append( mm )
        newcovs.append( cc )
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
    



def mvnpdf0I(X):
    dim=X.shape[1]
    M=np.zeros(dim)
    P=np.identity(dim)
    pp=multivariate_normal.pdf(X, M, P)

    return pp

def mvnpdf(X,M,P):
    
    invsqrtP=np.linalg.inv( sc.linalg.sqrtm(P) )
    Y=np.zeros(X.shape)
    for j in range(X.shape[0]):
        Y[j,:]=np.dot(invsqrtP,(X[j,:]-M).reshape(-1,1)).reshape(1,-1)
    pp=multivariate_normal.pdf(Y, np.zeros(X.shape[1]), np.identity(X.shape[1]))
    pp=pp*np.abs(np.linalg.det(invsqrtP))
    
    return pp


def mvnrnd(M,P,N=1):
    X=np.random.multivariate_normal(M.reshape(1,-1)[0], P, N)
    return X

def mvnrnd0I(dim,N=1):
    X=np.random.multivariate_normal(np.zeros(dim), np.identity(dim), N)
    return X


def getthresplots(pdiff,typ='Rel'):
    if typ=='Rel':
#         pdiff=100*np.abs( np.divide((py-pyest),py) )
        fig,ax=plt.subplots()
        h=ax.hist(np.abs(pdiff),np.linspace(np.percentile(np.abs(pdiff),20)/2,2*np.percentile(np.abs(pdiff),80),100))
        ax.set_title('Rel error')

    else:
#         pdiff=100*np.abs( np.divide((py-pyest),1) )
        fig,ax=plt.subplots()
        h=ax.hist(np.abs(pdiff),np.linspace(np.percentile(np.abs(pdiff),20)/2,2*np.percentile(np.abs(pdiff),80),100))
        ax.set_title('Abs error')


def getclusterIDs(X,means,ClusterIds=None):
    """
    given the means and their ids in clusterid, 
    .... cluster X into those ids
    """
    if ClusterIds is None:
        ClusterIds=np.arange(means.shape[0])

    Xids=np.zeros(X.shape[0])
    Y=None
    for j in range(means.shape[0]):
        pp=np.linalg.norm( X-np.tile(means[j,:],(X.shape[0],1)),axis=1 ).reshape(-1,1)
        if Y is None:
            Y=pp
        else:
            Y=np.hstack((Y,pp))

    dd=np.argmin(Y,axis=1)    

    Xids=ClusterIds[dd]
    return Xids