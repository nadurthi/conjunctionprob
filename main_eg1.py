import numpy as np
import Estimation as estmn
import matplotlib.pyplot as plt


"""
Example 1 of polar coordinate transformation
1. Everythinh is done in Y space
2. Past GMM covariance is not updated, only ht eweights are updated

"""

        
def Frev(X):
    r=np.linalg.norm(X)
    th=np.arctan2(X[1],X[0])
    return np.array([r,th])


def FJac(rth):
    th=rth[1]
    r=rth[0]
    return np.array( [[np.cos(th),-r*np.sin(th)],[np.sin(th),r*np.cos(th)]] )


def F_polar2cart(rth):
    x=rth[0]*np.cos(rth[1])
    y=rth[0]*np.sin(rth[1])
    return np.array([x,y])

def FJac_polar2cart(rth):
    th=rth[1]
    r=rth[0]
    return np.array( [[np.cos(th),-r*np.sin(th)],[np.sin(th),r*np.cos(th)]] )

    
Mu0=np.array([[30],[np.pi/2]])
P0=np.array( [[2**2,0],[0,(np.pi/6)**2]] )
N=1000
X0=estmn.mvnrnd(Mu0, P0, N)

N=X0.shape[0]
GMMX=estmn.GaussMixModel([1],[Mu0.reshape(1,-1)[0]],[P0])

Z=np.zeros((N,2))
px=GMMX.evalaute_probability(X0)
pz=np.zeros(px.shape)

for i in range(N):
    Z[i,:]=estmn.F_polar2cart(X0[i,:])
    pz[i]=px[i]/np.linalg.det( estmn.FJac_polar2cart(X0[i,:]))

# Y=Z.copy()
# py=pz.copy()

NM11=estmn.Normalize0I(Z)
Y,py=NM11.normalize(Z,probs=pz)


MuY,PY=estmn.GetMeanCov(Y)
GMMY=estmn.GaussMixModel([1],[MuY.reshape(1,-1)[0]],[PY])

GMMHistory=[]

# py=py/100
perctY=90
pyest=GMMY.evalaute_probability(Y)
pdiff=100*np.abs( np.divide((py-pyest),(py) ) )

# pdiff=py

thres=np.percentile(np.abs(pdiff),perctY)
print thres
fig,ax=plt.subplots()
h=ax.hist(np.abs(pdiff),np.linspace(0,2*np.percentile(np.abs(pdiff),80),100))

prevthresAbs=0
prevthresRel=0


Ncnt=200
Nsamples=[int(a) for a in np.linspace(100,10,Ncnt)]
for cnt in range(Ncnt):
    print "*****************************************************************"
    # first transform the initial X
    


    clf,nodepaths,leafnodes,Xleafids,Xclass,Xclass_errors= estmn.getclassifiedregions(pdiff,thres,Y,rndT=10+cnt,min_samples_leaf=Nsamples[cnt],max_depth=6)
    print "------Leaf errors --------- "
    print Xclass_errors
    print "--------------------------- "


    # get the leaf with max error
#     LeafmaxErrorId=0
#     mlf=0.0
#     for key,value in Xclass_errors.items():
#         if value >mlf:
#             LeafmaxErrorId=key
#             mlf=value
    
    LeafmaxErrorId=sorted(Xclass_errors.items(),key=lambda x: x[1],reverse=True)[0][0]
    
    # split this leaf using kmeans
    YY=Y[Xleafids==LeafmaxErrorId,:].copy()
    ppY=py[Xleafids==LeafmaxErrorId].copy()
    NYclusts,Yclustids,A=estmn.fitclusters(YY,NN=5,w=None)
    Ynewmeans,Ynewcovs=estmn.getpartitionMeanCovs(YY,Yclustids,probs=ppY)
    NYclusts=len(Ynewmeans)
    print "NYclusts = ",NYclusts
    print "prev N = ",GMMY.N
    wnew=np.ones(NYclusts)/NYclusts
#     print "append w = ",wnew
    print "appending means ",Ynewmeans
    print "old N = ",GMMY.N,len(GMMY.means)
    
#     print GMMY.weights
#     GMMHistory.append( (0,GMMY.MakeCopy()) )
#     print GMMY.weights
    
    GMMY.AppendComponents(wnew,Ynewmeans,Ynewcovs )
#     break
    print GMMY.optimize_weights_prob(Y,py)
#     GMMY.normalize_weights()
#     GMMHistory.append( (1,GMMY.MakeCopy()) )
#     print GMMY.optimize_weights_cvx(Y,py)
    
    print "append N = ",GMMY.N,len(GMMY.means)
    
    GMMY.Prune_LowRank()
    print "purged N = ",GMMY.N,len(GMMY.means)
    if cnt%5==0:
#         GMMY.Prune_byCov()
        GMMY.Prune_byweight()

#     GMMHistory.append( (2,GMMY.MakeCopy()) )
    
#     GMMY.weights=np.ones(GMMY.N)/GMMY.N
#     GMMY.optimize_weights_shitty(Y,py)
    
    


    pyest=GMMY.evalaute_probability(Y)
    if cnt%2==0: #
        pdiff=100*np.abs( np.divide((py-pyest),py) )
        estmn.getthresplots(pdiff,typ='Rel')
        thres=np.percentile(np.abs(pdiff),perctY)
        print "Rel thres= ",thres
#         if cnt>2 and np.abs(prevthresRel-thres)/prevthresRel>0.3 and thres>prevthresRel:
#             print "failed relative error"
#             break
            
        prevthresRel=thres

    else:
        pdiff=100*np.abs( np.divide((py-pyest),1) )
        estmn.getthresplots(pdiff,typ='Abs')
        thres=np.percentile(np.abs(pdiff),perctY)
        print "Abs thres= ",thres
#         if cnt>2 and np.abs(prevthresAbs-thres)/prevthresAbs>0.3 and thres>prevthresAbs:
#             print "failed absolute error"
#             break
        prevthresAbs=thres
    
    

    
fig,ax=plt.subplots(1,3,figsize=(20,7))
ax[0].plot(Y[:,0],Y[:,1],'bo')
ax[0].plot(YY[:,0],YY[:,1],'g*')
GMMY.plotcomp_ellipsoids(1,dims=[0,1],ax=ax[0])

ax[1].plot(X0[:,0],X0[:,1],'ro')
XX=X0[Xleafids==LeafmaxErrorId,:]
CC=estmn.coloriter()
ax[1].plot(XX[Yclustids==0,0],XX[Yclustids==0,1],'ks')
ax[1].plot(XX[Yclustids==1,0],XX[Yclustids==1,1],'g+')
    
#     ax[2].plot(YY[:,0],YY[:,2],'g*')
    



plt.show()
#     break