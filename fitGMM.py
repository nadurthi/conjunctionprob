import numpy as np
import scipy as sc
import Estimation as estmn
import quadratures as qdtr
import time

qdtr=reload(qdtr)
estmn=reload(estmn)


class GMM0I(object):
	"""
	Contract: give the dimension and the means for zero-identity
	you will get the weights,means,covs
	"""
	def __init__(self,dim=2):
		self.dim=dim



	def Getgmms(self,means):
		# first get the positive quadrant points
		for i in range(self.dim): 	
			means=means[ means[i]>=0,:]

	def Getgmms_bysamples_symmtericmeans(self,N,oldmeans,newmeans,useaxis='PA',thres=0.1):
		newmeans=np.abs(newmeans)
		


		Xnewmeans=None
		if useaxis=='PA':
			Xpa=qdtr.Generate_principalaxis(self.dim)
			for i in range(newmeans.shape[0]):
				pp=np.multiply( Xpa,np.tile(newmeans[i],(Xpa.shape[0],1)) )
				if Xnewmeans is None:
					Xnewmeans=pp
				else:
					Xnewmeans=np.vstack((Xnewmeans,pp))

		if useaxis=='CA':
			Xpa=qdtr.Generate_conjaxis(self.dim)
			for i in range(newmeans.shape[0]):
				pp=np.multiply( Xpa,np.tile(newmeans[i],(Xpa.shape[0],1)) )
				if Xnewmeans is None:
					Xnewmeans=pp
				else:
					Xnewmeans=np.vstack((Xnewmeans,pp))

		Xnewmeans=np.vstack((oldmeans,Xnewmeans))

		d=[]
		for i in range(Xnewmeans.shape[0]):
			for j in range(i+1,Xnewmeans.shape[0]):
				d.append((i,j,np.linalg.norm( Xnewmeans[j]-Xnewmeans[i] ),np.linalg.norm(Xnewmeans[i]),np.linalg.norm(Xnewmeans[j]) ))
		d=filter(lambda x: x[2]<=thres,d)
		delmeans=set([x[0] if x[3]>x[4] else x[1]  for x in d])
		Xnewmeans=np.array([Xnewmeans[i] for i in range(Xnewmeans.shape[0]) if i not in delmeans])

		return self.Getgmms_bysamples(1000,Xnewmeans)


	def Getgmms_bysamples(self,N,means):
		Ncomp=means.shape[0]
		while True:
			X0=estmn.mvnrnd0I(self.dim,N=N)
			starttime=time.time()
			Xids=estmn.getclusterIDs(X0,means,ClusterIds=None)
			ClusterM=[]
			ClusterP=[]
			flg=0
			for i in range(Ncomp):
				XX=X0[Xids==i,:]
				if XX.shape[0]<100:
					N=N*2
					flg=1
					break

				W=estmn.mvnpdf(XX,np.zeros(self.dim),np.identity(self.dim))
				W=W/np.sum(W)
				m,p=estmn.GetMeanCov(XX,w=None)
				ClusterM.append(m)
				ClusterP.append(p)

			if flg==0:
				break

		# now optimize the weights
		b=estmn.mvnpdf(X0,np.zeros(self.dim),np.identity(self.dim))
		A=None
		for i in range(Ncomp):
			pp=estmn.mvnpdf(X0,ClusterM[i],ClusterP[i]).reshape(-1,1)
			if i==0:
				A=pp
			else:
				A=np.hstack((A,pp))



		results=sc.optimize.lsq_linear(A,b,bounds= (np.zeros(Ncomp), np.ones(Ncomp)) )
		weights=results.x
		weights=weights/np.sum(weights)
		ClusterW=weights


		return (ClusterW,ClusterM,ClusterP		)