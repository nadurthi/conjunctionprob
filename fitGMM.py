import numpy as np
import scipy as sc


def getrawmoms(N,mu,sig):
	X=[]
	for i in range(N):
		C=sc.special.hermitenorm(i,True)
		

class GMM0I(object):
	def __init__(self,dim=2):
		self.dim=dim

		pass

	def Solvegmms(means):
		# first get the positive quadrant points
		for i in range(self.dim): 	
			means=means[ means[i]>=0,:]
