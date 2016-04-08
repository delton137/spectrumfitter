''' spectral_deconvolution.py : an obect-oriented framework for deconvolving dielectric spectral_deconvolution

'''
__author__  = "Daniel C. Elton"
__maintainer__ = "Daniel C. Elton"
__copyright__ = "Copyright 2016, Daniel C. Elton"
__license__ = "MIT"
__status__ = "Development"

from pylab import *
from scipy import optimize 
from scipy import special as sp

def deconvolve(self, dataX, datarp, datacp):
	
	def diffsq(params):
		self.setparams(params)
		(rp,cp) = self(dataX)
		diffrp = (datarp - rp)/datarp
		diffcp = (datacp - cp)/datacp 
		
		#Ldatacp = datacp/(datarp**2 + datacp**2)
		#Lfitcp =  cp/(rp**2 + cp**2)
		#diffLcp = (Ldatacp - Lfitcp)/Ldatacp
		return dot(diffcp,diffcp) + dot(diffrp,diffrp)  
	
		def costfun(params):
			"""Wrapper function neede for differential_evolution()

			Args: 
				params: a list of parameters for the model
			Returns: 
				The cost function
			"""
			
			
			
			
			Cost = diffsq(params) - 10*
			
			return Cost 

		params = self.getparams()
		bounds = self.getbounds()
		#optimize.fmin_tnc(diffsq, params, fprime=None,approx_grad=True,args=params,bounds=b,epsilon=1e-08,)
		optimize.differential_evolution(costfun,bounds) #perform optimization
		#mybounds = MyBounds(bounds=array(bounds))
				
		#ret = basinhopping(diffsq, params, niter=10,accept_test=mybounds)

		params = self.getparams() #get updated params
		#print params
		#params = ret.x
		#print params
		#self.setparams(params)

		self.RMS_error = sqrt(diffsq(params)/(2*len(dataX))) #Store RMS error  