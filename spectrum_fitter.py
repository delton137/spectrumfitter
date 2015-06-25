from pylab import *
from scipy import optimize 


class Debye:
	"""Debye lineshape object.
	
	Note: the tau parameter is assumed to have units of ps and frequencies are assumed to be in cm^-1
	"""

	def __init__(self,params=[1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "tau"]
		self.type = "Debye"
		 
	def __call__(self, w):
		convf = (2*3.14159)/33.35641 #ps -> cm with a 2 pi
		return self.p[0]*w*convf*self.p[1]/(1 + w**2*(convf*self.p[1])**2) 
       
class DHO:
	"""Damped harmonic oscillator object"""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wT","gamma"]
		self.type = "Damped Harmonic Oscillator"

		 
	def __call__(self, w):
		return self.p[0]*self.p[1]**2*self.p[2]*w/( (self.p[1]**2 - w**2)**2 + w**2*self.p[2]**2 ) 
	
class BRO:
	"""Briet-Rabbi oscillator"""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wT","gamma"]
		self.type = "Breit-Rabbi"
		 
	def __call__(self, w):
		return .5*self.p[0]*self.p[1]**2*self.p[2]*w*( 1/((self.p[1]**2 - w**2)**2 + w**2*self.p[2]**2)  + 1.00/((self.p[1]**2 + w**2)**2 + self.p[2]**2) ) 
	
class constant:
	"""this "lineshape" object is merely a constant term"""
	def __init__(self,params=[1],bounds=[(0,inf)],name="Eps Infinity"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["Eps Infinity"]
		self.type = "Constant"
	
	def __call__(self,w):
		return 0*w + self.p[0]

class spectralmodel: 
	"""A spectralmodel object is simply a list of lineshape objects"""
	def __init__(self,lineshapes=[]):
		self.lineshapes = lineshapes
		self.numlineshapes = 0
		self.RMS_error = None 
		
	def add(self,lineshape):
		"""add a new lineshape object to the spectral model's list of lineshapes"""
		self.lineshapes = self.lineshapes + [lineshape]
		self.numlineshapes += 1
		
	def setparams(self,params):
		"""set parameters in the model from a list of parameters for all lineshapes"""
		i = 0
		for l in self.lineshapes:
			for j in range(len(l.p)):
				l.p[j] = params[i]
				i += 1
	
	def getparams(self):
		"""get parameters and bounds for all the lineshapes in a model and return both as lists"""
		params = []
		bounds = []
		for lineshape in self.lineshapes:
			bounds = bounds + lineshape.bounds
			params = params + lineshape.p
		return params, bounds
				
	def eval(self,w):
		"""evaluate spectral_model model at frequencies in array w
			
		args: 
			w: an 1xN array with frequencies
		returns: 
			a 1xN array with the output of the model
		"""
		out = zeros(len(w))	
		for lineshape in self.lineshapes:
			out = out + lineshape(w)
		return out
	
	def print_model(self):
		"""Write out info about all of the parameters in the model """
		print "RMS error = ", self.RMS_error
		for lineshape in self.lineshapes:
			print lineshape.name, lineshape.p
					
	def fit_model(self, dataX, dataY):
		dataX = dataX
		dataY = dataY
		def diffsq(params):
			"""Wrapper function for differential_evolution()
	

			Args: 
				params: a list of parameters for the model
			Returns: 
				The sum of squared error (squared residuals) 
			"""
			self.setparams(params)
			diff = dataY - self.eval(dataX)
			return dot(diff,diff)

		params, b = self.getparams()

		#optimize.fmin_tnc(diffsq, params, fprime=None,approx_grad=True,args=params,bounds=b,epsilon=1e-08,)
		optimize.differential_evolution(diffsq,b) #perform optimization 
		
		self.RMS_error = sqrt(diffsq(params)/len(dataX)) #Store RMS error

