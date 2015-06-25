from pylab import *
from scipy import optimize 

class Debye:
	"""Debye lineshape object.
	
	Note: the wD parameter is assumed to be in units of cm^-1
	"""

	def __init__(self,params=[1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wD"]
		self.type = "Debye"
		 
	def __call__(self, w):
		#convf = (2*3.14159)/33.4 #ps -> cm with a 2 pi
		return self.p[0]*w*self.p[1]/(self.p[1]**2 + w**2) 
	
	def get_freq(self): 
		return self.p[1]
	
	def get_abs_freq(self):
		return self.p[1]
       
class DHO:
	"""Damped harmonic oscillator object"""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wT","gamma"]
		self.type = "DHO"

		 
	def __call__(self, w):
		return self.p[0]*self.p[1]**2*self.p[2]*w/( (self.p[1]**2 - w**2)**2 + w**2*self.p[2]**2 ) 
	
	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 ) 
	
class BRO:
	"""Briet-Rabbi oscillator (also called van Vleck-Weisskopf lineshape)"""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wT","gamma"]
		self.type = "BRO"
		 
	def __call__(self, w):
		return .5*self.p[0]*self.p[2]*w*( 1/((self.p[1] - w)**2 + self.p[2]**2)  + 1.00/((self.p[1] + w)**2 +  self.p[2]**2) ) 
	
	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 ) 
	
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
		"""get parameters for all the lineshapes in a model and return as list"""
		params = []
		for lineshape in self.lineshapes:
			params = params + lineshape.p
		return params
	
	def getbounds(self):
		"""get bounds for all the lineshapes in a model and return as list"""
		bounds = []
		for lineshape in self.lineshapes:
			bounds = bounds + lineshape.bounds
		return bounds
				
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
			
	def list(self):
		"""Write out info about all of the parameters in the model """
		print "RMS error = ", self.RMS_error
		for lineshape in self.lineshapes:
			print lineshape.name, lineshape.p
					
	def fit_model(self, dataX, dataY):

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

		params = self.getparams()
		bounds = self.getbounds()
		#optimize.fmin_tnc(diffsq, params, fprime=None,approx_grad=True,args=params,bounds=b,epsilon=1e-08,)
		optimize.differential_evolution(diffsq,bounds) #perform optimization 

		params = self.getparams() #get updated params
		self.RMS_error = sqrt(diffsq(params)/len(dataX)) #Store RMS error  
		
		
def plot_model(model,dataX,dataY,handle,plotmin,plotmax):
    """displays a pretty plot of the model and data using matplotlib
    
    args: 
        model: a spectral_model object
        dataX: a numpy array giving the experimental x-data
        dataY: a numpy array giving the experimental y-data
        handle: an integer giving the plot window number 
        plotmin: scalar, minimum frequency to plot 
        plotmax: scalar, maximum frequency to plot
    """
    figure(handle)

    plotomegas = linspace(10*plotmin, plotmax, 10000)

    plot(dataX, dataY, "ro", plotomegas, model.eval(plotomegas),'g')

    #plot all of the components
    for lineshape in model.lineshapes: 
        plot(plotomegas, lineshape(plotomegas),'g--')
    #ax.annotate('local max', xy=(3, 1),  xycoords='data')

    #xlim([plotmin,plotmax])
    #ylim([.001,max(dataY)])

