from pylab import *
from scipy import optimize 

#class MyBounds(object):
	#def __init__(self, bounds=array[[-10,10]] ):
		#self.xmax = bounds[:,1]
		#self.xmin = bounds[:,2]
	#def __call__(self, **kwargs):
		#x = kwargs["x_new"]
		#tmax = bool(all(x <= self.xmax))
		#tmin = bool(all(x >= self.xmin))
		#return tmax and tmin

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
		rp = self.p[0]*self.p[1]**2/(self.p[1]**2 + w**2) 	
		cp = rp*w/self.p[1]
		return (rp, cp)
	
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
		self.pnames = ["f", "wT", "gamma"]
		self.type = "DHO"

	def __call__(self, w):
		rp = self.p[0]*(self.p[1]**2)*(self.p[1]**2 - w**2)/( (self.p[1]**2 - w**2)**2 + w**2*self.p[2]**2 ) 
		cp = self.p[0]*(self.p[1]**2)*self.p[2]*w/( (self.p[1]**2 - w**2)**2 + w**2*self.p[2]**2 ) 
		return (rp, cp)
	
	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 ) 
	
class StretchedExp:
	"""Stretched Exponential lineshape """
	def __init__(self,params=[1,1,1],bounds=[(0,10000),(0,10000),(0,1)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "tau", "beta"]
		self.type = "StretchedExp"

	def __call__(self, w):
		(f, eps_omega) = self.calc_eps(w)
		rp = interp(w,f,real(eps_omega))
		cp = interp(w,f,imag(eps_omega))
		return (rp, cp) 
	
	def calc_eps(self,w):
		deltaeps = self.p[0]
		tau = self.p[1]
		beta = self.p[2]
		timestep = 1/(2*3.141*max(w)) #timestep in ps 
		maxt = 10*tau #length of corr function
		npts = int(maxt/timestep)
		times = array(range(0,npts))*timestep
		corr_fun = exp(-(times/tau)**beta)
		nextpow2 = int(2**ceil(log2(npts)))
		f = (1/(2.99*.01))*array(range(0,npts))/(timestep*nextpow2) #Frequencies (evenly spaced), inverse cm
		tderiv = -diff(corr_fun)/timestep #time derivative of corr fun
		eps_omega = deltaeps*ifft(tderiv,nextpow2) #inverse one-sided fourier 
		return (f, eps_omega[0:npts])
		
	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 ) 
	
		tderiv = diff(corr_funs[:,k])/timestep
class PowerLawDebye:
	"""Debye lineshape with additional "power law" wing. 
	See J. Phys. Chem. B, 2005, 109 (12), pp 6031-6035 
	
	Note: the wD parameter is assumed to be in units of cm^-1
	"""
	def __init__(self,params=[1,1],bounds=[(0,+inf),(0,+inf),(0,+inf),(1,10)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wD","A","q"]
		self.type = "PowerLawDebye"
		self.convfac = 1.0#/(2*3.141*2.99*.01)
	
	def __call__(self, w):
		A = self.p[2]
		q = self.p[3]
		tau = self.convfac/self.p[1]
		numomegas = len(w)
		start  = numomegas - len(w[w[:] > self.p[1]])
		HighFreqOmegas = 0.0*array(w)
		HighFreqOmegas[start:numomegas] = w[start:numomegas]
		TheWing = 1 + A*(HighFreqOmegas*tau)**q
		rp = TheWing*self.p[0]/(1 + (tau*w)**2)	
		cp = TheWing*self.p[0]*w*tau/(1 + (tau*w)**2)
		return (rp, cp)
	
	def get_freq(self): 
		return self.p[1]
	
	def get_abs_freq(self):
		return self.p[1]


		
	
class BRO:
	"""Briet-Rabbi oscillator (also called van Vleck-Weisskopf lineshape)"""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wT","gamma"]
		self.type = "BRO"
		
	#the real part has not been tested or checked to see if it is mathematically correct!!!
	def __call__(self, w):
		rp = .5*self.p[0]*(self.p[1]**2 + self.p[2]**2)( 1/((self.p[1] - w)**2 + self.p[2]**2)  + 1.00/((self.p[1] + w)**2 +  self.p[2]**2) ) 		 
		cp = .5*self.p[0]*self.p[2]*w*( 1/((self.p[1] - w)**2 + self.p[2]**2)  + 1.00/((self.p[1] + w)**2 +  self.p[2]**2) ) 
		return (rp, cp)

	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 ) 
	
class constant:
	"""this "lineshape" object is merely a constant term"""
	def __init__(self,params=[1],bounds=[(0,10000)],name="Eps Infinity"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["Eps inf."]
		self.type = "Constant"
	
	def __call__(self,w):
		rp = 0*w + self.p[0]
		cp = 0*w
		return (rp, cp)

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
	
	def getfreqs(self):
		"""get frequencies for all the lineshapes in a model and return as list"""
		freqs = zeros(self.numlineshapes)
		for i in range(self.numlineshapes):
			#if lineshape.type != "Constant": 
			print i
			freqs[i] = self.lineshapes[i].p[1]
		return freqs
	
	def fsum(self):
		"""evaluate the f-sum rule (sum the oscillator strengths)"""
		fsum = 0 
		for i in range(self.numlineshapes):
			fsum = fsum + self.lineshapes[i].p[0]
		return fsum
				
	def __call__(self,w):
		"""compute real and complex parts of the spectral_model model at frequencies in array w
			
		args: 
			w: an 1xN array with frequencies
		returns: 
			(rp, cp) a list with rp and cp as 1xN arrays  
		"""
		rp = zeros(len(w))	
		cp = zeros(len(w)) 
		for lineshape in self.lineshapes:
			(rpPart, cpPart) = lineshape(w)
			rp = rp + rpPart
			cp = cp + cpPart
		return (rp,cp)
	
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
					
	def fit_model(self, dataX, datarp, datacp):
	
		def diffsq(params):
			"""Wrapper function for differential_evolution()

			Args: 
				params: a list of parameters for the model
			Returns: 
				The sum of squared error (squared residuals) 
			"""
			self.setparams(params)
			(rp,cp) = self(dataX)
			diffrp = (datarp - rp)/datarp
			diffcp = (datacp - cp)/datacp 
			penalty = (78 - self.fsum())/78
			return dot(diffrp,diffrp) + dot(diffcp,diffcp) + penalty**2

		params = self.getparams()
		bounds = self.getbounds()
		#optimize.fmin_tnc(diffsq, params, fprime=None,approx_grad=True,args=params,bounds=b,epsilon=1e-08,)
		optimize.differential_evolution(diffsq,bounds) #perform optimization
		#mybounds = MyBounds(bounds=array(bounds))
				
		#ret = basinhopping(diffsq, params, niter=10,accept_test=mybounds)

		params = self.getparams() #get updated params
		#print params
		#params = ret.x
		#print params
		#self.setparams(params)

		self.RMS_error = sqrt(diffsq(params))/(2*len(dataX)) #Store RMS error  
		
def plot_model(model,dataX,dataYrp,dataYcp,Myhandle,plotmin,plotmax,scale='lin',show=False,Block=True):
	"""displays a pretty plot of the real and complex parts of the model and data using matplotlib

    args: 
        model: a spectral_model object
        dataX: a numpy array giving the experimental x-data
        dataYrp: a numpy array giving the real part of the experimental y-data
        dataYcp: a numpy array giving the complex part of the experimental y-data
        handle: an integer giving the plot window number 
        plotmin: scalar, minimum frequency to plot 
        plotmax: scalar, maximum frequency to plot
        show: logical, option to display plot
        blockoption: logical, option to block further processing after displaying window (Default True, False is experimental) 
        """
        if (scale == 'log'):
		plotomegas = logspace(log10(plotmin), log10(plotmax), 10000)
	else: 
		plotomegas = linspace(plotmin, plotmax, 10000)
		
	(rp, cp) = model(plotomegas)
	
	# Two subplots, unpack the axes array immediately
	f, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, sharey=False )
	
	ax1.plot(dataX, dataYrp, "ro", plotomegas, rp,'g')
	ax1.set_title('Real part')


	ax2.plot(dataX, dataYcp, "ro", plotomegas, cp,'g')
	ax2.set_title('Complex part')
	#plot all of the components
	for lineshape in model.lineshapes: 
		(rpPart, cpPart) = lineshape(plotomegas)
		ax1.plot(plotomegas, rpPart ,'g--')
		ax2.plot(plotomegas, cpPart ,'g--')

        
	if (scale == 'log'):
		ax1.set_xscale('log')
		ax2.set_xscale('log')
		xlim([plotmin,plotmax])

	ax1.set_yscale('log')
	ax1.set_ylim([.1*min(dataYrp),max(dataYrp)+1])
	ax2.set_yscale('log')
	ax2.set_ylim([min(dataYcp),max(dataYcp)+1])
	
	xlabel(r"$\omega$ cm$^{-1}$")
	ylabel(r"$\varepsilon(\omega)''$")

	LdataY = dataYcp/(dataYrp**2 + dataYcp**2)
	LfitY =  cp/(rp**2 + cp**2) 

	fig2 = plt.figure()
	ax3 = fig2.add_subplot(111)
	ax3.plot(dataX,LdataY,'ro',plotomegas,LfitY,'g')

 
    #ax.annotate('local max', xy=(3, 1),  xycoords='data')
	
   

