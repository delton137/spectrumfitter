''' Spectrum_fitter.py : an obect-oriented framework for fitting dielectric spectra. '''

__author__  = "Daniel C. Elton"
__maintainer__ = "Daniel C. Elton"
__copyright__ = "Copyright 2015, Daniel C. Elton"
__license__ = "MIT"
__status__ = "Development"

from pylab import *
from scipy import optimize 
from scipy import special as sp

class Debye:
	"""Debye lineshape object.
	Note: the wD parameter is assumed to be in units of cm^-1
	"""
	def __init__(self,params=[1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="Debye"):
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
	
	def print_params(self):
		print u"%11s f=%4.2f \u03C9 %6.2f 1/cm (%5.2f ps)" % (self.name, self.p[0], self.p[1], 33.34/(2*3.141*self.p[1]) )
       
class DHO:
	"""Damped harmonic oscillator object"""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="DHO"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wT", "gamma"]
		self.type = "DHO"

	def __call__(self, w):
		denom = (self.p[1]**2 - w**2)**2 + w**2*self.p[2]**2
		rp = self.p[0]*(self.p[1]**2)*(self.p[1]**2 - w**2)/denom
		cp = self.p[0]*(self.p[1]**2)*self.p[2]*w/denom
		return (rp, cp)
	
	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 ) 
	
	def print_params(self):
		print u"%11s f= %4.2f \u03C9= %6.2f + %6.2f i  %6.3f" % (self.name, self.p[0], self.p[1], self.p[2],  33.34/(2*3.141*self.p[2]))	

	
class BrendelDHO:
	""""Brendel model for amorphous materials - Gaussian distribution of DHOs. Ref: J. Appl. Phys. 71, 1 (1992)"""
	def __init__(self,params=[1,1,1,1],bounds=[(0,+inf),(0,+inf),(0,+inf),(0,+inf)],name="BrendelDHO"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["wp**2/w0**2", "wT", "gamma","sigma"]
		self.type = "BrendelDHO"
		self.f = 0 

	def __call__(self, w):
		sigma = self.p[3]
		x0 = self.p[1]
		g = self.p[2]
		a = sqrt(w**2 - 1j*g*w) 
		a = a.real - 1j*a.imag #we want the imaginary part to the root to be positive
		prefac = 1j*sqrt(3.14149)*self.p[0]*x0**2/(sqrt(22)*sigma)
		eps = prefac*exp(-.5)*(1/a)*( sp.erfcx(-1j*(a-x0)/sigma) +  sp.erfcx(-1j*(a+x0)/sigma) )
		#self.f = 2*prefac*exp(-x0**2/(2*sigma**2))*(1 + sp.erf(1j*x0/sigma))
		self.f = eps.real[1]
		return (eps.real, eps.imag)
	
	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 )
	
	def print_params(self):
		print u"%11s f =%4.2f \u03C9= %6.2f + %6.2f i  (%5.3f ps) \u03C3 = %6.2f cm^-1" % (self.name, self.p[0], self.p[1], self.p[2],  33.34/(2*3.141*self.p[2]), self.p[3])	
	
class StretchedExp:
	"""Stretched Exponential lineshape """
	def __init__(self,params=[1,1,1],bounds=[(0,10000),(0,10000),(0,1)],name="Str exp"):
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
		maxt = 3*tau #length of corr function
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
	
	def print_params(self):
		print u"%11s f =%4.2f \u03C9 %6.2f 1/cm (%5.2f ps) \u03B2 = %6.2f" % (self.name, self.p[0], self.p[1], 33.34/(2*3.141*self.p[1]), self.p[2] )
	
		
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
		HighFreqOmegas = w#0.0*array(w)
		#HighFreqOmegas[start:numomegas] = w[start:numomegas]
		TheWing = 1 + A*(HighFreqOmegas*tau)**q
		rp = TheWing*self.p[0]/(1 + (tau*w)**2)	
		cp = TheWing*self.p[0]*w*tau/(1 + (tau*w)**2)
		return (rp, cp)
	
	def get_freq(self): 
		return self.p[1]
	
	def get_abs_freq(self):
		return self.p[1]
	
	def print_params(self):
		print u"%11s f =%4.2f \u03C9 %6.2f 1/cm (%5.2f ps) A = %6.2f q = %6.2f" % (self.name, self.p[0], self.p[1], 33.34/(2*3.141*self.p[1]),self.p[2], self.p[3] )
	
class ColeCole:
	"""Cole-Cole lineshape object.
	
	Note: the wD parameter is assumed to be in units of cm^-1
	"""
	def __init__(self,params=[1,1,1],bounds=[(0,+inf),(0,+inf),(0,2)],name="unnamed"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wD","alpha"]
		self.type = "Debye"
	
	def __call__(self, w):
		rp = self.p[0]*self.p[1]**2/(self.p[1]**2 + w**2) 	
		cp = rp*w/self.p[1]
		return (rp, cp)
	
	def get_freq(self): 
		return self.p[1]
	
	def get_abs_freq(self):
		return self.p[1]
	
	def print_params(self):
		print u"%11s f =%4.2f \u03C9 %6.2f 1/cm (%5.2f ps) \u03B1 %6.2f" % (self.name, self.p[0], self.p[1], 33.34/(2*3.141*self.p[1]),self.p[2] )

class VanVleck:
	"""Van Vleck and Weisskopf lineshape Rev. Mod. Phys., 17:227 236, Apr 1945."""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="BRO"):
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

	def print_params(self):
		print u"%11s f =%4.2f \u03C9= %6.2f + %6.2f i  %6.3f" % (self.name, self.p[0], self.p[1], self.p[2],  33.34/(2*3.141*self.p[2]))	


class Gaussian:
	"""Gaussian peak in the imaginary part, ie. inertial absorption or homogenous broadening"""
	def __init__(self,params=[1,1,1],bounds=[(-inf,+inf),(-inf,+inf)],name="Gaussian"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["f", "wT","gamma"]
		self.type = "BRO"
		
	def __call__(self, w):
		# sp.dawsn
		rp = .5*self.p[0]*(self.p[1]**2 + self.p[2]**2)( 1/((self.p[1] - w)**2 + self.p[2]**2)  + 1.00/((self.p[1] + w)**2 +  self.p[2]**2) ) 		 
		cp = .5*self.p[0]*self.p[2]*w*( 1/((self.p[1] - w)**2 + self.p[2]**2)  + 1.00/((self.p[1] + w)**2 +  self.p[2]**2) ) 
		return (rp, cp)

	def get_freq(self):
		return self.p[1]
	
	def get_abs_freq(self):
		return sqrt( self.p[1]**2 + self.p[2]**2 ) 

	def print_params(self):
		print u"%11s f =%4.2f \u03C9= %6.2f + %6.2f i  %6.3f" % (self.name, self.p[0], self.p[1], self.p[2],  33.34/(2*3.141*self.p[2]))


	
class constant:
	"""this "lineshape" object is merely a constant term"""
	def __init__(self,params=[1],bounds=[(0,10000)],name="Eps Inf"):
		self.p = params
		self.bounds = bounds
		self.name = name
		self.pnames = ["Eps inf."]
		self.type = "Constant"
	
	def __call__(self,w):
		rp = 0*w + self.p[0]
		cp = 0*w
		return (rp, cp)
	
	def print_params(self):
		print "%11s f =%4.2f" % (self.name, self.p[0])

class spectralmodel: 
	"""A spectralmodel object is simply a list of lineshape objects"""
	def __init__(self,lineshapes=[]):
		self.lineshapes = lineshapes
		self.numlineshapes = 0
		self.RMS_error = 0 
		
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
			if (self.lineshapes[i].p[0] == "BrendelDHO"): 
				fsum = fsum + self.lineshapes[i].f
			else:
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
	
	def longeps(self,w):
		''' computes the longitudinal dielectric function for the spectral_model model at frequencies in array w'''
		(rp, cp) = self(w)
		denom = rp**2 + cp**2 
		return (rp/denom, cp/denom)
	
	def print_model(self):
		"""Write out info about all of the parameters in the model """
		for lineshape in self.lineshapes:
			lineshape.print_params()
		set_printoptions(precision=2)
		print ""
		print "Sum of trans. oscillator strengths = %6.2f" % self.fsum()
		print ""
		print "trans RMS error = %6.3f" % self.RMS_error
	
	def fit_model(self, dataX, datarp, datacp):
	
		def diffsq(params):
			self.setparams(params)
			(rp,cp) = self(dataX)
			diffrp = (datarp - rp)/datarp
			diffcp = (datacp - cp)/datacp 
			
			#Ldatacp = datacp/(datarp**2 + datacp**2)
			#Lfitcp =  cp/(rp**2 + cp**2)
			#diffLcp = (Ldatacp - Lfitcp)/Ldatacp
			return dot(diffrp,diffrp) + dot(diffcp,diffcp)
	
		def costfun(params):
			"""Wrapper function neede for differential_evolution()

			Args: 
				params: a list of parameters for the model
			Returns: 
				The cost function
			"""
			Error = diffsq(params) 
			penalty = (datarp[0] - self.fsum())/datarp[0] 
			
			return Error + 800*penalty**2 

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
		
def plot_model(model,dataX,dataYrp,dataYcp,Myhandle,xmin=None,xmax=None,xscale='lin',yscale='log',ymin=None,ymax=None,show=False,Block=True,longitudinal=False):
	"""displays a pretty plot of the real and complex parts of the model and data using matplotlib

    args: 
        model: a spectral_model object
        dataX: a numpy array giving the experimental x-data
        dataYrp: a numpy array giving the real part of the experimental y-data
        dataYcp: a numpy array giving the complex part of the experimental y-data
        handle: an integer giving the plot window number 
        xmin: scalar, minimum frequency to plot 
        xmax: scalar, maximum frequency to plot
        xscale: string, xscale type 'linear' or 'log'
        yscale: string, yscale type 'linear' or 'log' 
        show: logical, option to display plot
        blockoption: logical, option to block further processing after displaying window (Default True, False is experimental) 
        """
        if (xmin == None):
		xmin = min(dataX)
        if (xmax == None):
		xmax = max(dataX)
        if (ymin == None):
		ymin = min(dataYrp)
        if (ymax == None):
		ymax = max(dataYrp)
        
        if (xscale == 'log'):
		plotomegas = logspace(log10(xmin), log10(xmax), 10000)
	else: 
		plotomegas = linspace(xmin, xmax, 10000)
	
	if (longitudinal == True):
		(rp, cp) = model.longeps(plotomegas)
		denom = dataYrp**2 + dataYcp**2 	
		(dataYrp,dataYcp) = (dataYrp/denom, dataYcp/denom)
		(xmin,xmax) = (min(dataX),max(dataX))
		(ymin,ymax) = (min(dataYrp),max(dataYrp))
	else:
		(rp, cp) = model(plotomegas)
	
	# Two subplots, unpack the axes array immediately
	f, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, sharey=False )
	
	ax1.plot(dataX, dataYrp, "ro", plotomegas, rp,'g')
	ax1.set_title('Real part')


	ax2.plot(dataX, dataYcp, "ro", plotomegas, cp,'g')
	ax2.set_title('Complex part')
	#plot all of the components
	if longitudinal == False:
		for lineshape in model.lineshapes: 
			(rpPart, cpPart) = lineshape(plotomegas)
			ax1.plot(plotomegas, rpPart ,'g--')
			ax2.plot(plotomegas, cpPart ,'g--')

        
	ax1.set_xscale(xscale)
	ax1.set_xlim([xmin,xmax])
	ax1.set_yscale(yscale)
	ax1.set_ylim([ymin,ymax])

	ax2.set_xscale(xscale)
	ax2.set_xlim([xmin,xmax])
	ax2.set_yscale(yscale)
	ax2.set_ylim([ymin,ymax])
	
	xlabel(r"$\omega$ cm$^{-1}$")
	ylabel(r"$\varepsilon(\omega)''$")

   #ax.annotate('local max', xy=(3, 1),  xycoords='data')

   

