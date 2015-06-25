from pylab import *
from scipy import optimize 
from numpy import *
from spectrum_fitter import *


def plot_model(model,dataX,dataY):
    """displays a pretty plot of the model and data using matplotlib
    
    args: 
        model: a spectral_model object
        dataX: a numpy array giving the experimental x-data
        dataY: a numpy array giving the experimental y-data
    """
    
    plotmin = .1
    plotmax = max_freq_to_analyze
    plotomegas = linspace(10*plotmin, plotmax, 10000)

    loglog(dataX, dataY, "ro", plotomegas, model.eval(plotomegas),'g')

    #plot all of the components
    for lineshape in model.lineshapes: 
        loglog(plotomegas, lineshape(plotomegas),'g--')
    #ax.annotate('local max', xy=(3, 1),  xycoords='data')

    #xlim([plotmin,plotmax])
    #ylim([.001,max(dataY)])

    show()

#----------------------------------------------------------------------------------
#---------------------- Load dielectric data -------------------------------------
#---------------------------------------------------------------------------------- 

data = loadtxt(fname='Siegelstein.RI')
#data = loadtxt(fname='eps_omega_flex_F_TTM3F.dat')

#cspeed = 3*10**10
#data[:,0] = data[:,0]/cspeed #conv data to cm^-1
max_freq_to_analyze = 4500
maxw   = len(data[data[:,0] < max_freq_to_analyze,0])
omegas = data[0:maxw,0]

n      = data[0:maxw,1]
k      = data[0:maxw,2]
rp = n**2 + k**2
cp = 2*n*k

rp      = data[0:maxw,1]
cp      = -data[0:maxw,2]
  
dataX = omegas 
dataY = cp/(rp**2 + cp**2)

#----------------------------------------------------------------------------------
#---------------------- Define fit function & parameters -------------------------
#----------------------------------------------------------------------------------

#longitudinal model
print "Fitting longitudinal model..."
modelL = spectralmodel()
modelL.add(Debye([.2, .2],[(.001,1),(.0001,.8)], "Debye"))
modelL.add(DHO([2,222,200]  ,[(0,5)  ,(10,290)  ,(0,400)  ],"H-bond str"))
modelL.add(BRO([.5,400,100]  ,[(.01,2),(300,800)  ,(.1,400) ],"L1"))
modelL.add(BRO([.3,500,200]  ,[(.01,2),(300,800)  ,(.1,1000)],"L2"))
modelL.add(DHO([1,1500,100]  ,[(.01,2),(1300,1700),(.1,500) ],"v2"))
modelL.add(DHO([1.4,3500,100],[(.01,2),(3000,4000),(.1,500) ],"v1+v3"))

#modelL.fit_model(dataX,dataY)

def diffsq(params):
	"""Wrapper function for differential_evolution()
		Args: 
			params: a list of parameters for the model
		Returns: 
			The sum of squared error (squared residuals) 
	"""
	modelL.setparams(params)
	diff = dataY - modelL.eval(dataX)
	return dot(diff,diff)

params, b = modelL.getparams()
#optimize.fmin_tnc(diffsq, params, fprime=None,approx_grad=True,args=params,bounds=b,epsilon=1e-08,)
optimize.differential_evolution(diffsq,b) #perform optimization 
	
modelL.RMS_error = sqrt(diffsq(params)/len(dataX)) #Store RMS error

plot_model(modelL,dataX,dataY) 
modelL.print_model()

###transverse model
#dataY = cp
#print "Fitting transverse model..."
#modelT = spectralmodel()
#modelT.add(Debye([71, 8.2],[(70,72),(8,8.4)], "Debye"))
#modelT.add(Debye([1,  .1],[(.01,2),(.01,4)],  "2nd Debye"))
#modelT.add(BRO([1.7,200,20]  ,[(1,5)  ,(100,300)  ,(0,400)  ],"H-bond str"))
#modelT.add(BRO([.5,505,300]  ,[(.01,2),(300,800)  ,(.1,400) ],"L1"))
#modelT.add(BRO([.3,600,244]  ,[(.01,2),(300,800)  ,(.1,1000)],"L2"))
#modelT.add(DHO([.3,500,200]  ,[(.01,2),(300,800)  ,(.1,1000)],"L3"))
#modelT.add(DHO([1,1500,100]  ,[(.01,2),(1300,1700),(.1,500) ],"v2"))
#modelT.add(DHO([1.4,3500,100],[(.01,2),(3000,4000),(.1,500) ],"v1+v3"))
  
#def diffsqT(params):
	#"""Wrapper function for differential_evolution()
		#Args: 
			#params: a list of parameters for the model
		#Returns: 
			#The sum of squared error (squared residuals) 
	#"""
	#modelT.setparams(params)
	#diff = dataY - modelT.eval(dataX)
	#return dot(diff,diff)

#params, b = modelT.getparams()
##optimize.fmin_tnc(diffsq, params, fprime=None,approx_grad=True,args=params,bounds=b,epsilon=1e-08,)
#optimize.differential_evolution(diffsq,b) #perform optimization 
	
#modelT.RMS_error = sqrt(diffsq(params)/len(dataX)) #Store RMS error

##plot_model(modelT)

##----------------------------------------------------------------------------------
##---------------------- Print left hand side of LST relation ---------------------
##----------------------------------------------------------------------------------

#params, b = modelL.getparams()
#print params
#params, b = modelT.getparams()
#print params



