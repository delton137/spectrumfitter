from pylab import *
from scipy import optimize

class Parameter:
    def __init__(self, value, bounds=(-inf,+inf),name="unnamed"):
            self.value = value
            self.bounds = bounds
	    self.name = name

    def set(self, value):
            self.value = value

    def __call__(self):
            return self.value

#class model
    #def __init__(self,numDebye,numDHO)
	#self.numDebye = numDebye
	#self.numDHO = numDHO
	#self.lineshapes =   

    #def eval(w):
	#output = 0 
      	#for lineshape in lineshapes: 
	  #output = output + lineshape()
	    

def fit(function, parameters, x, y):
    
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return dot(y-function(x),y-function(x)) 

    if x is None: x = arange(y.shape[0])
    
    params = array([param() for param in parameters])
    b = [param.bounds for param in parameters]
        
    #optimize.fmin_tnc(f, p0, fprime=None,approx_grad=True,args=p0,bounds=b,epsilon=1e-08,)
    optimize.differential_evolution(f,b)
  
    
   
cspeed = 3*10**10
data = loadtxt(fname='Siegelstein.RI')

maxw = 1500
omegas = data[0:maxw,0] #assumed in cm^-1
n      = data[0:maxw,1]
k      = data[0:maxw,2]

rp = n**2 + k**2
cp = 2*n*k
  
convf = (2*3.14159)/33.35641 #ps -> cm with a 2 pi

dataX = omegas 
dataY = cp 

#----------------------------------------------------------------------------------
#---------------------- Define fit function & parameters -------------------------
#----------------------------------------------------------------------------------
fD1        = Parameter(71  ,(70,72), "fD1")
tauD1      = Parameter(8.2   ,(8,8.4) ,"tauD1") #in ps 

fD2        = Parameter(1  ,(.01,5), "fD2")
tauD2      = Parameter(1 ,(.01,2) ,"tauD2") #in ps 

f1        = Parameter(1.7  ,(1,5), "f1")
wT1       = Parameter(200 ,(100,300) , "wT1")  
gamma1    = Parameter(20   ,(0,400) , "gamma")  


f2        = Parameter(.5  ,(.01,2), "f2")
wT2       = Parameter(400  ,(300,500) , "wT2")  
gamma2    = Parameter(100   ,(1,400) , "gamma2")  


f4       = Parameter(.3  ,(.01,2), "f4")
wT4       = Parameter(500  ,(400,800) , "wT4")  
gamma4    = Parameter(10   ,(1,1000) , "gamma4")  


f3        = Parameter(1.4  ,(.01,5), "f3")
wT3       = Parameter(3500  ,(3000,4000) , "wT3")  
gamma3    = Parameter(100   ,(1,1000) , "gamma3")  



epsinfty   = Parameter(1.76,(0, 3),"epsinfty")

paramlist = [fD1, tauD1, fD2, tauD2,f1, wT1, gamma1, f2, wT2, gamma2, f3, wT3, gamma3, f4, wT4, gamma4, epsinfty]

#Debye function
def Debye(fD,tauD,w):
  return fD()*w*convf*tauD()/(1 + w**2*(convf*tauD())**2) 

#Damped harmonic oscillator
def DHO(f,wT,gamma,w):
  return f()*wT()**2*gamma()*w/( (wT()**2 - w**2)**2 + w**2*gamma()**2 ) 

#Breit-Rabbi lineshape

def DielectricFun(w):
  return Debye(fD1,tauD1,w) + Debye(fD2,tauD2,w) +  DHO(f1,wT1,gamma1,w) + DHO(f2,wT2,gamma2,w) + DHO(f3,wT3,gamma3,w) + DHO(f4,wT4,gamma4,w)  + epsinfty()

#---------------------- Execute fitting ------------------------------------------
fit(DielectricFun, paramlist, dataX, dataY)


#---------------------- Make pretty graph ----------------------------------------
plotmin = .1
plotmax = 4000
wfit = linspace(10*plotmin, plotmax, 10000)
loglog(dataX, dataY, "ro", wfit , DielectricFun(wfit), "r-", wfit, Debye(fD1,tauD1,wfit), "g--",wfit,Debye(fD2,tauD2,wfit), "g--",wfit,DHO(f1,wT1,gamma1,wfit), "g--",wfit,DHO(f2,wT2,gamma2,wfit), "g--",wfit, DHO(f3,wT3,gamma3,wfit), "g--",DHO(f4,wT4,gamma4,wfit), "g--") # Plot of the data and the fit

xlim([plotmin,plotmax])
ylim([.001,dataY.max()+2])


#---------------------- Write out info  ------------------------------------------
for param in paramlist:
  print "%s = %6.3f" %(param.name, param())

show()

