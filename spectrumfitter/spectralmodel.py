from scipy import optimize 
from numpy import *
import time

class SpectralModel: 
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
            print(i)
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
        print( "")
        print( "f-sum of oscillator strengths = %6.2f" % self.fsum())
        print( "")
        print( "RMS error = %6.3f" % self.RMS_error)
        
    def print_model_latex(self):
        """Write out info about all of the parameters in the model in LaTeX form"""
        print("\\begin{table}")
        print("     \\begin{tabular}{c c c c c c}")
        print("name            & $f$ & $\\omega_0$ (cm$^{-1}$) & $\\tau$ (ps) &  $\gamma$ (cm$^{-1}) $  &  $\sigma$ (cm$^{-1}$) \\\\ ")
        print("\hline") 
        for lineshape in self.lineshapes:
            lineshape.print_params_latex()
        print("     \\end{tabular}}")
        print("\\end{table}")
    
    def fit_model(self, dataX, datarp, datacp, differential_evolution=True, TNC=True, SLSQP=True, verbose=True):
        '''Fit the function using one or multiple optimization methods in serial''' 
    
        def diffsq(params):
            self.setparams(params)
            
            (rp,cp) = self(dataX)
            
            diffrp = (datarp - rp)/datarp
            diffcp = (datacp - cp)/datacp 
            
            #Ldatacp = datacp/(datarp**2 + datacp**2)
            #Lfitcp =  cp/(rp**2 + cp**2)
            #diffLcp = (Ldatacp - Lfitcp)/Ldatacp
            
            return dot(diffcp, diffcp) + dot(diffrp, diffrp) #+ dot(diffLcp,diffLcp)

        def costfun(params):
            """Wrapper function neede for the optimization method

            Args: 
                params: a list of parameters for the model
            Returns: 
                The cost (real scalar)
            """
            Error = diffsq(params) 
            
            fsumpenalty = datarp[0] - self.fsum()
            
            return Error + fsumpenalty**2 

        start_t = time.time()

        params = self.getparams()
        bounds = self.getbounds()

        if (differential_evolution == True):
            resultobject = optimize.differential_evolution(costfun,bounds,maxiter=2000)  
            if (verbose == True): print("diff. evolv. number of iterations = ", resultobject.nit)

        if (TNC == True):
            resultobject = optimize.minimize(costfun, x0=params, bounds=bounds, method='TNC')
            if (verbose == True): print("TNC number of iterations = ", resultobject.nit)
        
        if (SLSQP == True):
            resultobject = optimize.minimize(costfun, x0=params, bounds=bounds, method='SLSQP')
            if (verbose == True): print("SLSQP number of iterations = ", resultobject.nit)

        #mybounds = MyBounds(bounds=array(bounds))
        #ret = basinhopping(diffsq, params, niter=10,accept_test=mybounds)

        params = self.getparams() #get updated params
        
        end_t = time.time()
        m, s = divmod(end_t - start_t, 60)
        h, m = divmod(m, 60)
        print("Fit completed in %02d hr %02d min %02d sec" % (h, m, s))

        self.RMS_error = sqrt(diffsq(params)/(2*len(dataX))) #Store RMS error