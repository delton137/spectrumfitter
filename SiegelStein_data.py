import numpy as np
import matplotlib.pyplot as plt

cspeed = 3*10**10

data = np.loadtxt(fname='Siegelstein.RI')
omegas = 10000/data[:,0] #convert um into cm^-1
n = data[:,1]
k = data[:,2]
rp = n**2 + k**2
cp = 2*n*k
T = rp
L = cp#/(rp**2 + cp**2)

 
plt.loglog(omegas,T,'r',label=r"$\chi_T''(\omega)$" )
plt.loglog(omegas,L,'g',label=r"$\chi_L''(\omega)$" )
#plt.loglog(omegas,T,'r',label=r"$\varepsilon'(\omega)$"  )
#plt.loglog(omegas,L,'g',label=r"$\varepsilon''(\omega)$" )
#plt.loglog(omegas,T,'r',label=r"$n(\omega)$"  )
#plt.loglog(omegas,L,'g',label=r"$k(\omega)$" )
plt.title("Siegelstein Data",fontsize=17)
plt.xlabel(r"$\omega$ (cm$^{-1}$)",fontsize=24)
plt.legend(fontsize='x-large',frameon=False, loc='upper right')
f = plt.gcf()
plt.setp([a.get_yticklabels() for a in f.axes], fontsize=19)
plt.setp([a.get_xticklabels() for a in f.axes], fontsize=19)

 

plt.show() 
   