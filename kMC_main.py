#Include numerical python methods, keeping the 'np' namespace 
#to avoid conflicts in names
import numpy as np

#Include plot tools (e.g. pylab.plot(x, y, .....))
#Again, we keep the namespace (if not we can't make a function called e.g. 'plot')
import matplotlib as mpl
from matplotlib import pylab

#Standard python modules for reading commandlines and so on..
import sys, os

from math import pi, exp

#Python programmers use to underscore_separator convention usually
#unlike java / c++ programmers which use a newCapitalLetter convention
class kMC:

    #Declaring the system variables (not necessary, but transparent)
    x = None
    y = None
    system_size = None    
    T = None
    
    C = 1      #prop const stress
    r_0 = None #length constant stress
    
    def __init__(self, size, T, r_0)
        self.system_size = size
        self.T = T
        self.r_0 = r_0
        
        z_init_func = lambda x: (np.sin(2*pi*x/3333) + 1)
        
        self.x = np.arange(size) #0, 1, 2, 3, 4 ... 
        z = z_init_func(x)
        
        z[-1] = z[0] #Force initial periodicity
        
        plot(z)
        
        self.z = z
    
        
    def stress(self, r):
        return C*exp(r/r_0)
        
    def chem_pot_diff(self, x):
        
        r_s = get_surface(x)
        
        dr = r_s - self.z[x] 
        
        a_s = 1
        a_0 = 1
        
        return self.T*log(a_s/a_0) - self.stress(dr)
        
    
    
















#An if-test that executes the kMC function if and only if it is used 
#through this scope (to avoid it lauching if included in another file)
if __name__=="__main__":
    kMC_main()       
