#Include numerical python methods, keeping the 'np' namespace 
#to avoid conflicts in names
import numpy as np
from np_utils import find

#Include plot tools (e.g. pylab.plot(x, y, .....))
#Again, we keep the namespace (if not we can't make a function called e.g. 'plot')
import matplotlib as mpl
from matplotlib import pylab

#Standard python modules for reading commandlines and so on..
import sys, os

from math import pi, exp, sqrt

#Python programmers use to underscore_separator convention usually
#unlike java / c++ programmers which use a newCapitalLetter convention
class kMC:

    #Declaring the system variables (not necessary, but transparent)
    x = None
    z = None
    
    E = None
    R = None
    accu_all_rates = None
    
    E0_p = 1
    E0_v = 1
    E0_s = -1
    
    alpha_p = 1
    alpha_v = 1      
    r_0 = None #length constant stress
    
    system_size = None    
    inert_surface_loc = None
    beta = None
    
    
    def __init__(self, size, T, init_height, r_0):
        
        self.system_size = size
        self.beta = 1./T
        self.r_0 = r_0
    
        
        self.x = np.arange(size) #0, 1, 2, 3, 4 ... 
        self.z = np.zeros(size, dtype=int)
        self.E = np.zeros(size)
        self.R = np.zeros((size, 2))
        self.accu_all_rates = np.zeros(2*size)
        
        
        self.init_surface(init_height)
        
        
        self.E_vertical_saddle_single = self.E0_p/((sqrt(5)/2)**self.alpha_p)                
        
        
        self.get_all_energies();
        self.get_all_rates();
    
    
    def init_surface(self, init_height):
        
        self.z += init_height/2
        
#        z = 0.25*np.sin(2*pi*self.x/(self.system_size/3)) + 0.5
#        z[-1] = self.z[0] #Force initial periodicity
#
#        dz = z.max() - z.min()
#        for i, z_i in enumerate(z):
#            self.z[i] = int((z.min() + z_i/dz)*init_height)
            
        self.inert_surface_loc = self.z.max()*1.5
        
    
    def save_data(self):
        
        np.save("/tmp/1DKMC.npy", np.vstack((self.z, 
                                             self.E, 
                                             self.R[:, 0], 
                                             self.R[:, 1], 
                                             self.accu_all_rates[:self.system_size],
                                             self.accu_all_rates[self.system_size:])).transpose())
        
        
    def get_all_energies(self):
        
        for x_i in self.x:
            self.E[x_i] = self.get_energy(x_i)
          
          
    def get_all_rates(self):
        
        for x_i in self.x:
            self.R[x_i][0] = self.get_rate(x_i, -1)
            self.R[x_i][1] = self.get_rate(x_i,  1)
           
           
    def get_rate(self, x_i, dx):

        z_i = self.z[x_i]          

        #No particle present is given a zero rate.
        if z_i == 0:
            return 0
        
        x_r = (x_i + dx)%self.system_size
        x_l = (x_i - dx)%self.system_size
      
        z_r = self.z[x_r]
        z_l = self.z[x_l]
        
        
        #Clogging is also zero rate.        
        if (z_r >= self.inert_surface_loc - 1):
            return 0
        
        dz_r = z_r - z_i
        dz_l = z_l - z_i
        
        vertical_path = dz_r + 1
        
        R = 1 #rate of sliding on surface is always one.
        
        E_i = self.E[x_i]
        
        if vertical_path != 0:
        
            if (dz_r > 0):
                
                #Move upwards dz_r times.
                for i in xrange(dz_r):
                        
                    E_saddle = self.eval_vertical_saddlepoint(dz_r - i, dz_l - i)      
        
                    R *= self.get_specific_rate(E_i, E_saddle)
                    
                    self.z[x_i] = z_i + i + 1
                    E_i = self.get_energy(x_i) - self.E0_p #subtract the bottom connector.
                
                #Then we slide across the top.
                R *= self.get_specific_rate(E_i, self.E_vertical_saddle_single) #Rate of sliding over kink
            
                self.z[x_i] = z_i
                
                
            else:
                
                #First we slide across the kink
                R *= self.get_specific_rate(E_i, self.E_vertical_saddle_single) #Rate of sliding over kink
                
                self.z[x_i] = z_i - 1      
                
                x_rr = (x_i + 2*dx)%self.system_size
                z_rr = self.z[x_rr]
                
                dz_rr = z_rr - z_i
                
                #Move downwards -dz_r + 1 times
                for i in xrange(-dz_r + 1):
                    
                    self.z[x_r] = z_i - i
                    E_i = self.get_energy(x_r)
                    
                    #Since we moved the particle to the right, dz_rr is shortened by dz_r.
                    #The first time the left difference is -1, then each height increase by 1 each cycle.
                    E_saddle = self.eval_vertical_saddlepoint(dz_rr + dz_r + i, i - 1)
                
                    R *= self.get_specific_rate(E_i, E_saddle)                    
                    
                self.z[x_i] = z_i
                self.z[x_r] = z_r
                
        return R
       
       
    def eval_vertical_saddlepoint(self, dz_r, dz_l):

        E_saddle = 0
        
        if (dz_r >= 1):
            E_saddle += 2*self.E_vertical_saddle_single
        elif (dz_r == 0):
            E_saddle += self.E_vertical_saddle_single
            
            
        if (dz_l >= 1):
            E_saddle += 2*self.E_vertical_saddle_single
        elif (dz_l == 0):
            E_saddle += self.E_vertical_saddle_single
        
        return E_saddle
               
               
    def get_specific_rate(self, E, E_saddle):
        
        return exp(-self.beta*(E - E_saddle))
    
    
    def get_energy(self, x_i):
            
        E_i = 0
        
        z_i = self.z[x_i]
        
        #Particle-particle electro-static contribution
        if z_i != 0:
            
            E_i += self.E0_p
            
            z_im = self.z[(x_i - 1)%self.system_size]
            z_ip = self.z[(x_i + 1)%self.system_size]
            
            if (z_im >= z_i):
                E_i += self.E0_p
            if (z_ip >= z_i):
                E_i += self.E0_p
            
        
        #Particle-inert surface electro-static contribution                        
                
        E0_pv = sqrt(self.E0_p*self.E0_v)
        h = self.inert_surface_loc - z_i
        
        E_i += E0_pv/h
                
                
        #Particle mechanical stress energy contribution      
        E_i += self.E0_s*exp(-h/self.r_0)
                
        return E_i

    def move_particle(self, x_o, x_d):
        
        if (self.z[x_d] > self.inert_surface_loc - 1):
            print "warning: move requested above surface. R = %g, %g, from x = %d to x = %d" % (self.R[x_o][0], self.R[x_o][1], x_o,  x_d)
            sys.exit()
        if (self.z[x_o] == 0):
            print "warning: move requested below simulation zero. R = %g, %g, from x = %d to x = %d" % (self.R[x_o][0], self.R[x_o][1], x_o,  x_d)
            sys.exit()
            
        self.z[x_o] -= 1
        self.z[x_d] += 1
        
        affected = set([(x_o - 1)%self.system_size,
                         x_o,
                        (x_o + 1)%self.system_size,
                        (x_d - 1)%self.system_size,
                        x_d,
                        (x_d + 1)%self.system_size])
                        
        for x_i in affected:
            self.E[x_i]     = self.get_energy(x_i)
            self.R[x_i][0]  = self.get_rate(x_i, -1)
            self.R[x_i][1]  = self.get_rate(x_i,  1) 
        
    
    def run(self, n_c):

        cycle = 1
        
        while cycle <= n_c:
            
            self.accu_all_rates = np.cumsum(self.R.flatten())
            
            R_tot = self.accu_all_rates[-1]
            
            dart_hit = R_tot*np.random.uniform()

            R_tot = 0  
            choice = None
            
            result = find(self.accu_all_rates, lambda x : x >= dart_hit)
            (choice,), val = next(result)

            x_i = choice/2
            dx  = 2*(choice%2) - 1
            
            self.move_particle(x_i, (x_i + dx)%self.system_size)
            
            if (cycle%1000) == 0:
                self.save_data()
                print "cycle", cycle, "/", n_c
            
            cycle += 1
        




#An if-test that executes the kMC function if and only if it is used 
#through this scope (to avoid it lauching if included in another file)
if __name__=="__main__":
    
    system_size = 1000
    temperature = 1.0
    init_height = 10
    r_0         = 1 
    n_c         = 100000

    solver = kMC(system_size, temperature, init_height, r_0)
    solver.run(n_c)

        
       
