'''
Script to solve the spherical diffusion equation.
'''

import numpy as np
import matplotlib.pyplot as plt

class Diffusion():
    def __init__(self):
    # Put in constants here.
        self.R = 1e-05 # radius of particle
        self.C_max = 1200 # Cutoff conditions
        self.C_0 = 9500
        self.J_0 = 9.5e-6
        self.D_s = 1e-14 # Diffusion coefficient
        self.deltat = 1 # Time interval.
        self.deltar = (1/50) * self.R # radius_interval
        self.T_max = 2*3600 # End time
        self.t_array = np.linspace(0, self.T_max, deltat) # Time array
        self.r_array = np.linspace(0, self.R, deltar) # radius array
        self.C_array = np.ones((self.t_array.size,self.r_array.size))*self.C_0 # Concentration array, intialised.
#        self.j_array = np.zeros(self.T_max)
        
    def j(self,t_array):
    # Current as function of time. To be implemented later!
        self.j_array = np.sin(t_array)
            
    def boundary_conditions(self):
    # Script that defines the boundary conditions.
        self.C_array[:,0] = self.C_array[:,1]
        self.C_array[:,N] = self.C_array[:,N-1] - (self.deltar/self.D_s) * self.j(self.t_array)
        
    def solver(self):
    # This will solve the diffusion equation.
        for i in range(1,self.T_max):
            self.boundary_conditions()
            for k in range(0,self.R):
                self.C_array[i+1,k] = self.C_array[i,k] + self.deltat * self.D * (1 / self.r_array[k]) * \
                                      (self.C_array[i,k+1] - self.C_array[i,k-1]) / self.deltar) + \
                                      (self.C_array[i,k+1] - 2 * self.C_array[i,k] + self.C_array[i,k - 1]) / (self.deltar) ** 2 
                                                    

if __name__ == '__main__':
    class_init = Diffusion()
    solution  = class_init.solver()
    
