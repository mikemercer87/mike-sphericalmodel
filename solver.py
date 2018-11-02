'''
Script to solve the spherical diffusion equation.
'''

import numpy as np
import matplotlib.pyplot as plt

class Diffusion():
    def __init__(self):
    # Put in constants here.
        self.R = 1e-05 # radius of particle
        self.C_max = 12000 # Cutoff conditions
        self.C_0 = 9500
        self.J_0 = 9.5e-6
        self.D_s = 1e-14 # Diffusion coefficient
        self.deltat = 1 # Time interval.
        self.deltar = (1/50) * self.R # radius_interval
        self.no_rpoints = int(self.deltar / self.R)
        self.T_max = 2*3600 # End time
        self.no_tpoints = int(self.deltat / self.T_max)
        print(self.no_rpoints)
        print(self.no_tpoints)
        self.t_array = np.linspace(0, self.T_max, self.no_tpoints) # Time array
        self.r_array = np.linspace(0, self.R, self.no_rpoints) # radius array
        self.C_array = np.ones((self.t_array.size,self.r_array.size))*self.C_0 # Concentration array, intialised.
#        self.j_array = np.zeros(self.T_max)
        
    def j(self,t_array):
    # Current as function of time. To be implemented later!
        self.j_array = np.sin(0.5*np.pi*t_array/3600)
            
    def boundary_conditions(self):
    # Script that defines the boundary conditions.
        self.C_array[:,0] = self.C_array[:,1]
        self.j(self.t_array)
        self.C_array[:,-1] = self.C_array[:,-2] - self.J_0 * (self.deltar/self.D_s) * self.j_array
        
    def solver(self):
    # This will solve the diffusion equation.
        for i in range(1,self.T_max):
            self.boundary_conditions()
            for k in range(0,self.no_rpoints):
                self.numerator = self.C_array[i,k+1] - 2 * self.C_array[i,k] + self.C_array[i,k - 1]
                self.C_array[i+1,k] = self.C_array[i,k] + self.deltat * self.D * ((1 / self.r_array[k]) * (self.C_array[i,k+1] - self.C_array[i,k-1]) / self.deltar + self.numerator / (self.deltar) ** 2) 
                                                    

if __name__ == '__main__':
    diff_obj = Diffusion()
    solution  = diff_obj.solver()
    print(diff_obj.C_array)
#    for k in range(0,diff_obj.C_array.shape[0]):
    for k in range(0,3600,50):
        plt.plot(diff_obj.r_array,diff_obj.C_array[0,:],label='t=%.1f'%k)

#    plt.legend()    

    plt.show()
    
