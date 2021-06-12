import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


class FTCSDirichlet:
    
    def __init__(self, M_, N_, x0_ = 0, xL_ = 1, t0_ = 0, tF_ = 0.2, D_ = 0.1, alpha_ = -3):
        self.M = M_
        self.N = N_
        self.x0 = x0_
        self.xL = xL_
        self.t0 = t0_
        self.tF = tF_
        self.D = D_ # Diffusion coefficient
        self.alpha = alpha_ # Reaction rate
        self.f = open('tmp.txt', 'w')
        self.createGrid()
        self.params()
        self.solve()
        
    def __del__(self):
        # self.f.close()
        pass
        
    def createGrid(self):
        self.xspan = np.linspace(self.x0, self.xL, self.M)
        self.tspan = np.linspace(self.t0, self.tF, self.N)

    def params(self):
        # ----- Spatial discretization step -----
        dx = (self.xL - self.x0)/(self.M - 1)
        # ----- Time step -----
        dt = (self.tF - self.t0)/(self.N - 1)
        self.r = dt*self.D/dx**2
        self.s = dt*self.alpha
        
    def solve(self):
        # ----- Initializes matrix solution U -----
        self.U = np.zeros((self.M, self.N))
        # ----- Initial condition -----
        self.U[:,0] = 4*self.xspan - 4*self.xspan**2
        # ----- Dirichlet Boundary Conditions -----
        self.U[0,:] = 0.0
        self.U[-1,:] = 0.0
        
        for k in range(0, self.N-1):
            for i in range(1, self.M-1):
                self.U[i, k+1] = self.r*self.U[i-1, k] + (1-2*self.r+self.s)*self.U[i,k] + self.r*self.U[i+1,k] 
                
        np.savetxt("tmp_U.txt", self.U)

        self.T, self.X = np.meshgrid(self.tspan, self.xspan)
        
        dtS = int((self.xspan[-1] - self.xspan[0])/(self.D))
        tS = [self.xspan[0]+i*(self.D) for i in range(int(dtS)+1)]
        
        for elem in tS:
            self.f.write(str(elem) + ' ')
        self.f.write('\n')

        for i in range(dtS+1):
            ye = 4*tS[i] - 4*tS[i]**2
            self.f.write(str(ye) + ' ')
        self.f.write('\n')

        self.f.close()