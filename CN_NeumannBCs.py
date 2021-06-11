import numpy as np
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


class CrankNicolsonNeumann:
    
    def __init__(self, M_, N_, x0_ = 0, xL_ = 1, t0_ = 0, tF_ = 1, D_ = 0.1, alpha_ = -3):
        self.M = M_
        self.N = N_
        self.x0 = x0_
        self.xL = xL_
        self.t0 = t0_
        self.tF = tF_
        self.D = D_ # Diffusion coefficient
        self.alpha = alpha_ # Reaction rate
        self.f = open('tmp.txt', 'w')
        self.params()
        self.createGrid()
        self.matrices()
        self.solve()
        
    def __del__(self):
        # self.f.close()
        pass
        
    def params(self):
        # ----- Spatial discretization step -----
        self.dx = (self.xL - self.x0)/(self.M - 1)
        # ----- Time step -----
        dt = (self.tF - self.t0)/(self.N - 1)
        self.r = dt*self.D/(2*self.dx**2)
        self.s = dt*self.alpha/2
        self.a0 = 1 + 2*self.r - self.s
        self.c0 = 1 - 2*self.r + self.s
        
        
    def createGrid(self):
        self.xspan = np.linspace(self.x0, self.xL, self.M)
        self.tspan = np.linspace(self.t0, self.tF, self.N)

        for elem in self.tspan:
            self.f.write(str(elem) + ' ')
        self.f.write('\n')

        for elem in self.xspan:
            self.f.write(str(elem) + ' ')
        self.f.write('\n')
        
    def matrices(self):
        main_diag_a0 = self.a0*np.ones((1, self.M))
        off_diag_a0 = -self.r*np.ones((1, self.M-1))
        
        main_diag_c0 = self.c0*np.ones((1, self.M))
        off_diag_c0 = self.r*np.ones((1, self.M-1))
        
        # Left-hand side tri-diagonal matrix
        a = main_diag_a0.shape[1]
        diagonalsA = [main_diag_a0, off_diag_a0, off_diag_a0]
        self.A = sparse.diags(diagonalsA, [0,-1,1], shape=(a,a)).toarray()
        self.A[0, 1] = -2*self.r
        self.A[self.M-1, self.M-2] = -2*self.r
        
        # Right-hand side tri-diagonal matrix
        c = main_diag_c0.shape[1]
        diagonalsC = [main_diag_c0, off_diag_c0, off_diag_c0]
        self.A_rhs = sparse.diags(diagonalsC, [0,-1,1], shape=(c,c)).toarray()
        self.A_rhs[0, 1] = 2*self.r
        self.A_rhs[self.M-1, self.M-2] = 2*self.r
        
        
    def solve(self):
        
        # ----- Initializes matrix U -----
        self.U = np.zeros((self.M, self.N))
        
        #----- Initial condition -----
        self.U[:,0] = 4*self.xspan - 4*self.xspan**2
        
        #----- Neumann boundary conditions -----
        leftBC = np.arange(1, self.N+1)
        f = np.sin(leftBC*np.pi/2)

        rightBC = np.arange(1, self.N+1)
        g = np.sin(3*rightBC*np.pi/4)
        
        for k in range(1, self.N):
            ins = np.zeros((self.M-2, 1)).ravel()
            b1 = np.asarray([4*self.r*self.dx*f[k], 4*self.r*self.dx*g[k]])
            b1 = np.insert(b1, 1, ins)
            b2 = np.matmul(self.A_rhs, np.array(self.U[0:self.M, k-1]))
            b = b1 + b2  # Right hand side
            self.U[0:self.M, k] = np.linalg.solve(self.A,b)  # Solve x=A\b
        np.savetxt("tmp_U.txt", self.U)

        dtS = int((self.xspan[-1] - self.xspan[0])/(self.D))
        tS = [self.xspan[0]+i*(self.D) for i in range(int(dtS)+1)]
        
        for elem in tS:
            self.f.write(str(elem) + ' ')
        self.f.write('\n')

        for i in range(dtS+1):
            ye = 4*tS[i] - 4*tS[i]**2
            self.f.write(str(ye) + ' ')
        self.f.write('\n')
        # ----- Checks if the solution is correct:
        gc = np.allclose(np.dot(self.A,self.U[0:self.M,self.N-1]), b)

        self.f.close()