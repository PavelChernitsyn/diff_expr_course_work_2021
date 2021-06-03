import numpy as np
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


class BTCSNeumann:
    
    def __init__(self, M, N):
        self.M = M
        self.N = N
        self.x0 = 0
        self.xL = 1
        self.t0 = 0
        self.tF = 0.2
        self.D = 0.5  # Diffusion coefficient
        self.alpha = -5 # Reaction rate
        self.params()
        self.createGrid()
        self.solve()
        
        
    def params(self):
        # ----- Spatial discretization step -----
        self.dx = (self.xL - self.x0)/(self.M - 1)
        # ----- Time step -----
        dt = (self.tF - self.t0)/(self.N - 1)
        self.r = dt*self.D/self.dx**2
        self.s = dt*self.alpha
        
        
    def createGrid(self):
        self.xspan = np.linspace(self.x0, self.xL, self.M)
        self.tspan = np.linspace(self.t0, self.tF, self.N)
        
        
    def solve(self):
        maindiag = (1 + 2*self.r - self.s)*np.ones((1,self.M))
        offdiag = -self.r*np.ones((1, self.M-1))
        a = maindiag.shape[1]
        diagonals = [maindiag, offdiag, offdiag]
        A = sparse.diags(diagonals, [0,-1,1], shape=(a,a)).toarray()
        A[0, 1] = -2*self.r
        A[self.M-1, self.M-2] = -2*self.r
        
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
            c = np.zeros((self.M-2, 1)).ravel()
            b1 = np.asarray([2*self.r*self.dx*f[k], 2*self.r*self.dx*g[k]])
            b1 = np.insert(b1, 1, c)
            b2 = np.array(self.U[0:self.M, k-1])
            b = b1 + b2  # Right hand side
            self.U[0:self.M, k] = np.linalg.solve(A,b)  # Solve x=A\b
            
        # ----- Checks if the solution is correct:
        gc = np.allclose(np.dot(A,self.U[0:self.M,self.N-1]), b)
        print(gc)
        
    
    def plot_(self):
        # ----- Surface plot -----
        X, T = np.meshgrid(self.tspan, self.xspan)

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_surface(X, T, self.U, linewidth=0,
                       cmap=cm.coolwarm, antialiased=False)

        ax.set_xticks([0, 0.05, 0.1, 0.15, 0.2])

        ax.set_xlabel('Time')
        ax.set_ylabel('Space')
        ax.set_zlabel('U')
        plt.tight_layout()
        plt.show()


# def main():
#     sim = BTCSNeumann(50, 60)
#     sim.plot_()
    
# if __name__ == "__main__":
#     main()