import numpy as np
from scipy import sparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


class BTCSDirichlet:
    
    def __init__(self, M, N):
        self.M = M
        self.N = N
        self.x0 = 0
        self.xL = 1
        self.t0 = 0
        self.tF = 0.2
        self.D = 0.1 # Diffusion coefficient
        self.alpha = -3 # Reaction rate
        self.createGrid()
        self.params()
        self.solve()

        
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
        maindiag = (1 + 2*self.r - self.s)*np.ones((1, self.M-2))
        offdiag = -self.r*np.ones((1, self.M-3))
        a = maindiag.shape[1]
        diagonals = [maindiag, offdiag, offdiag]
        A = sparse.diags(diagonals, [0,-1,1], shape=(a,a)).toarray()
        
        # ----- Initializes matrix U -----
        self.U = np.zeros((self.M, self.N))

        #----- Initial condition -----
        self.U[:,0] = 4*self.xspan - 4*self.xspan**2
        
        #----- Dirichlet boundary conditions -----
        self.U[0,:] = 0.0 
        self.U[-1,:] = 0.0
        
        for k in range(1, self.N):
            c = np.zeros((self.M-4, 1)).ravel()
            b1 = np.asarray([self.r*self.U[0,k], self.r*self.U[-1,k]])
            b1 = np.insert(b1, 1, c)
            b2 = np.array(self.U[1:self.M-1, k-1])
            b = b1 + b2  # Right hand side
            self.U[1:self.M-1, k] = np.linalg.solve(A,b)  # Solve x=A\b

        # ----- Checks if the solution is correct:
        g = np.allclose(np.dot(A,self.U[1:self.M-1,self.N-1]), b)
        print(g)
     
        
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
#     sim = BTCSDirichlet(50, 60)
#     sim.plot_()
    
# if __name__ == "__main__":
#     main()


#M = 50 # GRID POINTS on space interval
#N = 60 # GRID POINTS on time interval