import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import os

class Plotter:
    def __init__(self, M_ = 50, N_ = 50, x0_ = 0, xL_ = 1, t0_ = 0, tF_ = 1):
        self.M = M_
        self.N = N_
        self.x0 = x0_
        self.xL = xL_
        self.t0 = t0_
        self.tF = tF_
        self.xspan = np.linspace(self.x0, self.xL, self.M)
        self.tspan = np.linspace(self.t0, self.tF, self.N)

        X, T = np.meshgrid(self.tspan, self.xspan)

        np.savetxt("tmp_meshgrid_X.txt", X)
        np.savetxt("tmp_meshgrid_T.txt", T)
        # self.plot_()

    # def __del__(self):
    #     os.remove("tmp_grid.txt")
    #     os.remove("tmp_meshgrid_X.txt")
    #     os.remove("tmp_meshgrid_T.txt")
    #     os.remove("tmp_U.txt")
    #     os.remove("tmp.txt")

    def plot_(self):
        x = [[],[]] #tS, yexact
        i = 0

        f = open('tmp.txt', 'r')
        for line in f:
            if (line == '\n'):
                    continue
            for elem in line.split(' '):
                if (elem == '\n'):
                    continue
                x[i].append(float(elem))
            i+=1
        f.close()

        U = np.loadtxt("tmp_U.txt")
        X = np.loadtxt("tmp_meshgrid_X.txt")
        T = np.loadtxt("tmp_meshgrid_T.txt")


        if (len(self.xspan) < len(self.tspan)):
            cut = len(self.xspan)
        else:
            cut = len(self.tspan)

        plt.plot(self.xspan[:cut], U[:, 1], 'r')
        plt.plot(x[0], x[1], 'b')
        plt.show()

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_surface(X, T, U, linewidth=0,
                               cmap=cm.coolwarm, antialiased=False)
        
        # ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        ax.set_xlabel('Time')
        ax.set_ylabel('Space')
        ax.set_zlabel('U')
        plt.tight_layout()
        plt.show()

        # os.remove("tmp.txt")
        # os.remove("tmp_U.txt")