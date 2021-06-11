import PySimpleGUI as sg
import BTCS_DirichletBCs
import FTCS_DirichletBCs
import BTCS_NeumannBCs
import CN_NeumannBCs
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

class MainApplication():
    def plotter(self):
        x = [[],[],[],[]] #tspan, xspan, tS, yexact
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
        # os.remove("tmp.txt")

        U = np.loadtxt("tmp_U.txt")
        
        X, T = np.meshgrid(x[0], x[1])

        if (len(x[0]) < len(x[1])):
            cut = len(x[0])
        else:
            cut = len(x[1])

        plt.plot(x[0][:cut], U[:, 1], 'r')
        plt.plot(x[2], x[3], 'b')
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

        os.remove("tmp.txt")
        os.remove("tmp_U.txt")

    def __init__(self):
        layout = [[sg.Button("BTCSDirichlet")], [sg.Button("FTCSDirichlet")], 
                  [sg.Button("BTCSNeumann")], [sg.Button("CrankNicolsonNeumann")]]

        # Create the window
        self.window = sg.Window("Diffusion", layout)

        # Create an event loop
        while True:
            event, values = self.window.read()
            # End program if user closes window
            if event == sg.WIN_CLOSED:
                break

            if event == "BTCSDirichlet":
                sim = BTCS_DirichletBCs.BTCSDirichlet(50, 50)
                self.plotter()

            if event == "FTCSDirichlet":
                sim = FTCS_DirichletBCs.FTCSDirichlet(40, 70)
                self.plotter()

            if event == "BTCSNeumann":
                sim = BTCS_NeumannBCs.BTCSNeumann(50, 50)
                self.plotter()

            if event == "CrankNicolsonNeumann":
                sim = CN_NeumannBCs.CrankNicolsonNeumann(50, 50)
                self.plotter()

    def __del__(self):
        self.window.close()