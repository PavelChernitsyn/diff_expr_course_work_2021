import PySimpleGUI as sg
import BTCS_DirichletBCs
import FTCS_DirichletBCs
import BTCS_NeumannBCs
import CN_NeumannBCs

class MainApplication():
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
                sim.plot_()

            if event == "FTCSDirichlet":
                sim = FTCS_DirichletBCs.FTCSDirichlet(40, 70)
                sim.plot_()

            if event == "BTCSNeumann":
                sim = BTCS_NeumannBCs.BTCSNeumann(50, 50)
                sim.plot_()

            if event == "CrankNicolsonNeumann":
                sim = CN_NeumannBCs.CrankNicolsonNeumann(50, 50)
                sim.plot_()

    def __del__(self):
        self.window.close()