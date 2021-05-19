import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg
import matplotlib
import backward_euler

fig = matplotlib.figure.Figure(figsize=(5, 4), dpi=100)
t = np.arange(1, 2, .01)
# fig.add_subplot(111).plot(ts, ys)
# fig.add_subplot(111).plot(t, yexact)
fig.clear()
matplotlib.use("TkAgg")

def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw_idle()
    #figure_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=1)
    return figure_canvas_agg

layout = [
    [sg.Text("Diffusion")],
    [sg.Canvas(key="-CANVAS-")],
    [sg.Button("Backward Euler"), sg.Button("Forward Euler"), 
     sg.Button("Runge-Kutta"), sg.Button("Heun")],
]

# Create the form and show it without the plot
window = sg.Window(
    "Course Work",
    layout,
    location=(0, 0),
    finalize=True,
    element_justification="center",
    font="Helvetica 18",
)

# Add the plot to the window
#can = draw_figure(window["-CANVAS-"].TKCanvas, fig)

while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED:
        break
    elif event == "Backward Euler":
        points = backward_euler.BackwardEuler().execute()
        fig.add_subplot(111).plot(points[0], points[1])
        fig.add_subplot(111).plot(points[2], points[3])
        draw_figure(window["-CANVAS-"].TKCanvas, fig)
        print("xuy1")
    elif event == "Forward Euler":
        print("xuy2")
    elif event == "Runge-Kutta":
        print("xuy3")
    elif event == "Heun":
        print("xuy4")

window.close()