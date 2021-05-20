from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import ttk 
import tkinter as tk
import numpy as np
import backward_euler
import forward_euler
import runge_kutta
import heun

class Plotter(FigureCanvasTkAgg):

    def __init__(self, master):

        self.figure = Figure(dpi=100)
        super().__init__(self.figure, master=master)
        self.axes = self.figure.add_subplot(111)
        self.get_tk_widget().grid(column=0, row=0, sticky='nsew')

    def draw_lists(self, points):

        self.axes.clear()
        # x_list = [x for x in range(0, 100)]
        # x_list = np.arange(1, 2, .01)
        # y_list = [x for x in x_list]
        self.axes.plot(points[0], points[1], color='y')
        self.axes.plot(points[2], points[3], color='b')
        self.draw_idle()

class MainApplication(ttk.Frame):

    def __init__(self, master, *args, **kwargs):

        super().__init__(master)
        self.grid(column=0, row=0, sticky='nsew')

        frame = ttk.Frame(self, borderwidth=8)
        frame.grid(column=0, row=0, sticky='nsew')
        frame.rowconfigure(0, weight=1)

        notes = ttk.Notebook(frame)
        notes.grid(column=0, row=0, sticky='nsew')
        notes.rowconfigure(0, weight=1)

        page = ttk.Frame(notes)
        notes.add(page, text='Picture')


        plot = Plotter(page)

        input_frame = ttk.Frame(self)
        input_frame.grid(column=1, row=0, sticky='nsew')

        # this binding doesn't update the plot
        button_BE = ttk.Button(input_frame, text='Backward Euler', 
                            command=lambda: 
                            plot.draw_lists(backward_euler.BackwardEuler().execute()))
        button_FE = ttk.Button(input_frame, text='Forward Euler', 
                            command=lambda: 
                            plot.draw_lists(forward_euler.ForwardEuler().execute()))
        button_RK = ttk.Button(input_frame, text='Runge Kutt', 
                            command=lambda: 
                            plot.draw_lists(runge_kutta.Runge_Kutt().execute()))
        button_H = ttk.Button(input_frame, text='Heun', 
                            command=lambda: 
                            plot.draw_lists(heun.Heun().execute()))
        button_BE.grid(column=0, row=5, columnspan=2, sticky='ew')
        button_FE.grid(column=0, row=7, columnspan=2, sticky='ew')
        button_RK.grid(column=0, row=9, columnspan=2, sticky='ew')
        button_H.grid(column=0, row=11, columnspan=2, sticky='ew')


# root = tk.Tk() 
# MainApplication(root)
# root.mainloop()