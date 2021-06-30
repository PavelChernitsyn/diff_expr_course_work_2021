from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import ttk 
import tkinter as tk
import numpy as np
import backward_euler
import forward_euler
import runge_kutta
import heun
import adams_bashforth_moulton
import os

class Plotter(FigureCanvasTkAgg):

    def __init__(self, master):

        self.figure = Figure(dpi=100)
        super().__init__(self.figure, master=master)
        self.axes = self.figure.add_subplot(111)
        self.get_tk_widget().grid(column=0, row=0, sticky='nsew')

    def draw_lists(self):

        self.axes.clear()
    
        x = [[],[],[],[]]
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

        self.axes.plot(x[0], x[1], color='y')
        self.axes.plot(x[2], x[3], color='b')
        self.axes.set_xlabel('T, sec')
        self.axes.set_ylabel('L, mm')
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


        self.plot = Plotter(page)

        input_frame = ttk.Frame(self)
        input_frame.grid(column=1, row=0, sticky='nsew')

        label_chain_lenght = ttk.Label(input_frame)
        self.slider_chain_lenght = ttk.Scale(input_frame, from_ = 0.1, to_ = 5,
                            command=lambda x:
                            label_chain_lenght.config(text = "Chain lenght = " + 
                            str("{0:.1f}".format(self.slider_chain_lenght.get())) + "m"))
        self.slider_chain_lenght.set(1)
        label_chain_lenght.config(text = "Chain lenght = " + str("{0:.1f}".format(self.slider_chain_lenght.get())) + "m")

        label_epsilon = ttk.Label(input_frame)
        self.slider_epsilon = ttk.Scale(input_frame, from_ = 1, to_ = 100,
                            command=lambda x:
                            label_epsilon.config(text = "Epsilon = " + str(int(self.slider_epsilon.get())) + "mm"))
        self.slider_epsilon.set(25)
        label_epsilon.config(text = "Epsilon = " + str(int(self.slider_epsilon.get())) + "mm")

        label_coef = ttk.Label(input_frame)
        self.slider_coef = ttk.Scale(input_frame, from_ = 0, to_ = 0.3,
                            command=lambda x:
                            label_coef.config(text = "Coefficient of friction = " + 
                            str("{0:.2f}".format(self.slider_coef.get()))))
        self.slider_coef.set(0.1)
        label_coef.config(text = "Coefficient of friction = " + str("{0:.1f}".format(self.slider_coef.get())))

        label_time = ttk.Label(input_frame)
        self.slider_time = ttk.Scale(input_frame, from_ = 0.1, to_ = 10,
                            command=lambda x:
                            label_time.config(text = "Time = " + 
                            str("{0:.1f}".format(self.slider_time.get())) + "sec"))
        self.slider_time.set(3)
        label_time.config(text = "Time = " + str(int(self.slider_time.get())) + "sec")

        self.label_err_msg = ttk.Label(input_frame)
        self.label_err_msg.config(text = "")

        button_BE = ttk.Button(input_frame, text='Backward Euler', command = self.button_BE_clicked)
                            
        button_FE = ttk.Button(input_frame, text='Forward Euler', command = self.button_FE_clicked)

        button_RK = ttk.Button(input_frame, text='Runge Kutt', command = self.button_RK_clicked)
        
        button_H = ttk.Button(input_frame, text='Heun', command = self.button_H_clicked)
        
        button_ADM = ttk.Button(input_frame, text='Adams-Bashforth-Moulton', command = self.button_ADM_clicked)

        button_discard_param = ttk.Button(input_frame, text='Discard params', command = self.button_discard_param_clicked)
                            
        button_BE.grid(column=0, row=0, columnspan=2, sticky='ew')
        button_FE.grid(column=0, row=1, columnspan=2, sticky='ew')
        button_RK.grid(column=0, row=2, columnspan=2, sticky='ew')
        button_H.grid(column=0, row=3, columnspan=2, sticky='ew')
        button_ADM.grid(column=0, row=4, columnspan=2, sticky='ew')
        label_chain_lenght.grid(column=0, row=5, columnspan=2, sticky='ew')
        self.slider_chain_lenght.grid(column=0, row=6, columnspan=2, sticky='ew')
        label_epsilon.grid(column=0, row=7, columnspan=2, sticky='ew')
        self.slider_epsilon.grid(column=0, row=8, columnspan=2, sticky='ew')
        label_coef.grid(column=0, row=9, columnspan=2, sticky='ew')
        self.slider_coef.grid(column=0, row=10, columnspan=2, sticky='ew')
        label_time.grid(column=0, row=11, columnspan=2, sticky='ew')
        self.slider_time.grid(column=0, row=12, columnspan=2, sticky='ew')
        button_discard_param.grid(column=0, row=13, columnspan=2, sticky='ew')
        self.label_err_msg.grid(column=0, row=14, columnspan=2, sticky='ew')

    def button_BE_clicked(self):
        err = backward_euler.BackwardEuler(
            0.001, self.slider_coef.get(), self.slider_chain_lenght.get(), 

            self.slider_epsilon.get(), self.slider_time.get()
            ).execute() 
        self.plot.draw_lists()
        self.label_err_msg.config(text = "Err. = " + str("{0:.8f}".format(err)))

    def button_FE_clicked(self):
        err = forward_euler.ForwardEuler(
            0.001, self.slider_coef.get(), self.slider_chain_lenght.get(), 

            self.slider_epsilon.get(), self.slider_time.get()
            ).execute()
        self.plot.draw_lists()
        self.label_err_msg.config(text = "Err. = " + str("{0:.8f}".format(err)))

    def button_RK_clicked(self):
        err = runge_kutta.Runge_Kutt(
            0.001, int(self.slider_coef.get()), self.slider_chain_lenght.get(), 
            self.slider_epsilon.get(), self.slider_time.get()
            ).execute()
        self.plot.draw_lists()
        self.label_err_msg.config(text = "Err. = " + str("{0:.8f}".format(err)))

    def button_H_clicked(self):
        err = heun.Heun(
            0.001, int(self.slider_coef.get()), self.slider_chain_lenght.get(), 
            self.slider_epsilon.get(), self.slider_time.get()
            ).execute()
        self.plot.draw_lists()
        self.label_err_msg.config(text = "Err. = " + str("{0:.8f}".format(err)))

    def button_ADM_clicked(self):
        err = adams_bashforth_moulton.ABM(
            0.001, int(self.slider_coef.get()), self.slider_chain_lenght.get(), 
            self.slider_epsilon.get(), self.slider_time.get()
            ).execute()
        self.plot.draw_lists()
        self.label_err_msg.config(text = "Err. = " + str("{0:.8f}".format(err)))

    def button_discard_param_clicked(self):
        self.slider_epsilon.set(2 * self.eps)
        self.slider_chain_lenght.set(100)
        self.slider_coef.set(200)
        self.check_sliders()

    def check_sliders(self):
        # print("hi")
        if (self.slider_chain_lenght.get() <= self.slider_epsilon.get()):
            self.label_err_msg.config(text = "\n\nStatus: ERROR!\nTop Temp <= Env Temp")
            return False
        else:
            self.label_err_msg.config(text = "\n\nStatus: DONE!")
            return True
    
    def __del__(self):
        # os.remove('tmp.txt')
        pass