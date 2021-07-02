from os import TMP_MAX
from tkinter.constants import X
import numpy as np
import matplotlib.pyplot as plt
import my_func as mf

class ABM:

    def __init__(self, _h = 0.001, _coef = 0.1, _chain_len = 1, _eps = 25, time_ = 3):
        self.h = _h
        self.coef = _coef
        self.chain_len = _chain_len
        self.eps = _eps
        self.x = np.array([0.0])
        self.time = time_
        self.f = open('tmp.txt', 'w')

    def __del__(self):
        self.f.close()
        pass

    def RungeKutta4thOrder(self, x):
        appr = int((x[-1] - x[0])/self.h)

        v = 0
        y0 = self.coef * self.chain_len / (1 + self.coef) + 2 * self.eps/1000
        X = np.array([y0, v])
        res_v = []
        res_x = []
        res_x = np.append(res_x, X[0])
        res_v = np.append(res_v, X[1])

        for i in range(appr):
            xdot1 = np.array(mf.myFunc(X, self.coef, self.chain_len))
            xp2 = X + xdot1 * (self.h/2)
            xdot2 = np.array(mf.myFunc(xp2, self.coef, self.chain_len))
            xp3 = X + xdot2 * (self.h/2)
            xdot3 = np.array(mf.myFunc(xp3, self.coef, self.chain_len))
            xp4 = X + xdot3 * self.h
            xdot4 = np.array(mf.myFunc(xp4, self.coef, self.chain_len))
            X = X + (self.h/6)*(xdot1  + 2 * xdot2  + 2 * xdot3 + xdot4)

            res_x = np.append(res_x, X[0])
            res_v = np.append(res_v, X[1])
        return [res_x, res_v]


    def ABM4thOrder(self):
        dx = int((self.time - 0) / self.h)

        xrk = [self.x[0] + k * self.h for k in range(dx + 1)]

        [res_x, res_v] = self.RungeKutta4thOrder((xrk[0], xrk[3]))

        t_arr = np.empty(0)

        x = np.array([[res_x[0], res_v[0]]])

        for i in range(1, len(res_x)):
            x = np.append(x, [[res_x[i], res_v[i]]], axis=0)
        xn = np.array([res_x[0], res_v[0]])
        print(x)

        res = []
        res = np.append(res, res_x)

        t = self.h * 4
        t_arr = np.append(t_arr, 0)
        t_arr = np.append(t_arr, self.h * 1)
        t_arr = np.append(t_arr, self.h * 2)
        t_arr = np.append(t_arr, self.h * 3)

        for i in range(3, dx):
            
            xdot1 = np.array(mf.myFunc(x[i], self.coef, self.chain_len))
            xdot2 = np.array(mf.myFunc(x[i-1], self.coef, self.chain_len))
            xdot3 = np.array(mf.myFunc(x[i-2], self.coef, self.chain_len))
            xdot4 = np.array(mf.myFunc(x[i-3], self.coef, self.chain_len))

            xp = x[i] + (self.h/24)*(55 * xdot1 - 59 * xdot2 + 37 * xdot3 - 9 * xdot4 )
            xdotp = np.array(mf.myFunc(xp, self.coef, self.chain_len))

            xn = x[i] + (self.h/24) * (9 * xdotp + 19 * xdot1 - 5 * xdot2 + xdot3)

            if (xn[0] > self.chain_len):
                xn[0] = self.chain_len

            t += self.h
            t_arr = np.append(t_arr, t)

            res = np.append(res, xn[0])
            x = np.append(x, [xn], axis = 0)

        return [t_arr, res]

    def execute(self):
        [ts, ys] = self.ABM4thOrder()

        self.f.write(str(0) + ' ')
        for index in ts:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        self.f.write(str(self.coef * self.chain_len / (1 + self.coef) + 2 * self.eps/1000) + ' ')        
        for index in ys:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        # #Глобальная ошибка в точке t = 0.5
        # y_res = ys[int(0.5 / self.h)]
        # y_res_an = T = 2 * self.eps/2000 * np.exp(np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * 0.5) + 2 * self.eps/2000 * np.exp(-np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * 0.5) + self.coef * self.chain_len / (1 + self.coef)
        # global_err = y_res - y_res_an

        t = np.arange(0, self.time + self.h, self.h)
        for index in t:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        for i in t:
            T = 2 * self.eps/2000 * np.exp(np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * i) + 2 * self.eps/2000 * np.exp(-np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * i) + self.coef * self.chain_len / (1 + self.coef)
            if (T > self.chain_len):
                T = self.chain_len
            self.f.write(str(T) + ' ')
        self.f.write('\n')

        # return global_err