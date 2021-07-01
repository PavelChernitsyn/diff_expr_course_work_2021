
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class Runge_Kutt:
    def __init__(self, _h = 0.001, _coef = 0.1, _chain_len = 1, _eps = 25, time_ = 3):
        self.h = _h
        self.coef = _coef
        self.chain_len = _chain_len
        self.eps = _eps
        self.time = time_
        self.f = open('tmp.txt', 'w')

    def __del__(self):
        self.f.close()
        pass

    def RKF45(self):
        appr = int((self.time - 0)/self.h)
        
        v = 0
        y0 = self.coef * self.chain_len / (1 + self.coef) + 2 * self.eps/1000
        x = np.array([y0, v])
        res = []
        
        self.f.write(str(0) + ' ')

        for i in range(appr):
            xdot1 = np.array(mf.myFunc(x, self.coef, self.chain_len))
            xp2 = x + xdot1 * self.h/5
            xdot2 = np.array(mf.myFunc(xp2, self.coef, self.chain_len))
            xp3 = x + xdot1 * (3*self.h/40) + xdot2 * (9*self.h/40)
            xdot3 = np.array(mf.myFunc(xp3, self.coef, self.chain_len))
            xp4 = x + xdot1 * (3*self.h/10) - xdot2 * (9*self.h/10) + xdot3 * (6*self.h/5)
            xdot4 = np.array(mf.myFunc(xp4, self.coef, self.chain_len))
            xp5 = x - xdot1 *(11*self.h/54) + xdot2 * (5*self.h/2) - xdot3 * (70*self.h/27) + xdot4 * (35*self.h/27)
            xdot5 = np.array(mf.myFunc(xp5, self.coef, self.chain_len))
            xp6 = x + xdot1 * (1631*self.h/55296) + xdot2 * (175*self.h/512) + xdot3 * (575*self.h/13824) + xdot4 * (44275*self.h/110592) + xdot5 * (253*self.h/4096)
            xdot6 = np.array(mf.myFunc(xp6, self.coef, self.chain_len))

            x = x + self.h * (37* xdot1 / 378 + 250 * xdot3 / 621 + 125 * xdot4 / 594 + 512 * xdot6 / 1771)

            if (x[0] > self.chain_len):
                x[0] = self.chain_len   
            res = np.append(res, x[0]) 
            self.f.write(str(i*self.h) + ' ')
        self.f.write('\n')  

        return res

    def execute(self):
        ys = self.RKF45()

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