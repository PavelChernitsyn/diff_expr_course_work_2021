
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
        
        j = 0

        v = 0
        x = 0
        y = {}
        y[j] = self.coef * 1 / (1 + self.coef) + 2 * self.eps/1000
        
        self.f.write(str(x) + ' ')

        for i in range(appr):
            a1, v1 = mf.myFunc(y[j], v, self.coef, self.chain_len)
            v1 += a1 * self.h
            yp2 = y[j] + v1 * self.h/5
            a2, v2 = mf.myFunc(yp2, v1, self.coef, self.chain_len)
            v2 += a2 * self.h
            yp3 = y[j] + v1 * (3*self.h/40) + v2 * (9*self.h/40)
            a3, v3 = mf.myFunc(yp3, v2, self.coef, self.chain_len)
            v3 += a3 * self.h
            yp4 = y[j] + v1 * (3*self.h/10) - v2 * (9*self.h/10) + v3 * (6*self.h/5)
            a4, v4 = mf.myFunc(yp4, v3, self.coef, self.chain_len)
            v4 += a4 * self.h
            yp5 = y[j] - v1 *(11*self.h/54) + v2 * (5*self.h/2) - v3 * (70*self.h/27) + v4 * (35*self.h/27)
            a5, v5 = mf.myFunc(yp5, v4, self.coef, self.chain_len)
            v5 += a5 * self.h
            yp6 = y[j] + v1 * (1631*self.h/55296) + v2 * (175*self.h/512) + v3 * (575*self.h/13824) + v4 * (44275*self.h/110592) + v5 * (253*self.h/4096)
            a6, v6 = mf.myFunc(yp6, v5, self.coef, self.chain_len)
            v6 += a6 * self.h

            j += 1
            y[j] = y[j-1] + self.h * (37* v1 / 378 + 250 * v3 / 621 + 125 * v4 / 594 + 512 * v6 / 1771)

            v = v1
            x += self.h
            if (y[j] > self.chain_len):
                y[j] = self.chain_len
            self.f.write(str(x) + ' ')
        self.f.write('\n')  

        return y

    def execute(self):
        ys = self.RKF45()
        for index in ys:
            self.f.write(str(ys[index]) + ' ')
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