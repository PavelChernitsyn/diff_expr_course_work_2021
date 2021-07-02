import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class BackwardEuler:
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

    def main_BE(self):
        appr = int((self.time - 0)/self.h)

        v = 0
        y0 = self.coef * self.chain_len / (1 + self.coef) + 2 * self.eps/1000
        x = np.array([y0, v])
        res = []
        xdot = np.array(mf.myFunc(x, self.coef, self.chain_len))
        
        self.f.write(str(0) + ' ')

        for i in range(appr - 1):
            xdot = np.array(mf.myFunc(x + xdot * self.h, self.coef, self.chain_len))

            x = x + xdot * self.h

            if (x[0] > self.chain_len):
                x[0] = self.chain_len   
            res = np.append(res, x[0])
            self.f.write(str(i*self.h) + ' ')
        self.f.write('\n')
        return res

    def execute(self):
        ys = self.main_BE()

        self.f.write(str(self.coef * self.chain_len / (1 + self.coef) + 2 * self.eps/1000) + ' ')
        for index in ys:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        # #Глобальная ошибка в точке t = 0.5
        # y_res = ys[int( / self.h)]
        # y_res_an = 2 * self.eps/2000 * np.exp(np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * 0.5) + 2 * self.eps/2000 * np.exp(-np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * 0.5) + self.coef * self.chain_len / (1 + self.coef)
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