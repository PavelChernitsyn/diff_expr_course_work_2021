
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf


class Heun:
    def __init__(self, _h = 0.01, _coef = 0.1, _chain_len = 1, _eps = 25, time_ = 3):
        self.h = _h
        self.coef = _coef
        self.chain_len = _chain_len
        self.eps = _eps
        self.time = time_
        self.f = open('tmp.txt', 'w')
    
    def __del__(self):
        self.f.close()
        pass

    def main_H(self):
        appr = int((self.time - 0)/self.h)
        
        j = 0

        x = 0
        y = {}
        y[j] = self.coef * self.chain_len / (1 + self.coef) + 2 * self.eps/1000
        
        self.f.write(str(x) + ' ')

        for i in range(appr):
            j += 1
            y[j] = y[j-1] + (self.h/2)*mf.myFunc(x, self.coef, self.chain_len) + (self.h/2)*mf.myFunc(x + mf.myFunc(x, self.coef, self.chain_len) * self.h, self.coef, self.chain_len)
            print(y[j])
            if (y[j] > 1):
                y[j] = 1

            x += self.h
            self.f.write(str(x) + ' ')
        self.f.write('\n')
        
        return y

    def execute(self):
        ys = self.main_H()
        for index in ys:
            self.f.write(str(ys[index]) + ' ')
        self.f.write('\n')

        #Глобальная ошибка в точке t = 0.5
        y_res = ys[int(0.5 / self.h)]
        y_res_an = T = 2 * self.eps/2000 * np.exp(np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * 0.5) + 2 * self.eps/2000 * np.exp(-np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * 0.5) + self.coef * self.chain_len / (1 + self.coef)
        global_err = y_res - y_res_an

        t = np.arange(0, self.time, self.h)
        for index in t:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        for i in t:
            T = 2 * self.eps/2000 * np.exp(np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * i) + 2 * self.eps/2000 * np.exp(-np.sqrt((1 + self.coef) * 9.8 / self.chain_len) * i) + self.coef * self.chain_len / (1 + self.coef)
            if (T > 1):
                T = 1
            self.f.write(str(T) + ' ')
        self.f.write('\n')

        return global_err