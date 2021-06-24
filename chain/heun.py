
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf


class Heun:
    def __init__(self, _h = 0.01, _coef = 0.1, _chain_len = 1, _eps = 25, _x = 200):
        self.h = _h
        self.coef = _coef
        self.chain_len = _chain_len
        self.eps = _eps
        self.x = np.array([0.0, _x])
        self.f = open('tmp.txt', 'w')
    
    def __del__(self):
        self.f.close()
        pass

    def main_H(self):
        appr = int((3 - 0)/self.h)
        
        j = 0

        x = 0
        y = {}
        y[j] = 0.1 * 1 / (1 + 0.1) + 50/1000
        
        self.f.write(str(x) + ' ')

        for i in range(appr):
            j += 1
            y[j] = y[j-1] + (self.h/2)*mf.myFunc(x) + (self.h/2)*mf.myFunc(x + mf.myFunc(x) * self.h)
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

        t = np.arange(0, 3, 0.1)
        for index in t:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        for i in t:
            T = 50/2000 * np.exp(np.sqrt((1 + 0.1) * 9.8 / 1) * i) + 50/2000 * np.exp(-np.sqrt((1 + 0.1) * 9.8 / 1) * i) + 0.1 * 1 / (1 + 0.1)
            if (T > 1):
                T = 1
            self.f.write(str(T) + ' ')
        self.f.write('\n')

        return 0