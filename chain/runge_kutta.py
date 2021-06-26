
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class Runge_Kutt:
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

    def RKF45(self):
        appr = int((self.time - 0)/self.h)
        
        j = 0

        x = 0
        y = {}
        y[j] = self.coef * 1 / (1 + self.coef) + 2 * self.eps/1000
        
        self.f.write(str(x) + ' ')

        for i in range(appr):
            yp2 = x + mf.myFunc(x, self.coef, self.chain_len)*(self.h/5)
            yp3 = x + mf.myFunc(x, self.coef, self.chain_len)*(3*self.h/40) + mf.myFunc(yp2, self.coef, self.chain_len)*(9*self.h/40)
            yp4 = x + mf.myFunc(x, self.coef, self.chain_len)*(3*self.h/10) - mf.myFunc(yp2, self.coef, self.chain_len)*(9*self.h/10) + mf.myFunc(yp3, self.coef, self.chain_len)*(6*self.h/5)
            yp5 = x - mf.myFunc(x, self.coef, self.chain_len)*(11*self.h/54) + mf.myFunc(yp2, self.coef, self.chain_len)*(5*self.h/2) - mf.myFunc(yp3, self.coef, self.chain_len)*(70*self.h/27) + mf.myFunc(yp4, self.coef, self.chain_len)*(35*self.h/27)
            yp6 = x + mf.myFunc(x, self.coef, self.chain_len)*(1631*self.h/55296) + mf.myFunc(yp2, self.coef, self.chain_len)*(175*self.h/512) + mf.myFunc(yp3, self.coef, self.chain_len)*(575*self.h/13824) + mf.myFunc(yp4, self.coef, self.chain_len)*(44275*self.h/110592) + mf.myFunc(yp5, self.coef, self.chain_len)*(253*self.h/4096)
            
            j += 1
            y[j] = y[j-1] + self.h*(37*mf.myFunc(x, self.coef, self.chain_len)/378 + 22 * self.eps*mf.myFunc(yp3, self.coef, self.chain_len)/621 + 125*mf.myFunc(yp4, self.coef, self.chain_len)/594 + 512*mf.myFunc(yp6, self.coef, self.chain_len)/1771)

            x += self.h
            if (y[j] > 1):
                y[j] = 1
            self.f.write(str(x) + ' ')
        self.f.write('\n')  

        return y

    def execute(self):
        ys = self.RKF45()
        for index in ys:
            self.f.write(str(ys[index]) + ' ')
        self.f.write('\n')

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

        return 0