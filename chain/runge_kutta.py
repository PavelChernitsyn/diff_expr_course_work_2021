
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class Runge_Kutt:
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

    def RKF45(self):
        appr = int((3 - 0)/self.h)
        
        j = 0

        x = 0
        y = {}
        y[j] = 0.1 * 1 / (1 + 0.1) + 50/1000
        
        self.f.write(str(x) + ' ')

        for i in range(appr):
            yp2 = x + mf.myFunc(x)*(self.h/5)
            yp3 = x + mf.myFunc(x)*(3*self.h/40) + mf.myFunc(yp2)*(9*self.h/40)
            yp4 = x + mf.myFunc(x)*(3*self.h/10) - mf.myFunc(yp2)*(9*self.h/10) + mf.myFunc(yp3)*(6*self.h/5)
            yp5 = x - mf.myFunc(x)*(11*self.h/54) + mf.myFunc(yp2)*(5*self.h/2) - mf.myFunc(yp3)*(70*self.h/27) + mf.myFunc(yp4)*(35*self.h/27)
            yp6 = x + mf.myFunc(x)*(1631*self.h/55296) + mf.myFunc(yp2)*(175*self.h/512) + mf.myFunc(yp3)*(575*self.h/13824) + mf.myFunc(yp4)*(44275*self.h/110592) + mf.myFunc(yp5)*(253*self.h/4096)
            
            j += 1
            y[j] = y[j-1] + self.h*(37*mf.myFunc(x)/378 + 250*mf.myFunc(yp3)/621 + 125*mf.myFunc(yp4)/594 + 512*mf.myFunc(yp6)/1771)

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