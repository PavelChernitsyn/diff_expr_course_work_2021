import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class BackwardEuler:
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

    def main_BE(self):
        appr = int((3 - 0)/self.h)

        j = 0

        x = 0
        y = {}
        y[j] = 0.1 * 1 / (1 + 0.1) + 50/1000

        self.f.write(str(x) + ' ')

        for i in range(appr):
            F_x_t = mf.myFunc(x)/(1+self.h)

            j += 1
            y[j] = y[j-1] + self.h*F_x_t
            print(y[j])

            x += self.h
            if (y[j] > 1):
                y[j] = 1
            self.f.write(str(x) + ' ')
        self.f.write('\n')

        return y

    def execute(self):
        ys = self.main_BE()
        
        for index in ys:
            self.f.write(str(ys[index]) + ' ')
        self.f.write('\n')

        t = np.arange(0, 3, 0.1)

        for index in t:
            self.f.write(str(index) + ' ')
        self.f.write('\n')
        
        for i in t:
            # T = np.sqrt(self.chain_len/((1+self.coef)*9.8)) * np.log(self.chain_len/((1+self.coef)*i/1000) 
            #     + np.sqrt(self.chain_len**2/(((1+self.coef)**2)*((i/1000)**2)) - 1))
            T = 50/2000 * np.exp(np.sqrt((1 + 0.1) * 9.8 / 1) * i) + 50/2000 * np.exp(-np.sqrt((1 + 0.1) * 9.8 / 1) * i) + 0.1 * 1 / (1 + 0.1)
            if (T > 1):
                T = 1
            # print("t = ", i , " , x(t) = ", T)
            self.f.write(str(T) + ' ')
        self.f.write('\n')

        return 0