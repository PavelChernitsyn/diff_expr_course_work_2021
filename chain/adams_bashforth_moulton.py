import numpy as np
import matplotlib.pyplot as plt
import my_func as mf

class ABM:
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

    def RungeKutta4thOrder(self, x):
        x_len = int((3 - 0)/self.h)

        x = 0
        y = 0.1 * 1 / (1 + 0.1) + 50/1000

        xsol = np.empty((0))
        xsol = np.append(xsol, x)

        y_res = np.empty((0))
        y_res = np.append(y_res, y)

        for i in range(x_len):
            k1 = mf.myFunc(x)

            yp2 = x + k1*(self.h/2)

            k2 = mf.myFunc(yp2)

            yp3 = x + k2*(self.h/2)

            k3 = mf.myFunc(yp3)

            yp4 = x + k3*self.h

            k4 = mf.myFunc(yp4)

            y = y + (self.h/6)*(k1 + 2*k2 + 2*k3 + k4)

            if (y > 1):
                y = 1

            x = x + self.h
            xsol = np.append(xsol, x)

            y_res = np.append(y_res, y) 

        return [xsol, y_res]


    def ABM4thOrder(self):
        dx = int((3 - 0) / self.h)

        xrk = [self.x[0] + k * self.h for k in range(dx + 1)]

        [xx, yy] = self.RungeKutta4thOrder((xrk[0], xrk[3]))
        print (xx)
        print (yy)

        self.x = xx
        xsol = np.empty(0)
        # xsol = np.append(xsol, self.x)

        y = yy
        yn = yy[0]
        y_res = np.empty(0)
        # y_res = np.append(y_res, y)

        for i in range(3, dx):
            y0prime = mf.myFunc(y[i])
            y1prime = mf.myFunc(y[i - 1])
            y2prime = mf.myFunc(y[i - 2])
            y3prime = mf.myFunc(y[i - 3])

            ypredictor = y[i] + (self.h/24)*(55*y0prime - 59*y1prime + 37*y2prime - 9*y3prime)
            ypp = mf.myFunc(ypredictor)

            yn = y[i] + (self.h/24)*(9*ypp + 19*y0prime - 5*y1prime + y2prime)

            if (yn > 1):
                yn = 1
            # print (yn)

            xs = xx[i] + self.h
            xsol = np.append(xsol, xs)

            self.x = xsol

            y_res = np.append(y_res, yn)
        return [xsol, y_res]

    def execute(self):
        [ts, ys] = self.ABM4thOrder()
        print(len(ts))

        for index in ts:
            self.f.write(str(index) + ' ')
        self.f.write('\n')
        
        for index in ys:
            self.f.write(str(index) + ' ')
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