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
        appr = int((self.time - 0)/self.h)

        v = 0
  
        x = 0
        y = self.coef * self.chain_len / (1 + self.coef) + 2 * self.eps/1000

        xsol = np.empty((0))
        xsol = np.append(xsol, x)

        y_res = np.empty((0))
        y_res = np.append(y_res, y)

        for i in range(appr):
            a1, v1 = mf.myFunc(y, v, self.coef, self.chain_len)
            v1 += a1 * self.h
            yp2 = y + v1 * (self.h/2)
            a2, v2 = mf.myFunc(yp2, v1, self.coef, self.chain_len)
            v2 += a2 * self.h
            yp3 = y + v2 * (self.h/2)
            a3, v3 = mf.myFunc(yp3, v2, self.coef, self.chain_len)
            v3 += a3 * self.h
            yp4 = y + v3 * self.h
            a4, v4 = mf.myFunc(yp4, v3, self.coef, self.chain_len)
            v4 += a4 * self.h
            y = y + (self.h/6)*(v1  + 2 * v2  + 2 * v3 + v4 )

            if (y > self.chain_len):
                y = self.chain_len

            x = x + self.h

            v = v1

            xsol = np.append(xsol, x)

            y_res = np.append(y_res, y) 

        return [xsol, y_res]


    def ABM4thOrder(self):
        dx = int((self.time - 0) / self.h)

        xrk = [self.x[0] + k * self.h for k in range(dx + 1)]

        [xx, yy] = self.RungeKutta4thOrder((xrk[0], xrk[3]))

        self.x = xx
        xsol = np.empty(0)

        y = yy
        yn = yy[0]
        y_res = np.empty(0)


        V = 0
        xs = 0

        for i in range(3, dx):
            A1, V1 = mf.myFunc(y[i], V, self.coef, self.chain_len)
            V1 += A1 * self.h
            A2, V2 = mf.myFunc(y[i - 1], V1, self.coef, self.chain_len)
            V2 += A2 * self.h
            A3, V3 = mf.myFunc(y[i - 2], V2, self.coef, self.chain_len)
            V3 += A3 * self.h
            A4, V4 = mf.myFunc(y[i - 3], V3, self.coef, self.chain_len)
            V4 += A4 * self.h

            ypredictor = y[i] + (self.h/24)*(55 * V1 - 59 * V2 + 37 * V3 - 9 * V4 )
            Ap, Vp = mf.myFunc(ypredictor, V, self.coef, self.chain_len)
            Vp += Ap * self.h

            yn = y[i] + (self.h/24) * (9 * Vp + 19 * V1 - 5 * V2 + V3)

            if (yn > self.chain_len):
                yn = self.chain_len
            # print (yn)
            V = V1

            xs = xx[i] + self.h
            xsol = np.append(xsol, xs)

            self.x = xsol

            y_res = np.append(y_res, yn)
        return [xsol, y_res]

    def execute(self):
        [ts, ys] = self.ABM4thOrder()

        for index in ts:
            self.f.write(str(index) + ' ')
        self.f.write('\n')
        
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