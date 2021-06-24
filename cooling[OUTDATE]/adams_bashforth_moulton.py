import numpy as np
import matplotlib.pyplot as plt
import my_func as mf

class ABM:
    def __init__(self, _h = 0.2, _x = 200, _y0 = 100, env_temp_ = 50):
        self.h = _h
        self.x = np.array([0.0, _x])
        self.y0 = np.array([_y0])
        self.start_temp = _y0
        self.env_temp = env_temp_
        self.f = open('tmp.txt', 'w')

    def __del__(self):
        self.f.close()
        pass

    def RungeKutta4thOrder(self, x):
        y_len = len(self.y0)
        x_len = int((x[-1] - x[0]) / self.h)

        x = x[0]
        y = self.y0

        xsol = np.empty((0))
        xsol = np.append(xsol, x)

        y_res = np.empty((0))
        y_res = np.append(y_res, y)

        for i in range(x_len):
            k1 = mf.myFunc(y, self.env_temp)

            yp2 = y + k1*(self.h/2)

            k2 = mf.myFunc(yp2, self.env_temp)

            yp3 = y + k2*(self.h/2)

            k3 = mf.myFunc(yp3, self.env_temp)

            yp4 = y + k3*self.h

            k4 = mf.myFunc(yp4, self.env_temp)

            for j in range(y_len):
                y[j] = y[j] + (self.h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

            x = x + self.h
            xsol = np.append(xsol, x)

            for r in range(len(y)):
                y_res = np.append(y_res, y[r]) 

        return [xsol, y_res]


    def ABM4thOrder(self):

        y_len = len(self.y0)

        dx = int((self.x[-1] - self.x[0]) / self.h)

        xrk = [self.x[0] + k * self.h for k in range(dx + 1)]

        [xx, yy] = self.RungeKutta4thOrder((xrk[0], xrk[3]))

        self.x = xx
        xsol = np.empty(0)
        xsol = np.append(xsol, self.x)

        y = yy
        yn = np.array([yy[0]])
        y_res = np.empty(0)
        y_res = np.append(y_res, y)

        for i in range(3, dx):
            x00 = self.x[i]; x11 = self.x[i-1]; x22 = self.x[i-2]; x33 = self.x[i-3]; xpp = self.x[i]+self.h

            y00 = np.array([y[i]])
            y11 = np.array([y[i - 1]])
            y22 = np.array([y[i - 2]])
            y33 = np.array([y[i - 3]])

            y0prime = mf.myFunc(y00, self.env_temp)
            y1prime = mf.myFunc(y11, self.env_temp)
            y2prime = mf.myFunc(y22, self.env_temp)
            y3prime = mf.myFunc(y33, self.env_temp)

            ypredictor = y00 + (self.h/24)*(55*y0prime - 59*y1prime + 37*y2prime - 9*y3prime)
            ypp = mf.myFunc(ypredictor, self.env_temp)

            for j in range(y_len):
                yn[j] = y00[j] + (self.h/24)*(9*ypp[j] + 19*y0prime[j] - 5*y1prime[j] + y2prime[j])

            xs = self.x[i] + self.h
            xsol = np.append(xsol, xs)

            self.x = xsol

            for r in range(len(yn)):
                y_res = np.append(y_res, yn)

            y = y_res

        return [xsol, y_res]

    def execute(self):
        [ts, ys] = self.ABM4thOrder()
        for index in ts:
            self.f.write(str(index) + ' ')
        self.f.write('\n')
        for index in ys:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        t = np.arange(0, self.x[-1], 1)
        for index in t:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        for i in t:
            y_math_res = self.env_temp + (self.start_temp - self.env_temp) * np.e**(-(np.log(2) * i/25))
            self.f.write(str(y_math_res) + ' ')
        self.f.write('\n')

        return 0