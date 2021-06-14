import numpy as np
import matplotlib.pyplot as plt
import my_func as mf

class ABM:
    def __init__(self, _h = 0.2, _x = np.array([1.0, 200.0]), _yinit = np.array([80.0])):
        self.h = _h
        self.x = _x
        self.y0 = np.array([80.0])
        self.f = open('tmp.txt', 'w')

    def __del__(self):
        self.f.close()
        pass

    def RungeKutta4thOrder(self, y0, x, h):
        y_len = len(y0)
        x_len = int((x[-1] - x[0]) / h)

        x = x[0]
        y = y0

        xsol = np.empty((0))
        xsol = np.append(xsol, x)

        y_res = np.empty((0))
        y_res = np.append(y_res, y)

        for i in range(x_len):
            k1 = mf.myFunc(y)

            yp2 = y + k1*(h/2)

            k2 = mf.myFunc(yp2)

            yp3 = y + k2*(h/2)

            k3 = mf.myFunc(yp3)

            yp4 = y + k3*h

            k4 = mf.myFunc(yp4)

            for j in range(y_len):
                y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

            x = x + h
            xsol = np.append(xsol, x)

            for r in range(len(y)):
                y_res = np.append(y_res, y[r]) 

        return [xsol, y_res]


    def ABM4thOrder(self, y0, x, h):

        y_len = len(y0)

        dx = int((x[-1] - x[0]) / h)

        xrk = [x[0] + k * h for k in range(dx + 1)]

        [xx, yy] = self.RungeKutta4thOrder(y0, (xrk[0], xrk[3]), h)

        x = xx
        xsol = np.empty(0)
        xsol = np.append(xsol, x)

        y = yy
        yn = np.array([yy[0]])
        y_res = np.empty(0)
        y_res = np.append(y_res, y)

        for i in range(3, dx):
            x00 = x[i]; x11 = x[i-1]; x22 = x[i-2]; x33 = x[i-3]; xpp = x[i]+h

            y00 = np.array([y[i]])
            y11 = np.array([y[i - 1]])
            y22 = np.array([y[i - 2]])
            y33 = np.array([y[i - 3]])

            y0prime = mf.myFunc(y00)
            y1prime = mf.myFunc(y11)
            y2prime = mf.myFunc(y22)
            y3prime = mf.myFunc(y33)

            ypredictor = y00 + (h/24)*(55*y0prime - 59*y1prime + 37*y2prime - 9*y3prime)
            ypp = mf.myFunc(ypredictor)

            for j in range(y_len):
                yn[j] = y00[j] + (h/24)*(9*ypp[j] + 19*y0prime[j] - 5*y1prime[j] + y2prime[j])

            xs = x[i] + h
            xsol = np.append(xsol, xs)

            x = xsol

            for r in range(len(yn)):
                y_res = np.append(y_res, yn)

            y = y_res

        return [xsol, y_res]

    def execute(self):
        [ts, ys] = self.ABM4thOrder(self.y0, self.x, self.h)
        for index in ts:
            self.f.write(str(index) + ' ')
        self.f.write('\n')
        for index in ys:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        t = np.arange(0, 200, 1)
        for index in t:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        for i in t:
            y_math_res = 40 + 40 * np.e**(-(np.log(2) * i/20))
            self.f.write(str(y_math_res) + ' ')
        self.f.write('\n')

        return 0