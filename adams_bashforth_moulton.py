import numpy as np
import matplotlib.pyplot as plt
import my_func as mf

class ABM:
    def __init__(self, _h = 0.2, _x = np.array([1.0, 2.0]), _yinit = np.array([4.0])):
        self.h = _h
        self.x = _x
        self.yinit = np.array([4.0])

    def __del__(self):
        pass

    # def feval(self, funcName, *args):
    #     return eval(funcName)(*args)


    def RungeKutta4thOrder(self, yinit, x, h):
        m = len(yinit)
        n = int((x[-1] - x[0]) / h)

        x = x[0]
        y = yinit

        xsol = np.empty((0))
        xsol = np.append(xsol, x)

        ysol = np.empty((0))
        ysol = np.append(ysol, y)

        for i in range(n):
            # k1 = self.feval(func, x, y)
            k1 = mf.myFunc(x, y)

            yp2 = y + k1*(h/2)

            # k2 = self.feval(func, x+h/2, yp2)
            k2 = mf.myFunc(x+h/2, yp2)

            yp3 = y + k2*(h/2)

            # k3 = self.feval(func, x+h/2, yp3)
            k3 = mf.myFunc(x+h/2, yp3)

            yp4 = y + k3*h

            # k4 = self.feval(func, x+h, yp4)
            k4 = mf.myFunc(x+h, yp4)

            for j in range(m):
                y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])

            x = x + h
            xsol = np.append(xsol, x)

            for r in range(len(y)):
                ysol = np.append(ysol, y[r]) 

        return [xsol, ysol]


    def ABM4thOrder(self, yinit, x, h):

        m = len(yinit)

        dx = int((x[-1] - x[0]) / h)

        xrk = [x[0] + k * h for k in range(dx + 1)]

        [xx, yy] = self.RungeKutta4thOrder(yinit, (xrk[0], xrk[3]), h)

        x = xx
        xsol = np.empty(0)
        xsol = np.append(xsol, x)

        y = yy
        yn = np.array([yy[0]])
        ysol = np.empty(0)
        ysol = np.append(ysol, y)

        for i in range(3, dx):
            x00 = x[i]; x11 = x[i-1]; x22 = x[i-2]; x33 = x[i-3]; xpp = x[i]+h

            y00 = np.array([y[i]])
            y11 = np.array([y[i - 1]])
            y22 = np.array([y[i - 2]])
            y33 = np.array([y[i - 3]])

            # y0prime = self.feval(func, x00, y00)
            # y1prime = self.feval(func, x11, y11)
            # y2prime = self.feval(func, x22, y22)
            # y3prime = self.feval(func, x33, y33)
            y0prime = mf.myFunc(x00, y00)
            y1prime = mf.myFunc(x11, y11)
            y2prime = mf.myFunc(x22, y22)
            y3prime = mf.myFunc(x33, y33)

            ypredictor = y00 + (h/24)*(55*y0prime - 59*y1prime + 37*y2prime - 9*y3prime)
            # ypp = self.feval(func, xpp, ypredictor)
            ypp = mf.myFunc(xpp, ypredictor)

            for j in range(m):
                yn[j] = y00[j] + (h/24)*(9*ypp[j] + 19*y0prime[j] - 5*y1prime[j] + y2prime[j])

            xs = x[i] + h
            xsol = np.append(xsol, xs)

            x = xsol

            for r in range(len(yn)):
                ysol = np.append(ysol, yn)

            y = ysol

        return [xsol, ysol]


    # def myFunc(self, x, y):
    #     dy = np.zeros((len(y)))
    #     dy[0] = 3*(1+x) - y[0]
    #     return dy

    def execute(self):
        [ts, ys] = self.ABM4thOrder(self.yinit, self.x, self.h)


        dt = int((self.x[-1]-self.x[0])/self.h)
        t = [self.x[0]+i*self.h for i in range(dt+1)]
        yexact = []
        for i in range(dt+1):
            ye = 3 * t[i] + np.exp(1 - t[i])
            yexact.append(ye)

        # plt.plot(ts, ys, 'r')
        # plt.plot(t, yexact, 'b')
        # #plt.xlim(x[0], x[1])
        # plt.legend(["ABM 4th Order ", "Exact solution"], loc=2)
        # plt.xlabel('x', fontsize=17)
        # plt.ylabel('y', fontsize=17)
        # plt.tight_layout()
        # plt.show()

        return ts, ys, t, yexact

    # diff = ys - yexact
    # print("Maximum difference =", np.max(abs(diff)))

    # plt.plot(ts, ys, 'r')
    # plt.plot(t, yexact, 'b')
    # plt.xlim(x[0], x[1])
    # plt.legend(["ABM 4th Order ", "Exact solution"], loc=2)
    # plt.xlabel('x', fontsize=17)
    # plt.ylabel('y', fontsize=17)
    # plt.tight_layout()
    # plt.show()

# x = ABM()
# x.execute()