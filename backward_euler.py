import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class BackwardEuler:
    def __init__(self, _h = 0.2, _x = np.array([1.0, 2.0]), _yinit = np.array([4.0])):
        self.h = _h
        self.x = _x
        self.yinit = np.array([4.0])

    def __del__(self):
        pass

    def main_BE(self, yinit, x_range, h):
        m = len(yinit)
        n = int((x_range[-1] - x_range[0])/h)

        x = x_range[0]
        y = yinit

        xsol = np.empty(0)
        xsol = np.append(xsol, x)

        global ysol
        ysol = np.empty(0)
        ysol = np.append(ysol, y)

        for i in range(n):
            # yprime = self.myFunc(x+h, y)/(1+h)
            yprime = mf.myFunc(x+h, y)/(1+h)

            for j in range(m):
                y[j] = y[j] + h*yprime[j]

            x += h
            xsol = np.append(xsol, x)

            for r in range(len(y)):
                ysol = np.append(ysol, y[r])  # Saves all new y's

        return [xsol, ysol]

    # def myFunc(self, x, y):
    #     dy = np.zeros((len(y)))
    #     dy[0] = 3*(1+x) - y[0]
    #     return dy

    def execute(self):
        [ts, ys] = self.main_BE(yinit=self.yinit, x_range=self.x, h=self.h)


        #--- Calculate the exact solution, for comparison ---#
        dt = int((self.x[-1] - self.x[0]) / self.h)
        t = [self.x[0]+i*self.h for i in range(dt+1)]
        yexact = []
        for i in range(dt+1):
            ye = 3 * t[i] + np.exp(1 - t[i])
            yexact.append(ye)

        return ts, ys, t, yexact


        # plt.plot(ts, ys, 'r')
        # plt.plot(t, yexact, 'b')
        # plt.xlim(self.x[0], self.x[1])
        # plt.legend(["Backward Euler method",
        #             "Exact solution"], loc=2)
        # plt.xlabel('x', fontsize=17)
        # plt.ylabel('y', fontsize=17)
        # plt.tight_layout()
        # plt.show()

# x = BackwardEuler()
# x.execute()