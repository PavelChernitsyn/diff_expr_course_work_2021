
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf


class Heun:
    def __init__(self, _h = 0.2, _x = np.array([1.0, 2.0]), _yinit = np.array([4.0])):
        self.h = _h
        self.x = _x
        self.yinit = np.array([4.0])
        self.f = open('tmp.txt', 'w')
    
    def __del__(self):
        self.f.close()
        pass

    # def feval(self, funcName, *args):
    #     return eval(funcName)(*args)


    def main_H(self, yinit, x_range, h):
        m = len(yinit)
        n = int((x_range[-1] - x_range[0])/h)
        
        x = x_range[0]
        y = yinit
        
        # Solution arrays
        # xsol = np.empty(0)
        # xsol = np.append(xsol, x)
        self.f.write(str(x) + ' ')
        
        ysol = np.empty(0)
        ysol = np.append(ysol, y)

        for i in range(n):
            # y0prime = self.feval(func, x, y)
            y0prime = mf.myFunc(x, y)

            k1 = y0prime * h

            ypredictor = y + k1

            # y1prime = self.feval(func, x+h, ypredictor)
            y1prime = mf.myFunc(x+h, ypredictor)

            for j in range(m):
                y[j] = y[j] + (h/2)*y0prime[j] + (h/2)*y1prime[j]

            x = x + h
            # xsol = np.append(xsol, x)
            self.f.write(str(x) + ' ')

            for r in range(len(y)):
                ysol = np.append(ysol, y[r])  # Saves all new y's

        self.f.write('\n')
        
        return ysol


    # def myFunc(seld, x, y):
    #     dy = np.zeros((len(y)))
    #     dy[0] = 3 * (1 + x) - y[0]
    #     return dy

    # -----------------------

    def execute(self):
        ys = self.main_H(self.yinit, self.x, self.h)
        for index in ys:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        dt = int((self.x[-1] - self.x[0]) / self.h)
        t = [self.x[0]+i*self.h for i in range(dt+1)]
        for index in t:
            self.f.write(str(index) + ' ')
        self.f.write('\n')
        yexact = []
        for i in range(dt+1):
            ye = 3 * t[i] + np.exp(1 - t[i])
            yexact.append(ye)
            self.f.write(str(ye) + ' ')
        self.f.write('\n')

        # return ts, ys, t, yexact
        return 0
            

#         plt.plot(ts, ys, 'r')
#         plt.plot(t, yexact, 'b')
#         plt.xlim(self.x[0], self.x[1])
#         plt.legend(["Heun's method", "Exact solution"], loc=2)
#         plt.xlabel('x', fontsize=17)
#         plt.ylabel('y', fontsize=17)
#         plt.tight_layout()
#         plt.show()

# x = Heun()
# x.execute()