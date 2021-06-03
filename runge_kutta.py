
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class Runge_Kutt:
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


    def RKF45(self, yinit, x_range, h):
        m = len(yinit)
        n = int((x_range[-1] - x_range[0])/h)
        
        x = x_range[0]
        y = yinit
        
        # xsol = np.empty(0)
        # xsol = np.append(xsol, x)
        self.f.write(str(x) + ' ')

        ysol = np.empty(0)
        ysol = np.append(ysol, y)

        for i in range(n):
            # k1 = self.feval(func, x, y)
            k1 = mf.myFunc(x, y)

            yp2 = y + k1*(h/5)

            # k2 = self.feval(func, x+h/5, yp2)
            k2 = mf.myFunc(x+h/5, yp2)

            yp3 = y + k1*(3*h/40) + k2*(9*h/40)

            # k3 = self.feval(func, x+(3*h/10), yp3)
            k3 = mf.myFunc(x+(3*h/10), yp3)

            yp4 = y + k1*(3*h/10) - k2*(9*h/10) + k3*(6*h/5)

            # k4 = self.feval(func, x+(3*h/5), yp4)
            k4 = mf.myFunc(x+(3*h/5), yp4)

            yp5 = y - k1*(11*h/54) + k2*(5*h/2) - k3*(70*h/27) + k4*(35*h/27)

            # k5 = self.feval(func, x+h, yp5)
            k5 = mf.myFunc(x+h, yp5)

            yp6 = y + k1*(1631*h/55296) + k2*(175*h/512) + k3*(575*h/13824) + k4*(44275*h/110592) + k5*(253*h/4096)

            # k6 = self.feval(func, x+(7*h/8), yp6)
            k6 = mf.myFunc(x+(7*h/8), yp6)


            for j in range(m):
                y[j] = y[j] + h*(37*k1[j]/378 + 250*k3[j]/621 + 125*k4[j]/594 + 512*k6[j]/1771)

            x = x + h
            # xsol = np.append(xsol, x)
            self.f.write(str(x) + ' ')

            for r in range(len(y)):
                ysol = np.append(ysol, y[r])

        self.f.write('\n')  

        return ysol


    # def myFunc(self, x, y):
    #     dy = np.zeros((len(y)))
    #     dy[0] = 3*(1+x) - y[0]
    #     return dy

    def execute(self):
        ys = self.RKF45(self.yinit, self.x, self.h)
        for index in ys:
            self.f.write(str(index) + ' ')
        self.f.write('\n')

        dt = int((self.x[-1]-self.x[0])/self.h)
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


#         y_diff = ys - yexact
#         print("max diff =", np.max(abs(y_diff)))


#         plt.plot(ts, ys, 'rs')
#         plt.plot(t, yexact, 'b')
#         plt.xlim(self.x[0], self.x[1])
#         plt.legend(["RKF45", "Exact solution"], loc=1)
#         plt.xlabel('x', fontsize=17)
#         plt.ylabel('y', fontsize=17)
#         plt.tight_layout()
#         plt.show()

#         # Uncomment the following to print the figure:
#         #plt.savefig('Fig_ex2_RK4_h0p1.png', dpi=600)


# x = Runge_Kutt()
# x.execute()
