
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class ForwardEuler:
    def __init__(self, _h = 0.2, _x = np.array([1.0, 2.0]), _yinit = np.array([4.0])):
        self.h = _h
        self.x = _x
        self.yinit = np.array([4.0])

    def __del__(self):
        pass

    # def feval(self, funcName, *args):
    #     return eval(funcName)(*args)


    def main_FE(self, yinit, x_range, h):
        m = len(yinit) # Number of ODEs
        n = int((x_range[-1] - x_range[0])/h) # Number of sub-intervals
        
        x = x_range[0] # Initializes variable x
        y = yinit # Initializes variable y
        
        xsol = np.empty(0) # Creates an empty array for x
        xsol = np.append(xsol, x) # Fills in the first element of xsol

        ysol = np.empty(0) # Creates an empty array for y
        ysol = np.append(ysol, y) # Fills in the initial conditions

        for i in range(n):
            # yprime = self.feval(func, x, y) # Evaluates dy/dx
            yprime = mf.myFunc(x, y)
            
            for j in range(m):
                y[j] = y[j] + h*yprime[j]
                
            x += h # Increase x-step
            xsol = np.append(xsol, x) # Saves it in the xsol array
            
            for r in range(len(y)):
                ysol = np.append(ysol, y[r]) # Saves all new y's 
                
        return [xsol, ysol]



    # def myFunc(self, x, y):
    #     dy = np.zeros((len(y)))
    #     dy[0] = 3*(1+x) - y[0]
    #     return dy


    def execute(self):
        [ts, ys] = self.main_FE(self.yinit, self.x, self.h)


        # Calculates the exact solution, for comparison
        dt = int((self.x[-1] - self.x[0]) / self.h)
        t = [self.x[0]+i*self.h for i in range(dt+1)]
        yexact = []
        for i in range(dt+1):
            ye = 3*t[i] + np.exp(1-t[i])
            yexact.append(ye)

        return ts, ys, t, yexact


#         plt.plot(ts, ys, 'r')
#         plt.plot(t, yexact, 'b')
#         plt.xlim(self.x[0], self.x[1])
#         plt.legend(["Forward Euler method", 
#                     "Exact solution"], loc=2)
#         plt.xlabel('x', fontsize=17)
#         plt.ylabel('y', fontsize=17)
#         plt.tight_layout()
#         plt.show()

# x = ForwardEuler()
# x.execute()