
import matplotlib.pyplot as plt
import numpy as np


def feval(funcName, *args):
    return eval(funcName)(*args)


def HeunsMethod(func, yinit, x_range, h):
    m = len(yinit)
    n = int((x_range[-1] - x_range[0])/h)
    
    x = x_range[0]
    y = yinit
    
    # Solution arrays
    xsol = np.empty(0)
    xsol = np.append(xsol, x)

    ysol = np.empty(0)
    ysol = np.append(ysol, y)

    for i in range(n):
        y0prime = feval(func, x, y)

        k1 = y0prime * h

        ypredictor = y + k1

        y1prime = feval(func, x+h, ypredictor)

        for j in range(m):
            y[j] = y[j] + (h/2)*y0prime[j] + (h/2)*y1prime[j]

        x = x + h
        xsol = np.append(xsol, x)

        for r in range(len(y)):
            ysol = np.append(ysol, y[r])  # Saves all new y's

    return [xsol, ysol]


def myFunc(x, y):
    dy = np.zeros((len(y)))
    dy[0] = 3 * (1 + x) - y[0]
    return dy

# -----------------------

h = 0.2
x = np.array([1, 2])
yinit = np.array([4.0])


[ts, ys] = HeunsMethod('myFunc', yinit, x, h)


dt = int((x[-1] - x[0]) / h)
t = [x[0]+i*h for i in range(dt+1)]
yexact = []
for i in range(dt+1):
    ye = 3 * t[i] + np.exp(1 - t[i])
    yexact.append(ye)
    

plt.plot(ts, ys, 'r')
plt.plot(t, yexact, 'b')
plt.xlim(x[0], x[1])
plt.legend(["Heun's method", "Exact solution"], loc=2)
plt.xlabel('x', fontsize=17)
plt.ylabel('y', fontsize=17)
plt.tight_layout()
plt.show()

