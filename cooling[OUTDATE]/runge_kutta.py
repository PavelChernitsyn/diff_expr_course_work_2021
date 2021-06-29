
import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class Runge_Kutt:
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

    def RKF45(self, y0, x_range, h):
        y_len = len(y0)
        x_len = int((x_range[-1] - x_range[0])/h)
        
        x = x_range[0]
        y = y0
        
        self.f.write(str(x) + ' ')

        y_res = np.empty(0)
        y_res = np.append(y_res, y)

        for i in range(x_len):
            k1 = mf.myFunc(y, self.env_temp)

            yp2 = y + k1*(h/5)

            k2 = mf.myFunc(yp2, self.env_temp)

            yp3 = y + k1*(3*h/40) + k2*(9*h/40)

            k3 = mf.myFunc(yp3, self.env_temp)

            yp4 = y + k1*(3*h/10) - k2*(9*h/10) + k3*(6*h/5)

            k4 = mf.myFunc(yp4, self.env_temp)

            yp5 = y - k1*(11*h/54) + k2*(5*h/2) - k3*(70*h/27) + k4*(35*h/27)

            k5 = mf.myFunc(yp5, self.env_temp)

            yp6 = y + k1*(1631*h/55296) + k2*(175*h/512) + k3*(575*h/13824) + k4*(44275*h/110592) + k5*(253*h/4096)

            k6 = mf.myFunc(yp6, self.env_temp)

            for j in range(y_len):
                y[j] = y[j] + h*(37*k1[j]/378 + 250*k3[j]/621 + 125*k4[j]/594 + 512*k6[j]/1771)

            x = x + h
            self.f.write(str(x) + ' ')

            for r in range(len(y)):
                y_res = np.append(y_res, y[r])

        self.f.write('\n')  

        return y_res

    def execute(self):
        ys = self.RKF45(self.y0, self.x, self.h)
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