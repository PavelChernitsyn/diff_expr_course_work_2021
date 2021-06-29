import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class ForwardEuler:
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

    def main_FE(self):
        y_len = len(self.y0)
        x_len = int((self.x[-1] - self.x[0])/self.h)
        
        x = self.x[0]
        y = self.y0
        
        self.f.write(str(x) + ' ')

        y_res = np.empty(0)
        y_res = np.append(y_res, y)

        for i in range(x_len):
            F_x_t = mf.myFunc(y, self.env_temp)
            
            for j in range(y_len):
                y[j] = y[j] + self.h*F_x_t[j]
                
            x += self.h
            self.f.write(str(x) + ' ')

            for r in range(len(y)):
                y_res = np.append(y_res, y[r])

        self.f.write('\n')

        return y_res

    def execute(self):
        ys = self.main_FE()
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