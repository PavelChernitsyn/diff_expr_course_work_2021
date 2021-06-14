import matplotlib.pyplot as plt
import numpy as np
import my_func as mf

class ForwardEuler:
    def __init__(self, _h = 0.2, _x = np.array([1.0, 200.0]), _yinit = np.array([80.0])):
        self.h = _h
        self.x = _x
        self.y0 = np.array([80.0])
        self.f = open('tmp.txt', 'w')

    def __del__(self):
        self.f.close()
        pass

    def main_FE(self, y0, x_range, h):
        y_len = len(y0)
        x_len = int((x_range[-1] - x_range[0])/h)
        
        x = x_range[0]
        y = y0
        
        self.f.write(str(x) + ' ')

        y_res = np.empty(0)
        y_res = np.append(y_res, y)

        for i in range(x_len):
            F_x_t = mf.myFunc(y)
            
            for j in range(y_len):
                y[j] = y[j] + h*F_x_t[j]
                
            x += h
            self.f.write(str(x) + ' ')

            for r in range(len(y)):
                y_res = np.append(y_res, y[r])

        self.f.write('\n')

        return y_res

    def execute(self):
        ys = self.main_FE(self.y0, self.x, self.h)
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