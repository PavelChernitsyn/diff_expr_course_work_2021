import numpy as np

def myFunc(x, mu, l):
    dy = 4.9 * x**2 * (mu + 1) / l - 9.8 * mu * x
    # dy = 5.self.time9 * (x)**2 - 0.98 * (x)
    return dy