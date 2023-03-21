
# A exp[−(t−t0)^2/(2σt^2)] cos [Ω(t − t0)]

import math
import numpy as np
import matplotlib.pyplot as plt

A = 2.0
t0 = 0
omega = 10.0
sigma = 4.0

t = np.arange(-20,20,0.1)
y = A*np.exp(-(t-t0)**2/(2*sigma**2)) * np.cos(omega*(t-t0))
plt.plot(t,y)
plt.show()

