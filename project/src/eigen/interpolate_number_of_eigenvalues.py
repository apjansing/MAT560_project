from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np

xnew = np.linspace(218, 3488, endpoint=True)
x = np.array( [ 218, 872, 3488 ] ) #cells in mesh
y = np.array( [ 252, 938, 3618 ] ) #eigenvalues
print(x, y)
f = interp1d(x, y)

plt.plot(x, y, '-', xnew, f(xnew))
plt.plot(x, y, 'x', xnew, f(xnew))
plt.show()