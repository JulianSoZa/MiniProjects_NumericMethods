
from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator
import matplotlib.pyplot as plt
import meshio
import numpy as np
import matplotlib.tri as tri

# Función para interpolación de vecino más cercano
def nr_interpolator_fun(x, y, values):
    nearest_interp = niveles[:,0]
    return nearest_interp(x, y)

# Función para interpolación lineal
def interpolator_fun(x, y, values):
    interp = LinearNDInterpolator(list(zip(x, y)), values)
    return interp(x, y)



niveles = np.load('niveles.npy')

name = 'Antioquia'
mesh = meshio.read(f"{name}.msh")
points = mesh.points

X = points[:,0]
Y = points[:,1]

interp = LinearNDInterpolator(list(zip(niveles[:,0], niveles[:,1])), niveles[:,2])
Z = interp(X, Y)

plt.pcolormesh(X, Y, Z, shading='auto')
#plt.plot(, y, "ok", label="input point")
plt.legend()
plt.colorbar()
plt.axis("equal")
plt.show()
