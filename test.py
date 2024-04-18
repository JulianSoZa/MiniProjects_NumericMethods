import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

# Crear algunos datos en una malla polar
r = np.linspace(0, 1, 4)
print(r)
theta = np.linspace(0, 2*np.pi, 100)
r, theta = np.meshgrid(r, theta)
z = np.sin(5*r) * np.cos(5*theta)

# Crear una malla rectangular para la interpolación
ri = np.array([0.5])
ti = np.linspace(0, 2*np.pi, 9)
xi, yi = np.meshgrid(ri, ti)

print('xi\n',xi)
print('yi\n',yi)

# Realizar la interpolación
zi = griddata((r.flatten(), theta.flatten()), z.flatten(), (xi, yi), method='cubic')

print(zi)

# Crear una figura y un conjunto de ejes
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

# Usar pcolormesh para graficar la malla 2D
c = ax.pcolormesh(theta, r, z, cmap='viridis')

# Agregar una barra de colores
fig.colorbar(c, ax=ax)

# Mostrar la gráfica
plt.show()