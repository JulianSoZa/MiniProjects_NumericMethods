import numpy as np
import meshio
import pyvista as pv

matriz = np.array([[1,2,3],[4,5,6]])
print(matriz)
np.save('matriz00.npy',matriz)

import os

# Reemplaza 'nombre_del_archivo' con el nombre de tu archivo
nombre_del_archivo = 'matriz00.npy'

if os.path.isfile(nombre_del_archivo):
    print('El archivo existe.')
else:
    print('El archivo no existe')

readtest = np.load("Matriz_de_masa.npy")

print (readtest)

