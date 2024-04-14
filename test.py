import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

nx = 6 #multiplo de 4
ny = 6

dx = 8/nx
dy = 6/ny

vt =  2
vb =  -1

data = []
row = []
col = []

ia = int(2/dx)
ib = int(6/dx)

jt = int(4/dy)
jb = int(2/dy)

nk = (nx + 1)*(ny + 1)
b = np.zeros(nk)

alpha = dy/dx
for k in range(nk):
    i = k%(nx+1)
    j = int(k/(nx+1))
    
    ### Fronteras superior o inferior
    if j==0 or j==ny:
        data.append(1)
        row.append(k)
        col.append(k)
        continue
    
    ### Fronteras laterales
    if i==0 or i==nx:
        data.append(1)
        row.append(k)
        col.append(k)
        continue
    
    ### Plato superior del capacitor
    if j==jt and (ia <= i <= ib):
        data.append(1)
        row.append(k)
        col.append(k)
        b[k] = vt
        continue

    ### Plato inferior del capacitor
    if j==jb and (ia <= i <= ib):
        data.append(1)
        row.append(k)
        col.append(k)
        b[k] = vb
        continue
    
    ### Derecha o izquierda
    J = j
    for di in [-1, 1]:
        I = i + di
        K = (nx+1)*J + I
        data.append( alpha**2)
        row.append(k)
        col.append(K)
        
    ### Arriba o abajo
    I = i
    for dj in [-1, 1]:
        J = j + dj
        K = (nx+1)*J + I
        data.append(1)
        row.append(k)
        col.append(K)
        
    ### Central
    # val = -2*( dx**2 + dy**2 )
    val = -2*( 1 + alpha**2 )
    data.append( val )
    row.append(k)
    col.append(k)
    
A = csr_matrix((data, (row, col)), shape=(nk, nk))
V = spsolve(A,b)

# Convertir la matriz densa a formato de diccionario, obtener un array de coordenadas xy y un array de valores
mtrx_dict = A.todok()
xy = np.array(list(mtrx_dict.keys()))
vals = np.array(list(mtrx_dict.values()))

fig, ax = plt.subplots()
# Crear un gráfico de dispersión
ax.scatter(xy[:,0], xy[:,1], s=20, c=vals)
# Establecer los ticks mayores cada 2 unidades y los ticks menores cada 0.5 unidades
ax.set_xticks(np.arange(0, nk, 1))
ax.set_yticks(np.arange(0, nk, 1))

# Habilitar el grid para los ticks mayores y menores
ax.grid(which='both')

# Establecer el estilo del grid para los ticks mayores
ax.grid(which='major', linewidth=1)

# Establecer el estilo del grid para los ticks menores
ax.grid(which='minor', linewidth=0.2)
plt.show()