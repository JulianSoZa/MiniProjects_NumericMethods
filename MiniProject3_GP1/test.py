import meshio
import gmsh
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator as Interpolator
from scipy.interpolate import NearestNDInterpolator as NRInterpolator
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


name = 'Antioquia'
mesh = meshio.read(f"{name}.msh")
niveles = np.load('niveles.npy')
data = []
row = []
col = []

vertices = mesh.points
num_vertices = len(vertices)
triangles = mesh.cells_dict['triangle']



print('Vamos k ')

x = niveles[:,0]
y = niveles[:,1]
z = niveles[:,2]

X = np.array(vertices[:,0], dtype=np.float32)
Y = np.array(vertices[:,1], dtype=np.float32)


print('Vamos kk ')
interp = Interpolator(list(zip(x,y)),z)
Z = interp(X,Y)

plt.pcolormesh(X, Y, Z, shading='auto')
plt.plot(x, y, "ok", label="input point")
plt.legend()
plt.colorbar()
plt.axis("equal")
plt.show()

'''for k in vertices:
    data.append(1)
    row.append(k)
    col.append(k)

for tri in triangles:
    #Puntos de cada triangulo (tri)
    PA = vertices[tri[0],:2]
    PB = vertices[tri[1],:2]
    PC = vertices[tri[2],:2]

    #Saca el jaconiano (?)
    G = np.matrix([PB-PA,PC-PA]).transpose()
    det_G = np.linalg.det(G)

    B = [ [PB[1]-PC[1], PC[1]-PA[1], PA[1]-PB[1]],
          [PC[0]-PB[0], PA[0]-PC[0], PB[0]-PA[0]]]
    B = np.matrix(B)

    K_local = 0.5*B.transpose()*B/det_G
    for i in range(3):
        I = tri[i]
        if I in vertices:
            continue
        for j in range(3):
            J = tri[j]

            ### Matrix K
            data.append( -K_local[i,j] )
            row.append(I)
            col.append(J)

K = csr_matrix((data, (row, col)), shape=(num_vertices, num_vertices))
b = np.zeros(num_vertices)

print(K)'''

#------------------------------------------------------------------------------------------------------------

'''
plt.figure()

plt.triplot(vertices[:,0], vertices[:,1], triangles, linewidth=1)
#plt.plot(vertices[:,0], vertices[:,1], 'o')

plt.grid()
plt.axis('equal')
plt.title("Malla")
plt.xlabel("x"); plt.ylabel("y")
#plt.show()
'''