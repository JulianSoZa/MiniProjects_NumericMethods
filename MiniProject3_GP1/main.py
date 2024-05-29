import meshio
import numpy as np
from scipy.sparse import csr_matrix
import pyvista as pv
from modules import initialConditions
import json

### Se carga la malla
mesh = meshio.read("Antioquia.msh")

vertices = mesh.points
num_vertices = len(vertices)

triangles = mesh.cells_dict['triangle']

#

dataK = []
dataM = []

row = []
col = []

for tri in triangles:

    # Puntos de los triangulos
    PA = vertices[tri[0],:2] 
    PB = vertices[tri[1],:2]
    PC = vertices[tri[2],:2]
    
    # Matriz dde masa elemental
    Ae = (1/2)*np.cross(PB - PA, PC - PA)
    Me = (Ae/3)*np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    
    #Matriz de rigidez
    Be = np.array([[-1, 1, 0], [-1, 0, 1]])
    
    Ja = np.array([[PB[0]-PA[0], PC[0]-PA[0]], [PB[1]-PA[1], PC[1]-PA[1]]])
    
    J_inv = (1/np.linalg.det(Ja))*np.array([[PC[1]-PA[1], -PC[0]+PA[0]], [-PB[1]+PA[1], PB[0]-PA[0]]])
    
    Ke = Be.transpose()@J_inv.transpose()@Ja@Be*Ae
        
    ### Matrices
    for i in range(3):
        I = tri[i]
        
        for j in range(3):
            J = tri[j]

            if J == I:
                dataK.append(Ke[i][j])
                dataM.append(Me[i][j])
                row.append(I)
                col.append(J)

#indices = np.where(np.array(row) == 0)[0]

K = csr_matrix((dataK, (row, col)))
M = csr_matrix((dataM, (row, col)))

M_diag = M.diagonal()

delt = 0.00001
Di = 300
Ds = 100
beta = 300
gamma = 10

with open("municipios.json", 'r') as openfile:
    towns_data = json.load(openfile)
towns = ["Sopetran", "San Luis"] ### Lista de municipios que har´an parte de los suceptibles
i_towns = ["Entrerrios", "San Luis", "Sopetran"] ### Lista de municipios que har´an parte de los infectados
### El array points corresponde a los puntos de la malla ya creada
Sn, In = initialConditions.create_initial_conditions(towns, i_towns, vertices, towns_data, 1, alpha=0.6)

Is = []
Ss = []

Is.append(In)
Ss.append(Sn)

for i in range(500):
    bI = M@In + delt*(-Di*K@In + beta*M@Sn@In - gamma*M@In)
    bS = M@Sn + delt*(-Ds*K@Sn - beta*M@Sn@In)

    In = bI/M_diag
    Sn = bS/M_diag
    
    Is.append(In)
    Ss.append(Sn)


malla = meshio.read("Antioquia.msh")

malla.point_data['Infectados_Inicial'] = Is[0]
malla.point_data['Infectados_Final'] = Is[-1]

plotter = pv.Plotter(shape=(1, 2))
#plotter = pv.Plotter()

plotter.subplot(0, 0)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Inicial')

plotter.subplot(0, 1)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Final')
plotter.show()