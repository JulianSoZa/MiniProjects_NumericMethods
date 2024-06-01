import meshio
import numpy as np
from scipy.sparse import csr_matrix
import pyvista as pv
from modules import initialConditions
import json

### Se carga la malla
mesh = meshio.read("AntioquiaConNiveles.vtk")

vertices = mesh.points
num_vertices = len(vertices)

triangles = mesh.cells_dict['triangle']

#

dataK = []
dataM = []

rowM = []
colM = []

rowK = []
colK = []

Aes = []

dh = np.inf

for tri in triangles:

    # Puntos de los triangulos
    PA = vertices[tri[0],:2] 
    PB = vertices[tri[1],:2]
    PC = vertices[tri[2],:2]
    
    lado_min = min(np.array([np.linalg.norm(PB - PA), np.linalg.norm(PC - PA), np.linalg.norm(PC - PB)]))
    
    if lado_min<dh:
        dh = lado_min
    
    # Matriz dde masa elemental
    Ae = (1/2)*np.cross(PB - PA, PC - PA)
    Aes.append(Ae)
    Me = (Ae/6)*np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    
    #Matriz de rigidez
    G = np.matrix([PB-PA,PC-PA]).transpose()
    
    B = np.matrix([[PB[1]-PC[1], PC[1]-PA[1], PA[1]-PB[1]],
                   [PC[0]-PB[0], PA[0]-PC[0], PB[0]-PA[0]]])

    Ke = 0.5*B.transpose()*B/np.linalg.det(G)

    ### Matrices
    for i in range(3):
        I = tri[i]
        
        for j in range(3):
            J = tri[j]
            
            if I == J:
                dataM.append(Me[i,j])
                rowM.append(I)
                colM.append(J)
                
            dataK.append(Ke[i,j])
            rowK.append(I)
            colK.append(J)

K = csr_matrix((dataK, (rowK, colK)))
M = csr_matrix((dataM, (rowM, colM)))

M_diag = M.diagonal()

Di = 1.5E6
Ds = 1E6
beta = 0.9
gamma = 0.2

dt_cri = [2*(dh**2)/(4*Di+1.5*gamma*dh), (dh**2)/(2*Ds)]

dt = np.min(dt_cri)

print('dh: ', dh)
print('dt: ', dt)

dt = 0.002

T = 3
Nt = int(T/dt)
t_save = np.round(7/dt)

with open("municipios.json", 'r') as openfile:
    towns_data = json.load(openfile)

towns = ["Entrerrios"] ### Lista de municipios que har´an parte de los suceptibles
i_towns = ["Entrerrios"] ### Lista de municipios que har´an parte de los infectados
### El array points corresponde a los puntos de la malla ya creada
Sn, In = initialConditions.create_initial_conditions(towns, i_towns, vertices, towns_data, 1, alpha=0.6)

Is = []
Ss = []

Is.append(In)
Ss.append(Sn)

Smax = np.max(Sn)

dt_span = np.repeat(dt, Nt)
#dt_span = np.concatenate((dt_span, np.repeat(dt*1.2, 7/(dt*1.2))))

for i in enumerate(dt_span, 1):
    bI = M*In + i[1]*(-Di*K*In + beta*M*Sn*In - gamma*M*In)
    bS = M*Sn + i[1]*(-Ds*K*Sn - beta*M*Sn*In)
    
    In = bI/M_diag
    Sn = bS/M_diag
    
    if np.max(Sn)>Smax:
        print('Posible inestabilidad')
    
    """if i%t_save==0:
        Is.append(In)
        Ss.append(Sn)"""
    
    if i[0]%(len(dt_span))==0:
        Is.append(In)
        Ss.append(Sn)
        
print(len(Is))

malla = meshio.read("AntioquiaConNiveles.vtk")

malla.point_data['Infectados_Inicial'] = Is[0]
malla.point_data['Infectados_Final'] = Is[-1]

plotter = pv.Plotter(shape=(1, 2))

plotter.subplot(0, 0)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Inicial')

plotter.subplot(0, 1)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Final')

plotter.show()

malla.point_data['Susceptibles_Inicial'] = Ss[0]
malla.point_data['Susceptibles_Final'] = Ss[-1]

plotter = pv.Plotter(shape=(1, 2))

plotter.subplot(0, 0)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_Inicial')

plotter.subplot(0, 1)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_Final')

plotter.show()

for i, j in zip(range(len(Is)), range(len(Ss))):
    infectados = 0
    susceptibles = 0

    for tri in enumerate(triangles):

        I_PA = Is[i][tri[1][0]] 
        I_PB = Is[i][tri[1][1]]
        I_PC = Is[i][tri[1][2]]
        
        infectados += (1/3)*(I_PA + I_PB + I_PC)*Aes[tri[0]]
        
        S_PA = Ss[j][tri[1][0]] 
        S_PB = Ss[j][tri[1][1]]
        S_PC = Ss[j][tri[1][2]]
        
        susceptibles += (1/3)*(S_PA + S_PB + S_PC)*Aes[tri[0]]
        
    print('Infectados: ', infectados)
    print('susceptibles: ', susceptibles)