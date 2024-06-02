import meshio
import numpy as np
from scipy.sparse import csr_matrix, save_npz, load_npz
import pyvista as pv
from modules import initialConditions
import json
import matplotlib.pyplot as plt

### Se carga la malla
malla = meshio.read("MiniProject3_GP1\data\meshes\AntioquiaConNiveles.vtk")

vertices = malla.points
num_vertices = len(vertices)

triangles = malla.cells_dict['triangle']

Ma = 'Matriz_de_masa'
Mr = 'Matriz_de_Rigidez'
d_h = 'dh'
A_es = 'Aes'

try: 
    print('Inicia la lectura\n')
    M = load_npz(f"MiniProject3_GP1/data/variables/{Ma}.npz")
    K = load_npz(f"MiniProject3_GP1/data/variables/{Mr}.npz")
    dh = np.load(f"MiniProject3_GP1/data/variables/{d_h}.npy")
    Aes = np.load(f"MiniProject3_GP1/data/variables/{A_es}.npy")
    print('Se leyó correctamente\n')
except:
    dh = np.inf
    Aes = []
    dataK = []
    dataM = []

    rowM = []
    colM = []

    rowK = []
    colK = []
    for tri in triangles:
        print('No ne leyó correctamente\n')
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
    save_npz(f'MiniProject3_GP1/data/variables/{Mr}.npz',K)
    save_npz(f'MiniProject3_GP1/data/variables/{Ma}.npz',M)
    np.save(f"MiniProject3_GP1/data/variables/{d_h}", dh)
    np.save(f"MiniProject3_GP1/data/variables/{A_es}", Aes)
    print('Se almacenaron las variables')

M_diag = M.diagonal()

Di = 0.6E6
Ds = 0.5E6
beta = 10
gamma = 0.002

dt_cri = [2*(dh**2)/(4*Di+1.5*gamma*dh), (dh**2)/(2*Ds)]

dt = np.min(dt_cri)

print('dh: ', dh)
print('dt: ', dt)

dt = 0.006

T = 7*2
Nt = int(T/dt)
dT = 7
t_save = np.round(dT/dt)

with open("MiniProject3_GP1/data/municipios.json", 'r') as openfile:
    towns_data = json.load(openfile)

towns = ['Guarne', 'Anza', 'Concordia', 'Heliconia', 'Girardota', 'Sabaneta', "Envigado", "Caldas", 
         "Bello", 'Entrerrios', 'Sopetran', 'Rionegro', 'San Roque', 'Don Matias', 'Barbosa', 'Cocorna',
         'Liborina', 'Belmira', 'Caicedo', 'Urrao', 'Guadalupe', 'Buritica', 'Giraldo', 'Abriaqui', 'Ebejico',
         'San Pedro De Los Milagros', 'Betulia', 'Titiribi', 'El Penol', 'Guatape', 'San Carlos', 'San Luis',
         'Yolombo', 'San Vicente', 'Itagui', 'Apartado', 'Turbo', 'Caucasia', 'Chigorodo', 'Copacabana', 
         'La Estrella', 'Necocli', 'La Ceja', 'Carepa', 'El Bagre', 'Puerto Berrio', 'Yarumal'] ### Lista de municipios que har´an parte de los suceptibles

i_towns = ["Medellin", 'Concepcion', 'Salgar', 'Armenia', 'Granada', 'Marinilla', 'Olaya', 'Gomez Plata', 'Santa Fe De Antioquia'] ### Lista de municipios que har´an parte de los infectados
### El array points corresponde a los puntos de la malla ya creada
Sn, In = initialConditions.create_initial_conditions(towns, i_towns, vertices, towns_data, 1, alpha=0.6)

malla.point_data[f'Infectados_dia_0'] = In
malla.point_data[f'Susceptibles_dia_0'] = Sn

Is = []
Ss = []

Is.append(In)
Ss.append(Sn)

Smax = np.max(Sn)

dt_span = np.repeat(dt, Nt+1)
#dt_span = np.concatenate((dt_span, np.repeat(dt*1.2, 7/(dt*1.2))))

for i in enumerate(dt_span, 1):
    bI = M*In + i[1]*(-Di*K*In + beta*M*Sn*In - gamma*M*In)
    bS = M*Sn + i[1]*(-Ds*K*Sn - beta*M*Sn*In)
    
    In = bI/M_diag
    Sn = bS/M_diag
    
    if np.max(Sn)>Smax:
        print('Posible inestabilidad')
    
    if (i[0]%t_save)==0:
        Is.append(In)
        Ss.append(Sn)
        
        print(f'Semana {int(i[0]*(dt)/7)} - Dia {int(i[0]*(dt))}')
        
        malla.point_data[f'Infectados_dia_{int(i[0]*(dt))}'] = In
        malla.point_data[f'Susceptibles_dia_{int(i[0]*(dt))}'] = Sn
    
    """if i[0]%(len(dt_span))==0:
        Is.append(In)
        Ss.append(Sn)"""
        
print(len(Is))

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

meshio.write('MiniProject3_GP1/data/meshes/VirusT2_Antioquia.vtk', malla)

InfT = []
SusT = []

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
        
    InfT.append(infectados)
    SusT.append(susceptibles)
    
    print('Infectados: ', infectados)
    print('susceptibles: ', susceptibles)
    
np.save("MiniProject3_GP1/data/variables/Infectados_instantes", InfT)
np.save("MiniProject3_GP1/data/variables/Susceptibles_instantes", SusT)

t_spanP = np.linspace(0, 13, 14)

plt.plot(t_spanP, InfT, label='Infectados')
plt.plot(t_spanP, SusT, label='Susceptibles')
plt.legend()

plt.show()