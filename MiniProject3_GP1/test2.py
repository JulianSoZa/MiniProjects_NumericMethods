import numpy as np
import meshio
import pyvista as pv
import matplotlib.pyplot as plt

malla = meshio.read("VirusT2_Antioquia.vtk")

plotter = pv.Plotter(shape=(1, 2))

plotter.subplot(0, 0)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_dia_343')

plotter.subplot(0, 1)
plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_dia_343')

triangles = malla.cells_dict['triangle']

InfT = []
SusT = []

Aes = np.load(f"Aes.npy")

for i, j in zip(range(50), range(50)):
    infectados = 0
    susceptibles = 0

    for tri in enumerate(triangles):

        I_PA = malla.point_data[f'Infectados_dia_{i*7}'][tri[1][0]] 
        I_PB = malla.point_data[f'Infectados_dia_{i*7}'][tri[1][1]]
        I_PC = malla.point_data[f'Infectados_dia_{i*7}'][tri[1][2]]
        
        infectados += (1/3)*(I_PA + I_PB + I_PC)*Aes[tri[0]]
        
        S_PA = malla.point_data[f'Susceptibles_dia_{i*7}'][tri[1][0]] 
        S_PB = malla.point_data[f'Susceptibles_dia_{i*7}'][tri[1][1]]
        S_PC = malla.point_data[f'Susceptibles_dia_{i*7}'][tri[1][2]]
        
        susceptibles += (1/3)*(S_PA + S_PB + S_PC)*Aes[tri[0]]
        
    InfT.append(infectados)
    SusT.append(susceptibles)
    
    print('Infectados: ', infectados)
    print('susceptibles: ', susceptibles)
    
np.save("Infectados_instantes", InfT)
np.save("Susceptibles_instantes", SusT)

t_span = np.linspace(0, 343, 50)

plt.plot(t_span, InfT, label='Infectados')
plt.plot(t_span, SusT, label='Susceptibles')
plt.legend()

plt.show()

