import numpy as np
import meshio
import pyvista as pv
import matplotlib.pyplot as plt

def plotter_analysis():

    malla = meshio.read("MiniProject3_GP1/data/meshes/VirusT2_Antioquia.vtk")

    plotter = pv.Plotter(shape=(1, 2))

    plotter.subplot(0, 0)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_dia_14')

    plotter.subplot(0, 1)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_dia_14')

    plotter.show()

    InfT = np.load('MiniProject3_GP1/data/variables/infectados_instantes.npy')
    SusT = np.load('MiniProject3_GP1/data/variables/Susceptibles_instantes.npy')

    """triangles = malla.cells_dict['triangle']

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
    np.save("Susceptibles_instantes", SusT)"""

    t_span = np.linspace(0, 2, 3)

    plt.plot(t_span, InfT, label='Infectados')
    plt.plot(t_span, SusT, label='Susceptibles')
    plt.legend()

    plt.show()

    plotter = pv.Plotter()
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Final')
    plotter.view_xy()
    plotter.open_gif('MiniProject3_GP1/data/gifs/Infectados.gif')
    print(dir(plotter))
    for i in range(3):
        I = malla.point_data[f'Infectados_dia_{int(i*7)}']
        Imin = np.min(I)
        Imax = np.max(I)
        I = np.where(I<0.8E-6, np.nan, I)
        plotter.update_scalars(I, render=False)
        plotter.add_title(f"Infectados día {int(i*7)}", font_size=14)
        plotter.update_scalar_bar_range(clim = [Imin, Imax])
        plotter.view_xy()
        plotter.write_frame()
        
        
    plotter.close()

    plotter = pv.Plotter()
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_Final')
    plotter.view_xy()
    plotter.open_gif('MiniProject3_GP1/data/gifs/Susceptibles.gif')
    for i in range(3):
        S = malla.point_data[f'Susceptibles_dia_{int(i*7)}']
        Smin = np.min(S)
        Smax = np.max(S)
        S = np.where(S<1.5E-5, np.nan, S)
        plotter.update_scalars(S, render=False)
        plotter.add_title(f"Susceptibles día {int(i*7)}", font_size=20, color='k')
        plotter.update_scalar_bar_range(clim = [Smin, Smax])
        plotter.view_xy()
        plotter.write_frame()
        
    plotter.close()
    
if __name__ == "__main__":
    plotter_analysis()