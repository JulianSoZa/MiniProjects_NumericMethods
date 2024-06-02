import numpy as np
import meshio
import pyvista as pv
import matplotlib.pyplot as plt

def plotter_analysis(t_span, dia, dias):
    malla = meshio.read("MiniProject3_GP1/data/meshes/VirusT2_Antioquia.vtk")

    plotter = pv.Plotter(shape=(1, 2))
    plotter.subplot(0, 0)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Inicial')
    plotter.subplot(0, 1)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Final')
    plotter.show()
    
    plotter = pv.Plotter(shape=(1, 2))
    plotter.subplot(0, 0)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_Inicial')
    plotter.subplot(0, 1)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_Final')
    plotter.show()
    
    plotter = pv.Plotter(shape=(1, 2))
    plotter.subplot(0, 0)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars=f'Infectados_dia_{dia}')
    plotter.subplot(0, 1)
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars=f'Susceptibles_dia_{dia}')
    plotter.show()

    InfT = np.load('MiniProject3_GP1/data/variables/infectados_instantes.npy')
    SusT = np.load('MiniProject3_GP1/data/variables/Susceptibles_instantes.npy')

    plt.plot(t_span, InfT, label='Infectados')
    plt.plot(t_span, SusT, label='Susceptibles')
    plt.legend()

    plt.show()

    plotter = pv.Plotter()
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Infectados_Final')
    plotter.view_xy()
    plotter.open_gif('MiniProject3_GP1/data/gifs/Infectados.gif')

    for i in dias:
        I = malla.point_data[f'Infectados_dia_{int(i)}']
        Imin = np.min(I)
        Imax = np.max(I)
        I = np.where(I<0.8E-6, np.nan, I)
        plotter.update_scalars(I, render=False)
        plotter.add_title(f"Infectados día {int(i)}", font_size=14)
        plotter.update_scalar_bar_range(clim = [Imin, Imax])
        plotter.view_xy()
        plotter.write_frame()
        
    plotter.close()

    plotter = pv.Plotter()
    plotter.add_mesh(malla, show_edges=False, cmap='jet', scalars='Susceptibles_Final')
    plotter.view_xy()
    plotter.open_gif('MiniProject3_GP1/data/gifs/Susceptibles.gif')
    
    for i in dias:
        S = malla.point_data[f'Susceptibles_dia_{int(i)}']
        Smin = np.min(S)
        Smax = np.max(S)
        S = np.where(S<1.5E-5, np.nan, S)
        plotter.update_scalars(S, render=False)
        plotter.add_title(f"Susceptibles día {int(i)}", font_size=20, color='k')
        plotter.update_scalar_bar_range(clim = [Smin, Smax])
        plotter.view_xy()
        plotter.write_frame()
        
    plotter.close()
    
if __name__ == "__main__":
    try: 
        print('Inicia la lectura\n')
        malla = meshio.read("MiniProject3_GP1/data/meshes/VirusT2_Antioquia.vtk")

        T = 7*9
        dT = 7
        
        t_span = np.linspace(0, T, int(T/dT)+1)
        
        dia = '7'
        
        dias = np.array([0])
        dias = np.append(dias, np.linspace(dT, T, int(T/dT)))
        print(dias)
        
        plotter_analysis(t_span, dia, dias)
        
        print('Se leyó correctamente\n')
    except:
        print('Por favor, asegurese que este ingresando los valores y la malla correcta \n')