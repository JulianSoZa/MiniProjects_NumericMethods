import pyvista as pv

mesh = pv.read('Malla_Campo_Polares.vtk')
pl1 = pv.Plotter()
pl1.add_mesh(mesh)
pl1.show_grid()
pl1.view_xy()
pl1.show()