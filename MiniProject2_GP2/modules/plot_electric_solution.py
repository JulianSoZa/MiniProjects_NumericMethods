import numpy as np
import meshio
import pyvista as pv
import electric_field, electric_potential
from domainDiscretization import cartesian as doCartesian 
from domainDiscretization import polar as doPolar

r_inf = 3
r_sup = 8

nth = 400
nr = 400

t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

V, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s)

En = electric_field.electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num, x_s, y_s)

#Grafica Polares Potencial Electrico ----------------------
ps = []
rs = []

for k in range(nk):
    j = int(k%(nth))
    i = int(k/(nth))

    x_dis = round(r_dis[i]*np.cos(t_dis[j]),4)
    y_dis = round(r_dis[i]*np.sin(t_dis[j]),4)
    ps.append([x_dis, y_dis])

for k in range(nth*nr):
    j = int(k%(nth))
    i = int(k/(nth))
    
    n1 = int(num(i,j))
    n2 = int(num(i+1,j))
    n3 = int(num(i+1,j+1))
    n4 = int(num(i,j+1))
    rs.append([n1, n2, n3, n4])

cells = [("quad", rs)]
original_mesh = meshio.Mesh(ps, cells)
original_mesh.point_data["Potencial"] = V

original_mesh_pv = pv.wrap(original_mesh)
    
pl1 = pv.Plotter()
pl1.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
pl1.show_grid()
pl1.view_xy()

#Grafica polares Campo Electrico ---------------------------

ps = []
rs = []

for k in range(nk):
    j = int(k%(nth))
    i = int(k/(nth))

    x_dis = round(r_dis[i]*np.cos(t_dis[j]),4)
    y_dis = round(r_dis[i]*np.sin(t_dis[j]),4)
    ps.append([x_dis, y_dis])

for k in range(nth*nr):
    j = int(k%(nth))
    i = int(k/(nth))
    
    n1 = int(num(i,j))
    n2 = int(num(i+1,j))
    n3 = int(num(i+1,j+1))
    n4 = int(num(i,j+1))
    rs.append([n1, n2, n3, n4])

cells = [("quad", rs)]
original_mesh = meshio.Mesh(ps, cells)
original_mesh.point_data["Potencial"] = En

original_mesh_pv = pv.wrap(original_mesh)
    
pl3 = pv.Plotter()
pl3.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
pl3.show_grid()
pl3.view_xy()

# --------------------------------------------------------------------------------

nx = 400
ny = 400

elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

V = electric_potential.electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, elementosIndices, nk)

En = electric_field.electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, delx, dely, puntosIndices, r_sup, r_inf, elementosIndices, nk)

#Grafica Cartesianas Potencial Electrico ----------------------

points = np.zeros((nk, 2))
            
for k in puntosIndices:
    i = k%(nx)
    j = int(k/(nx))
    pt = [x_dis[i], y_dis[j]]
    points[k] = pt

cells = [("quad", elementosIndices)]
original_mesh = meshio.Mesh(points, cells)
original_mesh.point_data["Potencial"] = V

original_mesh_pv = pv.wrap(original_mesh)
    
pl2 = pv.Plotter()
pl2.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
pl2.view_xy()
pl2.show_grid()

#Grafica cartesianas Campo Electrico ----------------------

points = np.zeros((nk, 2))
            
for k in puntosIndices:
    i = k%(nx)
    j = int(k/(nx))
    pt = [x_dis[i], y_dis[j]]
    points[k] = pt

cells = [("quad", elementosIndices)]
original_mesh = meshio.Mesh(points, cells)
original_mesh.point_data["Potencial"] = En

original_mesh_pv = pv.wrap(original_mesh)
    
pl4 = pv.Plotter()
pl4.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
pl4.view_xy()
pl4.show_grid()

pl1.show()
pl2.show()
pl3.show()
pl4.show()