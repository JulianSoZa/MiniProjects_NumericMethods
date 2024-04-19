import numpy as np
import meshio
import pyvista as pv

def ploter_finite_solutions(nk, nth, r_dis, t_dis, nr, num, V, name):
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
    original_mesh_pv.save(f'Malla_{name}_Polares.vtk')

    labels = dict(xlabel='X', ylabel='Y')


    pl1 = pv.Plotter()
    pl1.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
    pl1.show_grid()
    pl1.add_title(f"{name} electrico en polares")
    pl1.show_grid(**labels)
    pl1.view_xy()
    pl1.show()
    return

def ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, V):
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
    pl2.show()
    return

if __name__ == '__main__':
    import electric_field, electric_potential
    from domainDiscretization import cartesian as doCartesian 
    from domainDiscretization import polar as doPolar

    r_inf = 3
    r_sup = 8

    nth = 100
    nr = 100

    nx = 100
    ny = 100

    t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)
    V, V_sapce, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk)
    En, E_space = electric_field.electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num)

    ploter_finite_solutions(nk, nth, r_dis, t_dis, nr, num, V, 'Potencial')
    ploter_finite_solutions(nk, nth, r_dis, t_dis, nr, num, En, 'Campo')

    elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)
    V_c, V_space_c = electric_potential.electric_potential_solution_cartesian(nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, nk)
    En_c, E_space_c = electric_field.electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V_c, delx, dely, puntosIndices, r_sup, r_inf, nk)

    ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, V_c)
    ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, En_c)
    