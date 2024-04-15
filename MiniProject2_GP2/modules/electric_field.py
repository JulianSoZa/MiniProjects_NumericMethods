import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import meshio
import pyvista as pv

def electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num, x_s, y_s):
    Er = np.zeros(nk)
    Et = np.zeros(nk)
    En = np.zeros(nk)

    yy = sym.symbols('y')

    V_int = sym.sin(2*yy)
    V_ext = sym.exp(-0.2*(yy-3*sym.pi/4)**2) * (sym.sin(5*(yy**2)/(sym.pi)))**2

    Eint_th = sym.lambdify(yy, sym.diff(V_int,yy), 'numpy')
    Eext_th = sym.lambdify(yy, sym.diff(V_ext,yy), 'numpy')

    for k in range(nk):
        j = k%(nth)
        i = int(k/(nth))
        
        if (i==0):
            Er[k] = -(V[num(i+1,j)] - V[num(i,j)])/(dr)
            Et[k] = -Eint_th(t_dis[j])/r_dis[i]
            En[k] = np.sqrt(Er[k]**2 + Et[k]**2)
            continue
        if (i==nr):
            Er[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dr)
            Et[k] = -Eext_th(t_dis[j])/r_dis[i]
            En[k] = np.sqrt(Er[k]**2 + Et[k]**2)
            continue
        
        Er[k] = -(V[num(i+1,j)] - V[num(i-1,j)])/(2*dr)
        Et[k] = -(V[num(i,j+1)] - V[num(i,j-1)])/(2*dth*r_dis[i])
        En[k] = np.sqrt(Er[k]**2 + Et[k]**2)
    
    print(np.amax(En))
    print(np.amin(En))
    
    # Grafica -----------------
    
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
        
    pl = pv.Plotter()
    pl.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
    pl.show_grid()
    pl.view_xy()
    pl.show()
    
    return En
    
def electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, dx, dy, puntos, r_sup, r_inf, elementosIndices, nk):
    
    Ex = np.zeros(nk)
    Ey = np.zeros(nk)
    En = np.zeros(nk)
    
    num = lambda x,y: nx*y+x

    for k in range(nk):
        i = k%(nx)
        j = int(k/(nx))
        
        if ((y_dis[j]<0)&(x_dis[i]<0)):
            angulo = np.arctan(y_dis[j]/x_dis[i])+np.pi
        
        if ((y_dis[j]<0)&(x_dis[i]>0)):
            angulo = np.arctan(y_dis[j]/x_dis[i])+2*np.pi
        
        if (((y_dis[j]>0)&(x_dis[i]>0))|((y_dis[j]>0)&(x_dis[i]<0))):
            angulo =  np.arctan2(y_dis[j],x_dis[i])
            
        radio = np.sqrt(x_dis[i]**2 + y_dis[j]**2)
        
        if k in puntos:
            if(radio >= r_sup):
                if ((angulo>=0) & (angulo<=np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                continue
            
            if(radio <= r_inf):
                if ((angulo>=0) & (angulo<=np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                continue
            
            radio1 = np.sqrt(x_dis[i+1]**2 + y_dis[j]**2)
            radio2 = np.sqrt(x_dis[i]**2 + y_dis[j-1]**2)
            radio3 = np.sqrt(x_dis[i]**2 + y_dis[j+1]**2)
            radio4 = np.sqrt(x_dis[i-1]**2 + y_dis[j]**2)
            
            if((radio1>=r_sup)&(num(i+1,j) not in puntos))  |   ((radio2>=r_sup)&(num(i,j-1) not in puntos))   |   ((radio3>=r_sup)&(num(i,j+1) not in puntos))   |   ((radio4>=r_sup)&(num(i-1,j) not in puntos)):
                if ((angulo>=0) & (angulo<=np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                continue
            
            if((radio1<=r_inf)&(num(i+1,j) not in puntos)) |   ((radio2<=r_inf)&(num(i,j-1) not in puntos))   |   ((radio3<=r_inf)&(num(i,j+1) not in puntos))   |   ((radio4<=r_inf)&(num(i-1,j) not in puntos)):
                if ((angulo>=0) & (angulo<=np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[k] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[k] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[k] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
                    
                continue
            
            Ex[k] = -(V[num(i+1,j)] - V[num(i-1,j)])/(2*dx)
            Ey[k] = -(V[num(i,j+1)] - V[num(i,j-1)])/(2*dy)
            En[k] = np.sqrt(Ex[k]**2 + Ey[k]**2)
            
    print(np.amax(En))
    print(np.amin(En))
        
    points = np.zeros((nk, 2))
            
    for k in puntos:
        i = k%(nx)
        j = int(k/(nx))
        pt = [x_dis[i], y_dis[j]]
        points[k] = pt
    
    cells = [("quad", elementosIndices)]
    original_mesh = meshio.Mesh(points, cells)
    original_mesh.point_data["Potencial"] = En

    original_mesh_pv = pv.wrap(original_mesh)
        
    pl = pv.Plotter()
    pl.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
    pl.view_xy()
    pl.show_grid()
    pl.show()
    
    return En
    
if __name__ == "__main__":
    from domainDiscretization import cartesian as doCartesian
    import matplotlib.pyplot as plt
    from domainDiscretization import polar as doPolar
    from electric_potential import*
    
    nx = 400
    ny = 400
    
    r_inf = 3
    r_sup = 8
    
    elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)
    
    V = electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, elementosIndices)
    
    electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, delx, dely, puntosIndices, r_sup, r_inf, elementosIndices)