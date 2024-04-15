import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from pyvista import CellType
import meshio

def electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s):
    data = []
    row = []
    col = []

    b = np.zeros(nk)
    
    V_ext = lambda y: np.exp(-0.2*(y-3*np.pi/4)**2) * (np.sin(5*(y**2)/(np.pi)))**2
    V_int = lambda y: np.sin(2*y)

    num = lambda x,y: nth*x+y if (y>=0 and y<nth) else nth*x+nth-1 if y<0 else nth*x

    for k in range(nk):
        j = k%(nth)
        i = int(k/(nth))
        #print(f'{k} = {i}, {j}')
        
        rij = r_dis[i]
        #rij = 2
        clambda = dth/dr
        h = (clambda**2)*dr/rij
        alpha = 1/(rij**2)

        cA = -2*clambda**2 - h - alpha*2
        cB = clambda**2 + h
        cC = clambda**2
        cD = alpha
        cE = alpha
        
        #print(f'A V_{num(i,j)} + B V_{num(i+1,j)} + C V_{num(i-1,j)} + D V_{num(i,j+1)} + E V_{num(i,j-1)}')
        
        #Frontera:
        if(i == 0):
            data.append(1)
            row.append(k)
            col.append(k)
            b[k] = V_int(t_dis[j])
            #print(k)
            continue
        
        if(i == nr):
            data.append(1)
            row.append(k)
            col.append(k)
            b[k] = V_ext(t_dis[j])
            #print(k)
            continue
            
        data.append(cA)
        row.append(k)
        col.append(k)
        
        data.append(cB)
        row.append(k)
        col.append(int(num(i+1,j)))
        
        data.append(cC)
        row.append(k)
        col.append(int(num(i-1,j)))
        
        data.append(cD)
        row.append(k)
        col.append(int(num(i,j+1)))
        
        data.append(cE)
        row.append(k)
        col.append(int(num(i,j-1)))
        
    A = csr_matrix((data, (row, col)), shape=(nk, nk))
    V = spsolve(A,b)

    Vrt = np.zeros((nr+1, nth+1))

    for k in range(nk):
        j = k%(nth)
        i = int(k/(nth))
        Vrt[i,j] = V[k]
        
    # Grafica -----------------

    fig = plt.figure()
    ax = fig.add_subplot()

    cmap_T = 'viridis'

    cb1 = ax.pcolormesh(x_s, y_s, Vrt, shading='auto', cmap=cmap_T)
    fig.colorbar(cb1, ax=ax)

    plt.tight_layout()
    
    # Grafica 2 ------------
    
    r, t = np.meshgrid(r_dis, t_dis)
    
    fig3, ax3 = plt.subplots(subplot_kw={'projection': 'polar'})

    c = ax3.pcolormesh(t, r, Vrt.transpose(), shading='auto', cmap='viridis')

    fig3.colorbar(c, ax=ax3)
    
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
        
    pl = pv.Plotter()
    pl.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
    pl.show_grid()
    pl.view_xy()
    pl.show()
    
    return V, num

def electric_potential_solution_cartesian(elemento, coordenadas, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntos, elementosIndices):
    #print('Puntos: \n', puntos)
    data = []
    row = []
    col = []
    
    nk = nx*ny

    b = np.zeros(nk)
    
    V_ext = lambda y: np.exp(-0.2*(y-3*np.pi/4)**2) * (np.sin(5*(y**2)/(np.pi)))**2
    V_int = lambda y: np.sin(2*y)

    num = lambda x,y: nx*y+x

    alpha = dely/delx
    
    cA = -2*(1+alpha**2)
    cB = alpha**2
    cC = alpha**2
    cD = 1
    cE = 1
    
    for k in range(nk):
        i = k%(nx)
        j = int(k/(nx))
        
        #print(f'A V_{num(i,j)} + B V_{num(i+1,j)} + C V_{num(i-1,j)} + D V_{num(i,j+1)} + E V_{num(i,j-1)}')
        #Frontera:
        #print(k)
        #print(i,j)
        
        radio = np.sqrt(x_dis[i]**2 + y_dis[j]**2)
        
        if k in puntos:
            if(radio >= r_sup):
                data.append(1)
                row.append(k)
                col.append(k)
                if ((y_dis[j]<0)&(x_dis[i]<0)):
                    b[k] = V_ext(np.arctan(y_dis[j]/x_dis[i])+np.pi)
                
                if ((y_dis[j]<0)&(x_dis[i]>0)):
                    b[k] = V_ext(np.arctan(y_dis[j]/x_dis[i])+2*np.pi)
                
                if (((y_dis[j]>0)&(x_dis[i]>0))|((y_dis[j]>0)&(x_dis[i]<0))):
                    b[k] = V_ext(np.arctan2(y_dis[j],x_dis[i]))
                #print(k)
                continue
        
            if(radio <= r_inf):
                data.append(1)
                row.append(k)
                col.append(k)
                if ((y_dis[j]<0)&(x_dis[i]<0)):
                    b[k] = V_int(np.arctan(y_dis[j]/x_dis[i])+np.pi)
                    
                if ((y_dis[j]<0)&(x_dis[i]>0)):
                    b[k] = V_int(np.arctan(y_dis[j]/x_dis[i])+2*np.pi)
                
                if (((y_dis[j]>0)&(x_dis[i]>0))|((y_dis[j]>0)&(x_dis[i]<0))):
                    b[k] = V_int(np.arctan2(y_dis[j],x_dis[i]))
                #print(k)
                continue
            
            radio1 = np.sqrt(x_dis[i+1]**2 + y_dis[j]**2)
            radio2 = np.sqrt(x_dis[i]**2 + y_dis[j-1]**2)
            radio3 = np.sqrt(x_dis[i]**2 + y_dis[j+1]**2)
            radio4 = np.sqrt(x_dis[i-1]**2 + y_dis[j]**2)
            
            if((radio1>=r_sup)) |   ((radio2>=r_sup))   |   ((radio3>=r_sup))   |   ((radio4>=r_sup)):
                #print(k)
                data.append(1)
                row.append(k)
                col.append(k)
                if ((y_dis[j]<0)&(x_dis[i]<0)):
                    b[k] = V_ext(np.arctan(y_dis[j]/x_dis[i])+np.pi)
                
                if ((y_dis[j]<0)&(x_dis[i]>0)):
                    b[k] = V_ext(np.arctan(y_dis[j]/x_dis[i])+2*np.pi)
                
                if (((y_dis[j]>0)&(x_dis[i]>0))|((y_dis[j]>0)&(x_dis[i]<0))):
                    b[k] = V_ext(np.arctan2(y_dis[j],x_dis[i]))
                #print(k)
                continue
            
            if((radio1<=r_inf)) |   ((radio2<=r_inf))   |   ((radio3<=r_inf))   |   ((radio4<=r_inf)):
                #print(k)
                data.append(1)
                row.append(k)
                col.append(k)
                if ((y_dis[j]<0)&(x_dis[i]<0)):
                    b[k] = V_int(np.arctan(y_dis[j]/x_dis[i])+np.pi)
                    
                if ((y_dis[j]<0)&(x_dis[i]>0)):
                    b[k] = V_int(np.arctan(y_dis[j]/x_dis[i])+2*np.pi)
                
                if (((y_dis[j]>0)&(x_dis[i]>0))|((y_dis[j]>0)&(x_dis[i]<0))):
                    b[k] = V_int(np.arctan2(y_dis[j],x_dis[i]))
                    
                continue
            
            data.append(cA)
            row.append(k)
            col.append(k)
            
            data.append(cB)
            row.append(k)
            col.append(int(num(i+1,j)))
            
            data.append(cC)
            row.append(k)
            col.append(int(num(i-1,j)))
            
            data.append(cD)
            row.append(k)
            col.append(int(num(i,j+1)))
            
            data.append(cE)
            row.append(k)
            col.append(int(num(i,j-1)))
        else:
            data.append(1)
            row.append(k)
            col.append(k)
        
    A = csr_matrix((data, (row, col)))
    
    V = spsolve(A,b)
    
    points = np.zeros((nk, 2))
            
    for k in puntos:
        i = k%(nx)
        j = int(k/(nx))
        pt = [x_dis[i], y_dis[j]]
        points[k] = pt
    
    cells = [("quad", elementosIndices)]
    original_mesh = meshio.Mesh(points, cells)
    original_mesh.point_data["Potencial"] = V

    original_mesh_pv = pv.wrap(original_mesh)
        
    pl = pv.Plotter()
    pl.add_mesh(original_mesh_pv, show_edges=False, cmap='viridis', scalars="Potencial")
    pl.view_xy()
    pl.show_grid()
    pl.show()
    
    return V

if __name__ == "__main__":
    from domainDiscretization import cartesian as doCartesian
    import matplotlib.pyplot as plt
    from domainDiscretization import polar as doPolar

    
    nx = 400
    ny = 400
    
    r_inf = 3
    r_sup = 8
    
    elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)
    
    electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, elementosIndices)
    
    r_inf = 3
    r_sup = 8

    nth = 400
    nr = 200
    t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

    electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s)

    plt.show()