import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

def electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num, x_s, y_s):
    Er = np.zeros((nr+1, nth+1))
    Et = np.zeros((nr+1, nth+1))
    En = np.zeros((nr+1, nth+1))

    yy = sym.symbols('y')

    V_int = sym.sin(2*yy)
    V_ext = sym.exp(-0.2*(yy-3*sym.pi/4)**2) * (sym.sin(5*(yy**2)/(sym.pi)))**2

    Eint_th = sym.lambdify(yy, sym.diff(V_int,yy), 'numpy')
    Eext_th = sym.lambdify(yy, sym.diff(V_ext,yy), 'numpy')

    for k in range(nk):
        j = k%(nth)
        i = int(k/(nth))
        
        if (i==0):
            Er[i, j] = -(V[num(i+1,j)] - V[num(i,j)])/(dr)
            Et[i, j] = -Eint_th(t_dis[j])/r_dis[i]
            En[i, j] = np.sqrt(Er[i, j]**2 + Et[i, j]**2)
            continue
        if (i==nr):
            Er[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dr)
            Et[i, j] = -Eext_th(t_dis[j])/r_dis[i]
            En[i, j] = np.sqrt(Er[i, j]**2 + Et[i, j]**2)
            continue
        
        Er[i, j] = -(V[num(i+1,j)] - V[num(i,j)])/(dr)
        Et[i, j] = -(V[num(i,j+1)] - V[num(i,j)])/(dth*r_dis[i])
        En[i, j] = np.sqrt(Er[i, j]**2 + Et[i, j]**2)
        
    fig2 = plt.figure()
    ax2 = fig2.add_subplot()

    cmap_T = 'viridis'

    cb2 = ax2.pcolormesh(x_s, y_s, En, shading='auto', cmap=cmap_T)
    fig2.colorbar(cb2, ax=ax2)

    plt.tight_layout()
    
def electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, dx, dy, x_s, y_s, puntos, r_sup, r_inf):
    Ex = np.zeros((nx, ny))
    Ey = np.zeros((nx, ny))
    En = np.zeros((nx, ny))
    
    nk = nx*ny
    
    num = lambda x,y: nx*y+x

    for k in range(nk):
        j = k%(nx)
        i = int(k/(nx))
        
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
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                continue
            
            if(radio <= r_inf):
                if ((angulo>=0) & (angulo<=np.pi/2)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                continue
            
            radio1 = np.sqrt(x_dis[i+1]**2 + y_dis[j]**2)
            radio2 = np.sqrt(x_dis[i]**2 + y_dis[j-1]**2)
            radio3 = np.sqrt(x_dis[i]**2 + y_dis[j+1]**2)
            radio4 = np.sqrt(x_dis[i-1]**2 + y_dis[j]**2)
            
            if((radio1>=r_sup)&(radio2>=r_sup)) |   ((radio1>=r_sup)&(radio3>=r_sup))   |   ((radio4>=r_sup)&(radio3>=r_sup))   |   ((radio4>=r_sup)&(radio2>=r_sup)):
                if ((angulo>=0) & (angulo<=np.pi/2)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                continue
            
            if((radio1<=r_inf)&(radio2<=r_inf)) |   ((radio1<=r_inf)&(radio3<=r_inf))   |   ((radio4<=r_inf)&(radio3<=r_inf))   |   ((radio4<=r_inf)&(radio2<=r_inf)):
                if ((angulo>=0) & (angulo<=np.pi/2)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi/2) & (angulo<=np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j+1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>np.pi) & (angulo<=3*np.pi/2)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i-1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                if ((angulo>3*np.pi/2) & (angulo<2*np.pi)):
                    Ex[i, j] = -(V[num(i,j)] - V[num(i+1,j)])/(dx)
                    Ey[i, j] = -(V[num(i,j)] - V[num(i,j-1)])/(dy)
                    En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
                    
                continue
            
            Ex[i, j] = -(V[num(i+1,j)] - V[num(i-1,j)])/(2*dx)
            Ey[i, j] = -(V[num(i,j+1)] - V[num(i,j-1)])/(2*dy)
            En[i, j] = np.sqrt(Ex[i, j]**2 + Ey[i, j]**2)
        
    fig2 = plt.figure()
    ax2 = fig2.add_subplot()

    cmap_T = 'viridis'

    cb2 = ax2.pcolormesh(x_s, y_s, En, shading='auto', cmap=cmap_T)
    fig2.colorbar(cb2, ax=ax2)

    plt.tight_layout()
    
if __name__ == "__main__":
    from domainDiscretization import cartesian as doCartesian
    import matplotlib.pyplot as plt
    from domainDiscretization import polar as doPolar
    from electric_potential import*

    
    nx = 400
    ny = 400
    
    r_inf = 3
    r_sup = 8
    
    elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)
    
    V, x_s, y_s = electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices)
    
    electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, delx, dely, x_s, y_s, puntosIndices, r_sup, r_inf)

    plt.show()