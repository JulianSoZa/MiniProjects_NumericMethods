import numpy as np
import sympy as sym

def electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num):
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
    
    E_space = np.zeros((nk,4))
    
    for k in range(nk):
        j = k%(nth)
        i = int(k/(nth))
        E_space[k] = np.array([r_dis[i], t_dis[j], Er[k], Et[k]])
        
    return En, E_space
    
def electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, dx, dy, puntos, r_sup, r_inf, nk):
    
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
    
    E_space = np.zeros((nk,4))
    
    for k in range(nk):
        i = k%(nx)
        j = int(k/(nx))
        E_space[k] = np.array([x_dis[i], y_dis[j], Ex[k], Ey[k]])
    
    return En, E_space
    
if __name__ == "__main__":
    from domainDiscretization import cartesian as doCartesian
    from domainDiscretization import polar as doPolar
    import electric_potential, plot_electric_solution
    
    r_inf = 3
    r_sup = 8

    nth = 100
    nr = 100

    nx = 100
    ny = 100

    t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)
    V, V_sapce, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk)
    En, E_space = electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num)

    plot_electric_solution.ploter_finite_solutions(nk, nth, r_dis, t_dis, nr, num, En, 'Campo')

    elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)
    V_c, V_space_c = electric_potential.electric_potential_solution_cartesian(nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, nk)
    En_c, E_space_c = electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V_c, delx, dely, puntosIndices, r_sup, r_inf, nk)

    plot_electric_solution.ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, En_c)