if __name__ == "__main__":
    import electric_field, electric_potential
    from domainDiscretization import cartesian as doCartesian 
    from domainDiscretization import polar as doPolar
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import griddata

    r_inf = 3
    r_sup = 8
    
    """iteraciones = np.array([5, 10, 20, 50, 100, 500])
    
    fig, ax = plt.subplots()
    for j in iteraciones+1:
        nth = j
        nr = j
        
        t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

        V, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s)
        
        V_r7 = np.array([])
        t_r7 = np.array([])
        r_r7 = np.array([])
        
        Dis_r7 = np.min(abs(7-r_dis))
        
        for k in range(nk):
            j = k%(nth)
            i = int(k/(nth))
            
            if(abs(7-r_dis[i])<=Dis_r7):
                V_r7 = np.append(V_r7, V[k])
                t_r7 = np.append(t_r7, t_dis[j])
                r_r7 = np.append(r_r7, r_dis[i])
                
        V_r7 = np.append(V_r7, V_r7[0])
        t_r7 = np.append(t_r7, 2*np.pi)
        if j == iteraciones[-1]:
            ax.plot(t_r7, V_r7, label=f'n y nr = {j}')
            ax.legend()
        else:
            ax.plot(t_r7, V_r7, '--', label=f'nth y nr = {j}')
            ax.legend()
        print(f'Iteración para nr y nth = {j} terminada')
    
    #------------ interpolacion ---------------------------
    
    iteraciones = np.array([5, 10, 20, 50, 100, 500])
    
    fig2, ax2 = plt.subplots()
    
    for j in iteraciones+1:
        nth = j
        nr = j
    
        t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

        V, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s)

        ri = np.array([7])
        ti = np.linspace(0, 2*np.pi, nth+1)
        xi, yi = np.meshgrid(ri, ti)
        
        V_ij = np.zeros((nr+1, nth+1))
        
        for k in range(nk):
            j = k%(nth)
            i = int(k/(nth))
            V_ij[i,j] = V[k]
        
        zi = griddata((r.flatten(), th.flatten()), V_ij.flatten(), (xi, yi), method='cubic')
        
        if j == iteraciones[-1]:
            ax2.plot(t_dis, zi, label=f'Interpolacion nth y nr: {j}')
            ax2.legend()
        else:
            ax2.plot(t_dis, zi, '--', label=f'Interpolacion nth y nr: {j}')
            ax2.legend()
        print(f'Iteración con interpolación para nr y nth = {j} terminada')
           
    plt.show()"""
    
    #--------
    
    nx = 6
    ny = 6

    elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

    V = electric_potential.electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, elementosIndices, nk)
    
    # -------------------- Interpolacion ------------------
    
    t_dis = np.linspace(0, 2*np.pi, 100)
    x_r7 = 7*np.cos(-t_dis+np.pi/2)
    y_r7 = 7*np.sin(-t_dis+np.pi/2)
    
    x_s, y_s = np.meshgrid(x_dis, y_dis)
    
    V_ij = np.zeros((nx, ny))
    
    for k in range(nk):
        i = k%(nx)
        j = int(k/(nx))
        V_ij[i,j] = V[k]
    
    zi = griddata((x_s.flatten(), y_s.flatten()), V_ij.flatten(), (x_r7, y_r7), method='cubic')
    
    #print(f'x_r7:\n{x_r7}\ny_r7:\n{y_r7}\nt_dis:\n{t_dis}\n')
    
    plt.plot(t_dis, zi, label=f'Interpolacion nx = {nx} y ny = {ny}')
    plt.show()