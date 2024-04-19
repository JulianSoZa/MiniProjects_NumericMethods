if __name__ == "__main__":
    import electric_field, electric_potential
    from domainDiscretization import cartesian as doCartesian 
    from domainDiscretization import polar as doPolar
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import griddata

    r_inf = 3
    r_sup = 8
    
    disFina = 400
    iteraciones = np.array([5, 10, 20, 50, 100, disFina])
    
    fig, ax = plt.subplots()
    for h in iteraciones:
        nth = h
        nr = h
        
        t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

        V, V_space, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk)
        
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
        if h == iteraciones[-1]:
            ax.plot(t_r7, V_r7, label=f'nth y nr = {h}')
            ax.legend()
        else:
            ax.plot(t_r7, V_r7, '--', label=f'nth y nr = {h}')
            ax.legend()
        print(f'Iteración para nr y nth = {h} terminada')
    
    #------------ interpolacion ---------------------------
    
    fig2, ax2 = plt.subplots()
    
    for h in iteraciones:
        nth = h
        nr = h
    
        t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

        V, V_space, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk)

        ri = np.array([7])
        ti = np.linspace(0, 2*np.pi, nth+1)
        xi, yi = np.meshgrid(ri, ti)
        
        V_ij = np.zeros((nr+1, nth+1))
        
        for k in range(nk):
            j = k%(nth)
            i = int(k/(nth))
            V_ij[i,j] = V[k]
        
        zi = griddata((r.flatten(), th.flatten()), V_ij.flatten(), (xi, yi), method='cubic')
        
        if h == iteraciones[-1]:
            ax2.plot(t_dis, zi, label=f'Interpolacion nth y nr: {h}')
            ax2.legend()
        else:
            ax2.plot(t_dis, zi, '--', label=f'Interpolacion nth y nr: {h}')
            ax2.legend()
        print(f'Iteración con interpolación para nr y nth = {h} terminada')

    #-------- cartesianas ---------------------------------------
    iteraciones = np.array([5, 10, 15, 20, 35, 50, 75, 100, 200, disFina])
    
    fig3, ax3 = plt.subplots()
    
    MSE_ns = np.array([])
    nks = np.array([])
    
    for h in iteraciones:
        nx = h
        ny = h

        elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

        V_c, V_space_c = electric_potential.electric_potential_solution_cartesian(nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, nk)
        
        # -------------------- Interpolacion ------------------
        
        t_dis = np.linspace(0, 2*np.pi, iteraciones[-1]+1)
        x_r7 = 7*np.cos(-t_dis+np.pi/2)
        y_r7 = 7*np.sin(-t_dis+np.pi/2)
        
        x_s, y_s = np.meshgrid(x_dis, y_dis)
        
        V_ij = np.zeros((nx, ny))
        
        for k in range(nk):
            i = k%(nx)
            j = int(k/(nx))
            V_ij[i,j] = V_c[k]
        
        ziC = griddata((x_s.flatten(), y_s.flatten()), V_ij.flatten(), (x_r7, y_r7), method='cubic')
        
        #print(f'x_r7:\n{x_r7}\ny_r7:\n{y_r7}\nt_dis:\n{t_dis}\n')
        if h == iteraciones[-1]:
            ax3.plot(t_dis, ziC, label=f'Interpolacion nx = {nx} y ny = {ny}')
            ax3.legend()
        else:
            ax3.plot(t_dis, ziC, '--', label=f'Interpolacion nx = {nx} y ny = {ny}')
            ax3.legend()
        print(f'Iteración con interpolación para nx y ny = {h} terminada')
        
        # ---------- Error cuadrático medio --------------------
        suma_err = 0
        m = iteraciones[-1]+1
        
        for k in range(m):
            suma_err += (zi[k]-ziC[k])**2
        
        MSE_ns = np.append(MSE_ns, (1/m)*suma_err)
        nks = np.append(nks, h**2)
        
    print('MSE_ns:\n', MSE_ns)
    
    fig4, ax4 = plt.subplots()
    ax4.semilogx(nks, MSE_ns)

    plt.show()