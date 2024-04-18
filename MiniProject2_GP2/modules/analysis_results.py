if __name__ == "__main__":
    import electric_field, electric_potential
    from domainDiscretization import cartesian as doCartesian 
    from domainDiscretization import polar as doPolar
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import griddata

    r_inf = 3
    r_sup = 8
    
    iteraciones = np.array([5, 10, 20, 50, 100, 500])
    
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
                
        #print(f't_r7:\n{t_r7}\nV_r7:\n{V_r7}\nr_r7:\n{r_r7}\nr_dis:\n{r_dis}')
        #print('dr\n',dr)
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
           
    plt.show()
    
    #--------
    
    nx = 6
    ny = 6

    elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

    V = electric_potential.electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, elementosIndices, nk)
    
    V_r7 = np.array([])
    t_r7 = np.array([])
    r_r7 = np.array([])
    V_th_r7 = np.empty((0,2), int)
    
    Dis_r7 = np.min(abs(7-r_dis))
    
    xx, yy = np.meshgrid(x_dis, y_dis)
    radios = np.sqrt(xx**2 + yy**2)
    print('radios:\n', radios)
    Dis_r7 = np.min(abs(7-radios))
    print('Dis_r7:\n',Dis_r7)
    
    for k in range(nk):
        i = k%(nx)
        j = int(k/(nx))
        
        radio = np.sqrt(x_dis[i]**2 + y_dis[j]**2)
        if ((y_dis[j]<0)&(x_dis[i]<0)):
            angulo = np.arctan(y_dis[j]/x_dis[i])+np.pi
        
        if ((y_dis[j]<0)&(x_dis[i]>0)):
            angulo = np.arctan(y_dis[j]/x_dis[i])+2*np.pi
        
        if (((y_dis[j]>0)&(x_dis[i]>0))|((y_dis[j]>0)&(x_dis[i]<0))):
            angulo =  np.arctan2(y_dis[j],x_dis[i])
        
        if(abs(7-radio)<=Dis_r7):
            V_r7 = np.append(V_r7, V[k])
            t_r7 = np.append(t_r7, angulo)
            r_r7 = np.append(r_r7, radio)
            V_th_r7 = np.append(V_th_r7, np.array([[angulo, V[k]]]), axis = 0)
    
    print('Vrtsdw: \n',V_th_r7)
            
    V_th_r7 = V_th_r7[V_th_r7[:, 0].argsort()]
    t_r7 = V_th_r7[:, 0]
    V_r7 = V_th_r7[:, 1]
            
    print(f't_r7:\n{t_r7}\nV_r7:\n{V_r7}\nr_r7:\n{r_r7}')
    print('dr\n',dr)
    plt.plot(t_r7, V_r7)
    plt.show()