if __name__ == "__main__":
    import electric_field, electric_potential
    from domainDiscretization import cartesian as doCartesian 
    from domainDiscretization import polar as doPolar
    import matplotlib.pyplot as plt
    import numpy as np

    r_inf = 3
    r_sup = 8

    nth = 400
    nr = 400
    t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

    V, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s)
    
    V_r7 = np.zeros(nth+1)
    ks_r7 = []
    
    for k in range(nk):
        j = k%(nth)
        i = int(k/(nth))
        
        if ((r_dis[i]>=7)&(7<=r_dis[i]+dr)):
            V_r7[j] = V[k]
    
    plt.plot(t_dis, V_r7)
    plt.show()
    
"""nx = 400
ny = 400

elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

V = electric_potential.electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, elementosIndices, nk)
    
"""