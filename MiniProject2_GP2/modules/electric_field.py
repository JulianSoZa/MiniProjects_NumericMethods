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