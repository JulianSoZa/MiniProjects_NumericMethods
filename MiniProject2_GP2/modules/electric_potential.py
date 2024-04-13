import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

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
    
    return V, num
