import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

r_inf = 3
r_sup = 8

nth = 400 # Discretizacion de elementos en la parte angular
nr = 200 # Discretizacion de elementos en la parte radial

dth = 2*np.pi/nth
dr = (r_sup-r_inf)/nr

data = []
row = []
col = []

nk = nth*(nr+1) #Numero de puntos en numeracion global
b = np.zeros(nk)

t_dis = np.linspace(0, (2*np.pi), nth+1)           #Discretizacion de los Angulos
r_dis = np.linspace(r_inf, r_sup, nr+1)           #Discretizacion de los radios

#print(t_dis)

V_ext = lambda y: np.exp(-0.2*(y-3*np.pi/4)**2) * (np.sin(5*(y**2)/(2*np.pi)))**2
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

fig = plt.figure()
ax = fig.add_subplot()

th, r = np.meshgrid(t_dis, r_dis)            #Malla

x_s = r * np.cos(th)                         #Convercion a carteciano (grafica)
y_s = r* np.sin(th)                         #Convercion a carteciano (grafica)

cmap_T = 'viridis'

cb1 = ax.pcolormesh(x_s, y_s, Vrt, shading='auto', cmap=cmap_T)
fig.colorbar(cb1, ax=ax)

plt.tight_layout()

Er = np.zeros((nr+1, nth+1))
Et = np.zeros((nr+1, nth+1))
En = np.zeros((nr+1, nth+1))

yy = sym.symbols('y')

V_int = sym.sin(2*yy)
V_ext = sym.exp(-0.2*(yy-3*sym.pi/4)**2) * (sym.sin(5*(yy**2)/(2*sym.pi)))**2

Eint_th = sym.lambdify(yy, sym.diff(V_int,yy), 'numpy')
Eext_th = sym.lambdify(yy, sym.diff(V_ext,yy), 'numpy')

for k in range(nk):
    j = k%(nth)
    i = int(k/(nth))
    
    if (i==0):
        Er[i, j] = 0
        Et[i, j] = -Eint_th(t_dis[j])/r_dis[i]
        En[i, j] = np.sqrt(Er[i, j]**2 + Et[i, j]**2)
        print(En[i, j])
        continue
    if (i==nr):
        Er[i, j] = 0
        Et[i, j] = -Eext_th(t_dis[j])/r_dis[i]
        En[i, j] = np.sqrt(Er[i, j]**2 + Et[i, j]**2)
        continue
    
    Er[i, j] = -(V[num(i+1,j)] - V[num(i,j)])/(dr)
    Et[i, j] = -(V[num(i,j+1)] - V[num(i,j)])/(dth*r_dis[i])
    En[i, j] = np.sqrt(Er[i, j]**2 + Et[i, j]**2)
    
fig2 = plt.figure()
ax2 = fig2.add_subplot()

cmap_T = 'viridis'

cb2 = ax2.pcolormesh(x_s, y_s, Et, shading='auto', cmap=cmap_T)
fig.colorbar(cb2, ax=ax2)

plt.tight_layout()
plt.show()
