from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon 

nx = 80
ny = 80

r_inf = 3
r_sup = 8

x_new = []
y_new = []
elemento = []
P_totales = []

puntos = []
puntos2 = []

x_dis = np.linspace(-r_sup, r_sup, nx)           #Discretizacion del eje X
y_dis = np.linspace(-r_sup, r_sup, ny)           #Discretizacion del eje Y 
delx = np.abs(x_dis[0] - x_dis[1])               #Paso de x
dely = np.abs(y_dis[0] - y_dis[1])               #paso de y 

print('Delta X:',delx)
print(f'Delta Y: {dely} \n')



for i in range(ny):
    for j in range(nx):
        P_totales.append([x_dis[j], y_dis[i]])


for i in range(len(P_totales)-nx-1):
    if ( (i+1)%(nx) != 0 and r_inf<=np.sqrt((((P_totales[i][0]+P_totales[i+1][0])/2)**2) + (((P_totales[i][1] + P_totales[i+nx][1])/2)**2))<=r_sup):
        
        A = [P_totales[i][0], P_totales[i][1]]
        B = [P_totales[i+1][0], P_totales[i+1][1]]
        C = [P_totales[i+1+nx][0], P_totales[i+1+nx][1]]
        D = [P_totales[i+nx][0], P_totales[i+nx][1]]
        elemento.append([A,B,C,D])
        puntos.extend([A,B,C,D])
puntos = list(OrderedDict.fromkeys(map(tuple, puntos)))
puntos.sort()



 

print('------------------------------------------------------------')
print('Seccion de pruebas\n')

print('Puntos', puntos)
print('segundo elemento',elemento[1])

print('\n------------------------------------------------------------')



#-------------------------------------------------------------------------------Grafica-------------------------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(6, 6))
circle_inf = plt.Circle((0, 0), r_inf, color='blue', fill=False)
circle_sup = plt.Circle((0, 0), r_sup, color='blue', fill=False)

for element in elemento:
    square = Polygon(element, edgecolor='black', facecolor='brown', alpha=0.5)
    plt.gca().add_patch(square)


ax2.set_xlim(-r_sup - 1, r_sup + 1)
ax2.set_ylim(-r_sup - 1, r_sup + 1)
ax2.add_artist(circle_inf)
ax2.add_artist(circle_sup)
ax2.grid(color='lightgray', linestyle='none', linewidth=1)  
plt.show()