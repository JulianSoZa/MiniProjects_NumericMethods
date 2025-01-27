import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from collections import OrderedDict

def cartesian_discretization(nx, ny, r_inf, r_sup):
    elemento = []
    P_totales = []
    
    elementosIndices = []
    puntosIndices = []

    coordenadas = []
    
    nk = nx*ny

    x_dis = np.linspace(-r_sup, r_sup, nx)           #Discretizacion del eje X
    y_dis = np.linspace(-r_sup, r_sup, ny)           #Discretizacion del eje Y 
    delx = np.abs(x_dis[0] - x_dis[1])               #Paso de x
    dely = np.abs(y_dis[0] - y_dis[1])               #paso de y 

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
            coordenadas.extend([A,B,C,D])
            
            n1 = i
            n2 = i+1
            n3 = i+1+nx
            n4 = i+nx
            elementosIndices.append([n1, n2, n3, n4])
            puntosIndices.append(n1)
            puntosIndices.append(n2)
            puntosIndices.append(n3)
            puntosIndices.append(n4)
            
    puntosIndices = np.unique(puntosIndices)
            
    coordenadas = list(OrderedDict.fromkeys(map(tuple, coordenadas)))
    coordenadas.sort()
    
    print('Se ha discretizado el domino cartesianamente\n')
    
    return elemento, coordenadas, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk
    
def cartesian_discretization_Bonus(nx, ny, r_inf, r_sup, p):
    elemento = []
    P_totales = []
    elementosIndices = []
    puntosIndices = []
    coordenadas = []

    puntos = []
    nk = nx*ny

    x_dis = np.linspace(-r_sup, r_sup, nx)           #Discretizacion del eje X
    y_dis = np.linspace(-r_sup, r_sup, ny)           #Discretizacion del eje Y 
    delx = np.abs(x_dis[0] - x_dis[1])               #Paso de x
    dely = np.abs(y_dis[0] - y_dis[1])               #paso de y 

    for i in range(ny):
        for j in range(nx):
            P_totales.append([x_dis[j], y_dis[i]])


    for i in range(len(P_totales)-nx-1):

        A = [P_totales[i][0], P_totales[i][1]]
        B = [P_totales[i+1][0], P_totales[i+1][1]]
        C = [P_totales[i+1+nx][0], P_totales[i+1+nx][1]]
        D = [P_totales[i+nx][0], P_totales[i+nx][1]]
        filtro = [A, B, C, D]

        # Contador para contar cuántos puntos cumplen la condición
        puntos_cumplen = sum(1 for punto in filtro if punto[0]**2 + punto[1]**2 >= r_inf ** 2 and punto[0] ** 2 + punto[1] ** 2 <= r_sup ** 2)

        # Si al menos dos puntos cumplen la condición, se añade el elemento
        if (i+1)%(nx) != 0 and puntos_cumplen >= p:
            elemento.append([A, B, C, D])
            puntos.extend([A, B, C, D])

            n1 = i
            n2 = i+1
            n3 = i+1+nx
            n4 = i+nx
            elementosIndices.append([n1, n2, n3, n4])
            puntosIndices.append(n1)
            puntosIndices.append(n2)
            puntosIndices.append(n3)
            puntosIndices.append(n4)

    puntosIndices = np.unique(puntosIndices)
            
    coordenadas = list(OrderedDict.fromkeys(map(tuple, coordenadas)))
    coordenadas.sort()

    print('Se ha discretizado el domino cartesianamente\n')

    return elemento, coordenadas, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk

if __name__ == "__main__":
    nx = 10
    ny = 10
    p = 4
    r_inf = 3
    r_sup = 8

    elemento, coordenada, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = cartesian_discretization(nx, ny, r_inf, r_sup)
    #elemento, coordenada, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = cartesian_discretization_Bonus(nx, ny, r_inf, r_sup, p)
    
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