import numpy as np
import matplotlib.pyplot as plt

def polar_discretization(nth, nr):
    r_inf = 3
    r_sup = 8

    t_dis = np.linspace(0, (2*np.pi), nth)           #Discretizacion de los Angulos
    r_dis = np.linspace(r_inf, r_sup, nr)           #Discretizacion de los radios

    th, r = np.meshgrid(t_dis, r_dis)            #Malla

    x_s = r * np.cos(th)                         #Convercion a carteciano (grafica)
    y_s = r* np.sin(th)                         #Convercion a carteciano (grafica)
    
    return t_dis, r_dis, th, r, x_s, y_s

if __name__ == "__main__":
    
    nth = 15
    nr = 12
    
    t_dis, r_dis, th, r, x_s, y_s = polar_discretization(nth, nr)
    
    #-------------------------------------------------------------------------------Grafica-------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 6))

    for i in range(nr):
        ax.plot(x_s[i], y_s[i], color='firebrick')

    for j in range(nth):
        ax.plot(x_s[:, j], y_s[:, j], color='firebrick')
    ax.set_aspect('equal')  
    ax.grid(color='lightgray', linestyle='--', linewidth=1)
    plt.show()