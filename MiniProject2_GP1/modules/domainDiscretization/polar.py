import numpy as np
import matplotlib.pyplot as plt

def polar_discretization(nth, nr):
    r_inf = 3
    r_sup = 8
    
    dth = 2*np.pi/nth
    dr = (r_sup-r_inf)/nr
    
    nk = nth*(nr+1) #Numero de puntos en numeracion global

    t_dis = np.linspace(0, (2*np.pi), nth+1)           #Discretizacion de los Angulos
    r_dis = np.linspace(r_inf, r_sup, nr+1)           #Discretizacion de los radios

    th, r = np.meshgrid(t_dis, r_dis)            #Malla

    x_s = r * np.cos(th)                         #Convercion a carteciano (grafica)
    y_s = r* np.sin(th)                         #Convercion a carteciano (grafica)
    
    print('Se ha discretizado el dominio polarmente\n')
    
    return t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk

if __name__ == "__main__":
    
    nth = 15
    nr = 12
    
    t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = polar_discretization(nth, nr)
    
    fig, ax = plt.subplots(figsize=(6, 6))

    for i in range(nr+1):
        ax.plot(x_s[i], y_s[i], color='firebrick')

    for j in range(nth):
        ax.plot(x_s[:, j], y_s[:, j], color='firebrick')
    ax.set_aspect('equal')  
    ax.grid(color='lightgray', linestyle='--', linewidth=1)
    plt.show()