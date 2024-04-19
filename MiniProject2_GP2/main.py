import matplotlib.pyplot as plt
import pandas as pd
from modules import electric_field, electric_potential
from modules.domainDiscretization import cartesian as doCartesian 
from modules.domainDiscretization import polar as doPolar

r_inf = 3
r_sup = 8

nth = 100
nr = 100

t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

V, V_space, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk)

df = pd.DataFrame(V_space)
df.columns = ['r', 'theta ', 'V']
df.to_csv('Potencial_Electrico_Polares.txt', sep='\t')

En, E_space = electric_field.electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num)

df = pd.DataFrame(E_space)
df.columns = ['r', 'theta', 'Er', 'Et']
df.to_csv('Campo_Electrico_Polares.txt', sep='\t')

nx = 100
ny = 100

elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

V, V_space = electric_potential.electric_potential_solution_cartesian(nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, nk)

df = pd.DataFrame(V_space)
df.columns = ['x', 'y', 'V']
df.to_csv('Potencial_Electrico_Cartesianas.txt', sep='\t')

En, E_space = electric_field.electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, delx, dely, puntosIndices, r_sup, r_inf, nk)

df = pd.DataFrame(E_space)
df.columns = ['x', 'y', 'Ex', 'Ey']
df.to_csv('Campo_Electrico_Cartesianas.txt', sep='\t')