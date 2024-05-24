import matplotlib.pyplot as plt
import pandas as pd
from modules import electric_field, electric_potential, plot_electric_solution
from modules.domainDiscretization import cartesian as doCartesian 
from modules.domainDiscretization import polar as doPolar
import analysis_results

# ----------------- Coordenadas Polares -------------------------------

r_inf = 3
r_sup = 8

nth = 400
nr = 400

t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

V, V_space, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk)
plot_electric_solution.ploter_finite_solutions(nk, nth, r_dis, t_dis, nr, num, V, 'Potencial')

df = pd.DataFrame(V_space)
df.columns = ['r', 'theta ', 'V']
df.to_csv('Potencial_Electrico_Polares.txt', sep='\t')

En, E_space = electric_field.electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num)
plot_electric_solution.ploter_finite_solutions(nk, nth, r_dis, t_dis, nr, num, En, 'Campo')

df = pd.DataFrame(E_space)
df.columns = ['r', 'theta', 'Er', 'Et']
df.to_csv('Campo_Electrico_Polares.txt', sep='\t')

# -------------------------------------------------------------------------

# ----------------- Coordenadas Cartesianas -------------------------------

nx = 400
ny = 400

elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

V_c, V_space = electric_potential.electric_potential_solution_cartesian(nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, nk)
plot_electric_solution.ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, V_c,  'Potencial')

df = pd.DataFrame(V_space)
df.columns = ['x', 'y', 'V']
df.to_csv('Potencial_Electrico_Cartesianas.txt', sep='\t')

En_c, E_space = electric_field.electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V_c, delx, dely, puntosIndices, r_sup, r_inf, nk)
plot_electric_solution.ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, En_c, 'Campo')

df = pd.DataFrame(E_space)
df.columns = ['x', 'y', 'Ex', 'Ey']
df.to_csv('Campo_Electrico_Cartesianas.txt', sep='\t')

# -------------------------------------------------------------------------

# --------------- Grafias del analisis de resultados -----------------------

analysis_results.results_analysis_graphs()