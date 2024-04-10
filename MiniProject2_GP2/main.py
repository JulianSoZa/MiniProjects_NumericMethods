import matplotlib.pyplot as plt
from modules import electric_field, electric_potential, plot_electric_solution
from modules.domainDiscretization import cartesian as doCartesian 
from modules.domainDiscretization import polar as doPolar

nx = 10
ny = 20
r_inf = 3
r_sup = 8
elemento, puntos = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

nth = 400
nr = 400
t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

V, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s)

electric_field.electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num, x_s, y_s)

plt.show()