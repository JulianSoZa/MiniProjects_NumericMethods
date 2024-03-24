from modules import electric_field, electric_potential, plot_electric_solution
from modules.domainDiscretization import cartesian_discretization, polar_discretization

nx = 10
ny = 20
r_inf = 3
r_sup = 8
elemento, puntos = cartesian_discretization.cartesian_discretization(nx, ny, r_inf, r_sup)

nth = 30
nr = 30
t_dis, r_dis, th, r, x_s, y_s = polar_discretization.polar_discretization(nth, nr)