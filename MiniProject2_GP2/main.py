import matplotlib.pyplot as plt
from modules import electric_field, electric_potential
from modules.domainDiscretization import cartesian as doCartesian 
from modules.domainDiscretization import polar as doPolar

r_inf = 3
r_sup = 8

nth = 400
nr = 400

t_dis, r_dis, th, r, x_s, y_s, dth, dr, nk = doPolar.polar_discretization(nth, nr)

V, num = electric_potential.electric_potential_solution(nth, nr, dth, dr, t_dis, r_dis, nk, x_s, y_s)

electric_field.electric_field_solution(nth, nr, dth, dr, t_dis, r_dis, nk, V, num, x_s, y_s)

nx = 400
ny = 400

elemento, puntos, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization(nx, ny, r_inf, r_sup)

V = electric_potential.electric_potential_solution_cartesian(elemento, puntos, nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, elementosIndices, nk)

electric_field.electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V, delx, dely, puntosIndices, r_sup, r_inf, elementosIndices, nk)