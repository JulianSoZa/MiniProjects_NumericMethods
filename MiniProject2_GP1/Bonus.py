import matplotlib.pyplot as plt
import pandas as pd
from modules import electric_field, electric_potential, plot_electric_solution
from modules.domainDiscretization import cartesian as doCartesian 
from modules.domainDiscretization import polar as doPolar
import analysis_results
import numpy as np
from collections import OrderedDict
from matplotlib.patches import Polygon

r_inf = 3
r_sup = 8

nx = 10
ny = 10

nk = nx*ny

p = [1,2,3,4]     #puntos dentro de los limites 1<p<4

for i in p:
    elemento, coordenada, x_dis, y_dis, delx, dely, puntosIndices, elementosIndices, nk = doCartesian.cartesian_discretization_Bonus(nx, ny, r_inf, r_sup, i)

    V_c, V_space = electric_potential.electric_potential_solution_cartesian(nx, ny, x_dis, y_dis, delx, dely, r_inf, r_sup, puntosIndices, nk)
    plot_electric_solution.ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, V_c,  'Potencial')

    En_c, E_space = electric_field.electric_field_solution_cartesian(nx, ny, x_dis, y_dis, V_c, delx, dely, puntosIndices, r_sup, r_inf, nk)
    plot_electric_solution.ploter_finite_solutions_cartesian(nk, nx, x_dis, y_dis, puntosIndices, elementosIndices, En_c, 'Campo')



