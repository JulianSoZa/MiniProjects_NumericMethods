# Import modules:
import meshio
import gmsh
import numpy as np
import os
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

import pyvista as pv

   
def create_initial_conditions(towns, i_towns, points, towns_data, a =1, alpha = 0.8):
    
    num_points = len(points)
    
    S_vec = np.zeros(num_points)
    den = 2*np.pi*np.log( 1/(1-alpha) )
    
    ### Población inicial de Suceptibles
    for town in towns:
        
        area = towns_data[town]["area"]
        xc,yc = towns_data[town]["centroide"]
        pob = towns_data[town]["poblacion"]
        
        sigma = np.sqrt( area/den )
        r2 = (points[:,0]-xc)**2 + (points[:,1]-yc)**2
        S_vec += a*pob*(1/(2*np.pi*sigma**2))*np.exp(-r2/(2*sigma**2))
        
    ### Población inicial de Ifectados
    I_vec = np.zeros(num_points)
    for town in i_towns:
        
        area = towns_data[town]["area"]
        xc,yc = towns_data[town]["centroide"]
        pob = towns_data[town]["poblacion"]
    
        sigma = np.sqrt( area/den )
        r2 = (points[:,0]-xc)**2 + (points[:,1]-yc)**2
        I_vec += a*pob*(1/(2*np.pi*sigma**2))*np.exp(-r2/(2*sigma**2))
    
    return S_vec, I_vec
        
"""if __name__ == "__main__":
    import json
    with open("datos_proyecto/municipios.json", 'r') as openfile:
        towns_data = json.load(openfile)
    towns = [...] ### Lista de municipios que har´an parte de los suceptibles
    i_towns = [...] ### Lista de municipios que har´an parte de los infectados
    ### El array points corresponde a los puntos de la malla ya creada
    S0, I0 = create_initial_conditions(towns, i_towns, points, towns_data, 1, alpha=0.99)"""