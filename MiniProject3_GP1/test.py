import meshio
import gmsh
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator as Interpolator
from scipy.interpolate import NearestNDInterpolator as NRInterpolator
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


name = 'Antioquia'
mesh = meshio.read(f"{name}.msh")
niveles = np.load('niveles.npy')


vertices = mesh.points
num_vertices = len(vertices)
triangles = mesh.cells_dict['triangle']

xy = niveles[:,:2]
z = niveles[:,2]

interpolator_fun = Interpolator(xy, z)
nr_interpolator_fun = NRInterpolator(xy,z)

print('Vamos k ')
mesh_level = interpolator_fun(vertices[:,0],vertices[:,1])
nans_values = np.isnan(mesh_level)
mesh_level[nans_values] = nr_interpolator_fun(vertices[nans_values,0],vertices[nans_values,1])



