import meshio
import gmsh
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator as Interpolator
from scipy.interpolate import NearestNDInterpolator as NRInterpolator

def anymesh2D(file, tm, name):
    gmsh.initialize()
    gmsh.clear()

    map = np.load(file)
    for i in range(len(map)):
        gmsh.model.geo.addPoint(map[i][0], map[i][1], 0, tm, i)
    poly = gmsh.model.geo.add_polyline(np.append(np.arange(1, len(map)),1))    #Linea cerrada de los puntos
    curve = gmsh.model.geo.add_curve_loop( [poly])
    s = gmsh.model.geo.add_plane_surface([curve])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.SurfaceFaces",1)
    gmsh.option.setNumber("Mesh.Points",1)
    gmsh.write(f"{name}.msh")
    gmsh.fltk.run()
    gmsh.finalize()

    return

if __name__ == "__main__":
    file = "mapa_antioquia.npy"
    tm = 1e3
    name = 'Antioquia'
    anymesh2D(file, tm, name)
