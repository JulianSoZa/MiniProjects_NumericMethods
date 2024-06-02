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
    polys = []
    
    map = np.load(file)
    for i in range(len(map)):
        gmsh.model.geo.addPoint(map[i][0], map[i][1], 0, tm)
        if i > 0:
            poly = gmsh.model.geo.addLine(i, i+1)
            polys.append(poly)
    
    polys.append(gmsh.model.geo.addLine(len(map), 1))
    #poly = gmsh.model.geo.add_polyline(np.append(np.arange(1, len(map)),1))    #Linea cerrada de los puntos
    curve = gmsh.model.geo.add_curve_loop(polys)
    s = gmsh.model.geo.add_plane_surface([curve])
    gmsh.model.mesh.embed(1, [poly], 1, s)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.SurfaceFaces",1)
    gmsh.option.setNumber("Mesh.Points",1)
    gmsh.write(f"MiniProject3_GP1/data/meshes/{name}.msh")
    #gmsh.fltk.run()
    gmsh.finalize()

    return polys  

if __name__ == "__main__":
    try: 
        print('Inicia la lectura\n')
        malla = meshio.read("MiniProject3_GP1\data\meshes\Antioquia.msh")
        print('Se leyó correctamente\n')
    except:
        print('El archivo no estaba creado \n')
        file = "MiniProject3_GP1/data/meshes/mapa_antioquia.npy"
        tm = 1e3                        #tamaño promedio del elemento
        name = 'Antioquia'
        polys = anymesh2D(file, tm, name)
        print('Archivo creado \n')
