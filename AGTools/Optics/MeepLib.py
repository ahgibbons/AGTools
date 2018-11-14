import meep as mp
import numpy as np
import pandas as pd

def geo1D_thin_film_si(film_n, si_n, ps_thickness, si_thickness, pos=0, bare_substrate=False):
    ps_block = mp.Block(size=mp.Vector3(1,1,ps_thickness),
                               center = mp.Vector3(0,0,pos - (si_thickness + ps_thickness/2 )),
                               material = mp.Medium(epsilon = film_n**2))

    si_block = mp.Block(size=mp.Vector3(1,1,si_thickness),
                        center = mp.Vector3(0,0,pos - si_thickness/2),
                        material = mp.Medium(epsilon = si_n**2))

    if bare_substrate:
        geometry = [si_block]
    else:
        geometry = [ps_block, si_block]
    
    return geometry

def geo2D_thin_film(n_matrix, film_thickness, film_height,
                    film_base_pos):
    
    matrix_block = mp.Block(size=mp.Vector3(film_thickness, film_height, 1e20),
                            center = mp.Vector3(film_base_pos - film_thickness/2, 0, 0),
                            material = mp.Medium(epsilon = n_matrix**2))

    return [matrix_block]

def geo2D_photonic_crystal(n_matrix, n_other, film_thickness, film_height, si_thickness, n_si,
                           n_layers, film_base_pos, layer_d):
    """
    Generate a layered structure with alternating layers. Design is to make equivalent structure to 
    'geo2D_spherical_pc'. 
    """

    matrix_block = mp.Block(size=mp.Vector3(film_thickness, film_height, 1e20),
                            center = mp.Vector3(film_base_pos-film_thickness/2),
                            material = mp.Medium(epsilon = n_matrix**2))

    si_block = mp.Block(size=mp.Vector3(si_thickness, film_height, 1e20),
                        center = mp.Vector3(film_base_pos+si_thickness/2),
                        material = mp.Medium(epsilon = n_si**2))

    layers = [si_block, matrix_block]

    for n in range(n_layers):
        x = film_base_pos - film_thickness*(n + 0.5)/n_layers
        pc_layer = mp.Block(size=mp.Vector3(layer_d, film_height, 1e20),
                            center = mp.Vector3(x, 0, 0),
                            material = mp.Medium(epsilon = n_other**2))

        layers.append(pc_layer)

    return layers



    
def geo3D_photonic_crystal(n_matrix, n_other, film_thickness, filmx, filmy, si_thickness, n_si,
                           n_layers, film_base, layer_d):
    """
    Generate a layered structure with alternating layers. Design is to make equivalent structure to 
    'geo2D_spherical_pc'. 
    """

    matrix_block = mp.Block(size=mp.Vector3(filmx, filmy, film_thickness),
                            center = mp.Vector3(0,0,film_base-film_thickness/2),
                            material = mp.Medium(epsilon = n_matrix**2))

    si_block = mp.Block(size=mp.Vector3(filmx, filmy, si_thickness),
                        center = mp.Vector3(0,0,film_base+si_thickness/2),
                        material = mp.Medium(epsilon = n_si**2))

    layers = [si_block, matrix_block]

    for n in range(n_layers):
        z = film_base - film_thickness*(n + 0.5)/n_layers
        pc_layer = mp.Block(size=mp.Vector3(filmx, filmy, layer_d),
                            center = mp.Vector3(0,0,z),
                            material = mp.Medium(epsilon = n_other**2))

        layers.append(pc_layer)

    return layers

def geo2D_ellipsoid_pc(n_matrix, n_substrate, film_xsize, substrate_xsize, film_ysize,
                        n_layers, film_base_x, 
                        ellipsex, ellipsey, ellipse_spacing, layer_offset, ellipse_x_offset=0 ):
    """
    Generate a layered structure with spheres. Within a layer the spheres are hexagonally packed
    """
    matrix_block = mp.Block(size=mp.Vector3(film_xsize, film_ysize, 1e20),
                            center = mp.Vector3(film_base_x-film_xsize/2),
                            material = mp.Medium(epsilon = n_matrix**2))

    si_block = mp.Block(size=mp.Vector3(substrate_xsize, film_ysize, 1e20),
                        center = mp.Vector3(film_base_x+substrate_xsize/2,0,0),
                        material = mp.Medium(epsilon = n_substrate**2))

    print(n_substrate)
    
    spheres = [si_block, matrix_block]
    #spheres = []
    
    for n in range(n_layers):
        x = film_base_x - film_xsize*(n + 0.5)/n_layers + ellipse_x_offset
        offset = layer_offset[n]
        num_spheres = int(film_ysize / ellipse_spacing)
        sphere_layer = [mp.Ellipsoid(size=mp.Vector3(ellipsex, ellipsey, 1e20), center=mp.Vector3(x, y+offset - film_ysize/2, 0),
                                  material=mp.Medium(epsilon=1)) for y in np.linspace(0, film_ysize,
                                                                                      num_spheres)]
        spheres = spheres + sphere_layer

    return spheres

    
def geo2D_spherical_pc(n_matrix, n_substrate, film_xsize, substrate_xsize, film_ysize,
                        n_layers, film_base_x, 
                        sph_radius, sph_spacing, layer_offset, sph_x_offset=0 ):
    """
    Generate a layered structure with spheres. Within a layer the spheres are hexagonally packed
    """
    matrix_block = mp.Block(size=mp.Vector3(film_xsize, film_ysize, 1e20),
                            center = mp.Vector3(film_base_x-film_xsize/2),
                            material = mp.Medium(epsilon = n_matrix**2))

    si_block = mp.Block(size=mp.Vector3(substrate_xsize, film_ysize, 1e20),
                        center = mp.Vector3(film_base_x+substrate_xsize/2,0,0),
                        material = mp.Medium(epsilon = n_substrate**2))

    print(n_substrate)
    
    spheres = [si_block, matrix_block]
    #spheres = []
    
    for n in range(n_layers):
        x = film_base_x - film_xsize*(n + 0.5)/n_layers + sph_x_offset
        offset = layer_offset[n]
        num_spheres = int(film_ysize / sph_spacing)
        sphere_layer = [mp.Sphere(radius=sph_radius, center=mp.Vector3(x, y+offset - film_ysize/2, 0),
                                  material=mp.Medium(epsilon=1)) for y in np.linspace(0, film_ysize,
                                                                                      num_spheres)]
        spheres = spheres + sphere_layer

    return spheres

def print_settings(params):

    print("Simulation User Parameters\n\n")

    for k in params.__dict__.keys():
        print("{:s}:\t{:s}".format(k, str(params.__dict__[k])))
    """
    print("")
    print("CELL_WIDTH:\t\t{:d}".format(params.cell_width))
    print("CELL_HEIGHT:\t\t{:d}".format(params.cell_height))
    print("RESOLUTION:\t\t{:f}".format(params.resolution))
    print("DPML:\t\t{:d}".format(params.dpml))
    print("TIME:\t\t{:d}".format(params.time))
    print("DT:\t\t{:d}".format(params.dt))
    print("Frequency:\t\t{:f}".format(params.freq))
    print("Wavelength:\t\t{:f}".format(params.wavelength))

    print("PS Thickness:\t\t{:f}".format(params.ps_thickness))
    print("Si Thickness:\t\t{:f}".format(params.si_thickness))
    """

    
class SimParams:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
