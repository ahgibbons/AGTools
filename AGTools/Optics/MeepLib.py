import meep as mp
import numpy as np

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

def geo2D_spherical_pc(n_matrix, n_substrate, film_thickness, film_height,
                        substrate_thickness, n_layers, substrate_base,
                        sph_radius, sph_spacing, layer_offset):
    """
    Generate a layered structure with spheres. Within a layer the spheres are hexagonally packed
    """
    matrix_block = mp .Block(size=mp.Vector3(film_thickness, film_height, 1e20),
                            center = mp.Vector3(substrate_base-substrate_thickness-film_thickness/2),
                            material = mp.Medium(epsilon = n_matrix**2))

    sub_block = mp.Block(size = mp.Vector3(substrate_thickness, film_height, 1e20),
                         center = mp.Vector3(substrate_base-substrate_thickness/2,0,0),
                         material = n_substrate**2)

    spheres = [matrix_block]
    #spheres = []
    
    for n in range(n_layers):
        x = substrate_base - substrate_thickness - film_thickness*(n/n_layers + 1/12)
        offset = layer_offset[n]
        num_spheres = int(film_height / sph_spacing)
        sphere_layer = [mp.Sphere(radius=sph_radius, center=mp.Vector3(x, y+offset - film_height/2, 0),
                                  material=mp.Medium(epsilon=1)) for y in np.linspace(0, film_height,
                                                                                      num_spheres)]
        spheres = spheres + sphere_layer

    return spheres 

t_check = "Hello"
    
class SimParams:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
