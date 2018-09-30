import meep as mp

def geo1D_thin_film_si(film_n, si_n, ps_thickness, si_thickness, pos=0):
    ps_block = mp.Block(size=mp.Vector3(1,1,ps_thickness),
                               center = mp.Vector3(0,0,pos - (si_thickness + ps_thickness/2 )),
                               material = mp.Medium(epsilon = film_n**2))

    si_block = mp.Block(size=mp.Vector3(1,1,si_thickness),
                        center = mp.Vector3(0,0,pos - si_thickness/2),
                        material = mp.Medium(epsilon = si_n**2))

    
    return [ps_block, si_block]