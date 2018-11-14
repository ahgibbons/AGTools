import meep as mp
import argparse
import configparser

import numpy as np
import matplotlib.pyplot as plt

import sys

import AGTools.Optics.BraggStackModel as BSM
import AGTools.Optics.MeepLib as mplib
import random
ps_mat = BSM.SellmeierMedium(BSM.ps_params)
si_mat = BSM.SellmeierMedium(BSM.si_params)

## 1D Photonic crystal Model
## distances are in nm

params = mplib.SimParams()
params.cell_width = 2400
params.cell_height = 1000
params.resolution = 0.5
params.dpml = 50
params.time = 20000
params.dt = 250

params.wavelength = 400
params.freq = 1/params.wavelength
params.df = 2*params.freq

params.si_n = 5.070096878062704  ## n at 400 nm
params.ps_n = 1.6285771373334388 ## n at 400 nm
params.air_n = 1
params.other_n = 1

params.ps_thickness = 800
params.si_thickness = 250

params.sph_radius  = 50
params.sph_spacing = 100
params.num_layers  = 4
params.xoffset = 0

old_offsets = [10, 20, 30,15,5,8]
params.offset = [random.randint(0,params.sph_spacing) for i in range(params.num_layers)]

def generate_model(params, args):

    
    cell = mp.Vector3(params.cell_width, params.cell_height, 0)

    spheres = mplib.geo2D_spherical_pc(params.ps_n, params.si_n,
                                       params.ps_thickness, params.si_thickness, params.cell_height, params.num_layers,
                                       params.cell_width/2-200, params.sph_radius,
                                       params.sph_spacing, params.offset, params.xoffset)

    thin_film = mplib.geo2D_thin_film(params.ps_n, params.ps_thickness, params.cell_height,
                                      params.cell_width/2-100)

    pc_ideal = mplib.geo2D_photonic_crystal(params.ps_n, params.other_n,
                                            params.ps_thickness, params.cell_height,
                                            params.si_thickness, params.si_n,
                                            5, params.cell_width/2-200, 50)

    


    geometry = pc_ideal
    
    if args.background:
        geometry = []


    ## Gaussian source
    source_pos = -1*params.cell_width/2 + params.dpml+5
    sources = [mp.Source(mp.GaussianSource(frequency=params.freq, fwidth=params.df),
                         mp.Ez, center=mp.Vector3(source_pos,0,0),
                         size=mp.Vector3(0,params.cell_height,0))]
    
    
    pml_layers = [mp.PML(params.dpml, direction=mp.X)]


    sim = mp.Simulation(cell_size = cell,
                        boundary_layers = pml_layers,
                        geometry = geometry,
                        filename_prefix=args.filename,
                        sources = sources,
                        k_point = mp.Vector3(0,0,0),
                        resolution = params.resolution)

    freg_trans = mp.FluxRegion(center = mp.Vector3(0.5*params.cell_width - params.dpml - 1,0,0),
                         size = mp.Vector3(0, params.cell_height, 0))

    freg_ref    = mp.FluxRegion(center = mp.Vector3(source_pos + 10, 0, 0),
                                size = mp.Vector3(0, params.cell_height), weight=1.0)


    trans_flux = sim.add_flux(params.freq, params.df , 500, freg_trans)
    ref_flux    = sim.add_flux(params.freq, params.df, 500, freg_ref)
    vol = mp.Volume(mp.Vector3(0), size = mp.Vector3(params.cell_width,0,0))

    if not args.reflectflux:
        sim.load_minus_flux("pulse_bg_flux_{:d}".format(params.wavelength),ref_flux)

    if args.outdir:
        print("Output directory: {:s}".format(args.outdir))
        sim.use_output_directory(args.outdir)

    if args.geometry:
        sim.run(mp.at_beginning(mp.output_epsilon),
                until=1)
    else:   
        sim.run(mp.at_beginning(mp.output_epsilon),
                mp.to_appended("ez", mp.at_every(params.dt, mp.output_efield_z)),
                #mp.to_appended("ep", mp.at_every(params.dt, mp.output_dpwr)),
                mp.in_volume(vol, mp.to_appended("ez_slice", mp.at_every(params.dt, mp.output_efield_z))),
                until=params.time)

    if args.reflectflux:
        sim.save_flux("pulse_bg_flux_{:d}".format(params.wavelength),ref_flux)
    
    #sim.save_flux("bg_flux_other", trans_flux)
    sim.display_fluxes(trans_flux, ref_flux)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--thickness", type=float)
    parser.add_argument("-o", "--outdir", type=str)
    parser.add_argument("-w", "--wavelength", type=float)
    parser.add_argument("-f", "--filename", type=str, default=None)
    parser.add_argument("-g", "--geometry", action="store_true")
    parser.add_argument("-b", "--background", action="store_true")
    parser.add_argument("-r", "--reflectflux", action="store_true")
    parser.add_argument("--time", type=int, default=params.time)
    parser.add_argument("--radius", type=float, default = params.sph_radius)
    parser.add_argument("--spacing", type=float, default = params.sph_spacing)
    parser.add_argument("--nlayers", type=float, default = params.num_layers)
    parser.add_argument("--xoffset", type=float, default = params.xoffset)
    args = parser.parse_args()

    
    if args.thickness:
        params.ps_thickness = args.thickness    

    if args.wavelength:
        params.wavelength = args.wavelength
        params.freq = 1/params.wavelength

    params.sph_radius = args.radius
    params.sph_spacing = args.spacing
    params.num_layers = args.nlayers
    params.xoffset = args.xoffset
    params.time = args.time
    

    generate_model(params, args)

    
    mplib.print_settings(params)

    print(args)
