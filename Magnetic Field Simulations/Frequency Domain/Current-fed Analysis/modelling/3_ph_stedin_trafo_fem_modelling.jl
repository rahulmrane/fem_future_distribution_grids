# #######################################################################################
# FEM Modelling of the STEDIN Transformer (3 Phase)
# Rahul Rane
# #######################################################################################

printstyled(" ▸ Loading Packages .... \r", color = :red)
start = time_ns()
using LinearAlgebra
using fem_future_distribution_grids
elapsed = round((time_ns() - start)/10^9, digits=2)
printstyled(" ✓ Packages loaded ("*string(elapsed)*" seconds)                               \n", color = :green)

function initialise_params(mesh)
    ## Time
    printstyled(" ▸ Setting initial parameters .... \r", color = :red)
    start = time_ns()
    
    S = 400e3                    # Power rating
    Vp = 10750 * sqrt(2)         # Primary peak phase voltage
    Vs = 420 * sqrt(2/3)         # Secondary peak phase voltage
    Ip = (S/10750) * sqrt(2/9)   # Primary peak phase current
    Is = (S/420) * sqrt(2/3)     # Secondary peak phase current
    Np = 266                     # Primary turns
    Ns = 6                       # Secondary turns
    
    omega = 2*pi*50              # Frequency
    
    ## HV winding dimensions (all phases left/right are identical)
    wwhv = 3e-2
    hwhv = 74e-2
    mwhv = 14.75e-2
    Awhv = wwhv * hwhv
    
    ## LV winding dimensions (all phases left/right are identical)
    wwlv = 2e-2
    hwlv = 74e-2
    mwlv = 11.25e-2
    Awlv = wwlv * hwlv
    
    ## Calculate current density in the windings
    Jp = Np * Ip / Awhv
    Js = Ns * Is / Awlv
    
    ## Source current density J
    ## One term for each of the windings, with a positive and negative part
    ## Note the phase shift between the phases
    sourcefunction(group_id, Jp, Js) = 0 * exp(-1im * 2pi/3) * (1 * (group_id==3) - 1 * (group_id==4)) + 
                                        0 * (1 * (group_id==5) - 1 * (group_id==6)) + 
                                        0 * exp(1im * 2pi/3) * (1 * (group_id==7) - 1 * (group_id==8)) + 
                                        Js * exp(-1im * 2pi/3) * (-1 * (group_id==9) + 1 * (group_id==10)) +
                                        Js * (-1 * (group_id==11) + 1 * (group_id==12)) + 
                                        Js * exp(1im * 2pi/3) * (-1 * (group_id==13) + 1 * (group_id==14))
    sourceperelement = sourcefunction.(mesh.e_group, Jp, Js)
    
    ## Relative permeability model
    mu0 = 4e-7 * pi
    mur = 2500                              # Relative permeability of the core
    reluctivityfunction(group_id, mu0, mur) = (1 / mu0) + (1/(mu0*mur) - 1/mu0) * (group_id == 2)
    reluctivityperelement = reluctivityfunction.(mesh.e_group, mu0, mur)
    
    ## Conductivity
    conductivityfunction(group_id) = 0
    conductivityperelement = map(conductivityfunction, mesh.e_group)
    
    ## Time
    elapsed = round((time_ns() - start)/10^9, digits=2)
    printstyled(" ✓ Initial parameters set ("*string(elapsed)*" seconds)                               \n", color = :green)

    return sourceperelement, reluctivityperelement, conductivityperelement, omega
end

## Load and generate mesh data
mesh = fem_future_distribution_grids.Mesh_Data_STEDIN.get_mesh_data_tri_1e()

## Initialise parameters
sourceperelement, reluctivityperelement, conductivityperelement, omega = initialise_params(mesh)

## Calculate the vector potential
u = fem_future_distribution_grids.FEM_Tri_1e.fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega)

## Post-process for magnetic field and current density
Bx, By, B, Hx, Hy, H, mag_energy = fem_future_distribution_grids.Post_Process_Frequency.post_process(mesh, u, reluctivityperelement)

## Save as plots
printstyled(" ▸ Saving Plots .... \r", color = :red)
start = time_ns()
## Contour plot of the magnetic flux density
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(B))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_B", @__DIR__)
## Contour plot of the magnetic field strength
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(H))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_H", @__DIR__)
## Contour plot of the magnetic energy
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(mag_energy))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_energy", @__DIR__)
elapsed = round((time_ns() - start)/10^9, digits=2)
printstyled(" ✓ Plots saved ("*string(elapsed)*" seconds)                               \n", color = :green)

## Save as VTK file for Paraview visualization
printstyled(" ▸ Saving VTK file .... \r", color = :red)
start = time_ns()
fem_future_distribution_grids.Save_VTK.save_vtk(mesh, norm.(u), norm.(B), norm.(H), norm.(mag_energy), reluctivityperelement, "3_ph_stedin_trafo", @__DIR__)
elapsed = round((time_ns() - start)/10^9, digits=2)
printstyled(" ✓ VTK file saved ("*string(elapsed)*" seconds)                               \n", color = :green)

printstyled(" ✓ Code execution completed !!!!", color = :magenta)