# #######################################################################################
# FEM Modelling of the STEDIN Transformer (3 Phase) [Open Circuit Test]
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

    ## Calculate turn density in the windings
    Tp = Np / Awhv
    Ts = Ns / Awlv

    ## External resistance Rext
    Rp = 1.8131 + 600
    Rs = 1.2999e-3 + 0.2
    CFFp = 0.3
    CFFs = 0.3
    
    ## External Inductance Lext
    Lext = 1e-6
    
    ## Z-axis length
    z_len_p = 0.4
    z_len_s = 0.4
    z_length = [z_len_s; z_len_s; z_len_s]

    ## Source voltage V
    ## One term for each of the windings, with a positive and negative part
    coil_voltage = [Vs*exp(1im*-2pi/3); Vs; Vs*exp(1im*2pi/3)]

    ## External resistance Rext
    ## One term for each of the windings, with a positive and negative part
    ext_resistance = [Rs/CFFs; Rs/CFFs; Rs/CFFs]

    ## External inductance Lext
    ## One term for each of the windings, with a positive and negative part
    ext_inductance = Lext .* [1; 1; 1]
    
    ## Source turn density T
    ## One term for each of the windings, with a positive and negative part
    sourcefunction(group_id, Tp, Ts) = [Ts*(-1*(group_id==9) + 1*(group_id==10)),
                                        Ts*(-1*(group_id==11) + 1*(group_id==12)),
                                        Ts*(-1*(group_id==13) + 1*(group_id==14))]
    sourceperelement = sourcefunction.(mesh.e_group, Tp, Ts)
    
    ## Relative permeability model
    mu0 = 4e-7 * pi
    mur = 2500                              # Relative permeability of the core
    reluctivityfunction(group_id, mu0, mur) = (1 / mu0) + (1/(mu0*mur) - 1/mu0) * (group_id == 2)
    reluctivityperelement = reluctivityfunction.(mesh.e_group, mu0, mur)
    
    ## Conductivity
    conductivityfunction(group_id) = 0
    conductivityperelement = map(conductivityfunction, mesh.e_group)

    # Find the linear indices of elements that are part of the core (e_group == 2)
    core_elements = findall(x -> x == 2, mesh.e_group)
    
    ## Time
    elapsed = round((time_ns() - start)/10^9, digits=2)
    printstyled(" ✓ Initial parameters set ("*string(elapsed)*" seconds)                               \n", color = :green)

    return sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, core_elements
end

## Load and generate mesh data
mesh = fem_future_distribution_grids.Mesh_Data_STEDIN.get_mesh_data_tri_1e()

## Initialise parameters
sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, core_elements = initialise_params(mesh)

## Calculate the vector potential
max_itr = 1e5
u = fem_future_distribution_grids.FEM_VoltageFed_Tri_1e.fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, core_elements, max_itr)

## Post-process for magnetic field and current density
Bx, By, B, Hx, Hy, H, mag_energy = fem_future_distribution_grids.Post_Process_Frequency.post_process(mesh, u, reluctivityperelement)

# Calculate source current density
Jel = fem_future_distribution_grids.Post_Process_Frequency.source_current_density(mesh, u, sourceperelement)

# Calculate core loss
Pv, Pcore = fem_future_distribution_grids.Post_Process_Frequency.core_loss(mesh, B, z_length)
printstyled("Core loss = "*string(Pcore)*" W\n", color = :blue)

## Save as plots
printstyled(" ▸ Saving Plots .... \r", color = :red)
start = time_ns()
## Contour plot of the magnetic flux density
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(B))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_B_OC", @__DIR__)
## Contour plot of the magnetic field strength
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(H))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_H_OC", @__DIR__)
## Contour plot of the magnetic energy
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(mag_energy))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_energy_OC", @__DIR__)
# Contour plot of the source current density
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(Jel))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_Jel_OC", @__DIR__)
# Contour plot of the Pv
msh, colors = fem_future_distribution_grids.Mesh_Data_STEDIN.get_cell_mesh_tri_1e(norm.(Pv))
fem_future_distribution_grids.Makie_Plots.plot_surface_tri_1e(msh, colors, "3_ph_stedin_trafo_Pv_OC", @__DIR__)
elapsed = round((time_ns() - start)/10^9, digits=2)
printstyled(" ✓ Plots saved ("*string(elapsed)*" seconds)                               \n", color = :green)

## Save as VTK file for Paraview visualization
printstyled(" ▸ Saving VTK file .... \r", color = :red)
start = time_ns()
fem_future_distribution_grids.Save_VTK.save_vtk(mesh, imag(u[1:mesh.nnodes]), norm.(B), norm.(H), norm.(mag_energy), reluctivityperelement, norm.(Jel), norm.(Pv), "3_ph_stedin_trafo_OC", @__DIR__)
elapsed = round((time_ns() - start)/10^9, digits=2)
printstyled(" ✓ VTK file saved ("*string(elapsed)*" seconds)                               \n", color = :green)

printstyled(" ✓ Code execution completed !!!!", color = :magenta)