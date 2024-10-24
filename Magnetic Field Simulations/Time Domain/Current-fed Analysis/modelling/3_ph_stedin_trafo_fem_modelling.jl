# #######################################################################################
# FEM Modelling of the STEDIN Transformer (3 Phase)
# Rahul Rane
# #######################################################################################

printstyled(" ▸ Loading Packages .... \r", color = :red)
start = time_ns()
using LinearAlgebra
using Plots
using DataFrames
using CSV
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
    
    ## Conductivity
    conductivityfunction(group_id) = 0
    conductivityperelement = map(conductivityfunction, mesh.e_group)

    ## Specify time start, end and step
    init_time = 0
    n_cycles = 5
    final_time = n_cycles * (2 * pi / omega)
    dt = (final_time - init_time) / (60 * n_cycles)
    time_steps = Vector(init_time:dt:final_time)
    
    ## Relative permeability model
    mu0 = 4e-7 * pi
    mur = 2500                              # Relative permeability of the core
    reluctivityfunction(group_id, mu0, mur) = (1 / mu0) + (1/(mu0*mur) - 1/mu0) * (group_id == 2)
    reluctivityperelementstep = reluctivityfunction.(mesh.e_group, mu0, mur)
    reluctivityperelement = Vector{Vector{Float64}}(undef, length(time_steps))
    for k in eachindex(time_steps)
        reluctivityperelement[k] = reluctivityperelementstep
    end

    ## Harmonic content
    har_content = [1]
    har_phase = [0]
    # # Third harmonic addition
    # har_content = [1, 3]
    # har_phase = [0, pi/3]
    
    ## Time
    elapsed = round((time_ns() - start)/10^9, digits=2)
    printstyled(" ✓ Initial parameters set ("*string(elapsed)*" seconds)                               \n", color = :green)

    return sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, har_content, har_phase
end

## Load and generate mesh data
mesh = fem_future_distribution_grids.Mesh_Data_STEDIN.get_mesh_data_tri_1e()

## Initialise parameters
sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, har_content, har_phase = initialise_params(mesh)

## Calculate the vector potential
u = fem_future_distribution_grids.FEM_Transient_Tri_1e.fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, har_content, har_phase)

## Post-process for magnetic field and current density
Bx, By, B, Hx, Hy, H, mag_energy = fem_future_distribution_grids.Post_Process_Time.post_process(mesh, u, reluctivityperelement, time_steps)

## Plots
selected_node = 4276
selected_element = 7051

# Select single element or node
u_wave = [u_curr[selected_node] for u_curr in u]
Bx_wave = [Bx_curr[selected_element] for Bx_curr in Bx]
By_wave = [By_curr[selected_element] for By_curr in By]
B_wave = [B_curr[selected_element] for B_curr in B]
Hx_wave = [Hx_curr[selected_element] for Hx_curr in Hx]
Hy_wave = [Hy_curr[selected_element] for Hy_curr in Hy]
H_wave = [H_curr[selected_element] for H_curr in H]
mag_energy_wave = [mag_energy_curr[selected_element] for mag_energy_curr in mag_energy]
reluctivityperelement_wave = [reluctivityperelement_curr[selected_element] for reluctivityperelement_curr in reluctivityperelement]

plot(time_steps, u_wave, label = "u")
xlabel!("time (sec)")
ylabel!("u")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_u.png")

plot(time_steps, Bx_wave, label = "Bx")
xlabel!("time (sec)")
ylabel!("Bx")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_Bx.png")

plot(time_steps, By_wave, label = "By")
xlabel!("time (sec)")
ylabel!("By")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_By.png")

plot(time_steps, B_wave, label = "u")
xlabel!("time (sec)")
ylabel!("B")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_B.png")

plot(time_steps, Hx_wave, label = "Hx")
xlabel!("time (sec)")
ylabel!("Hx")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_Hx.png")

plot(time_steps, Hy_wave, label = "Hy")
xlabel!("time (sec)")
ylabel!("Hy")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_Hy.png")

plot(time_steps, H_wave, label = "H")
xlabel!("time (sec)")
ylabel!("H")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_H.png")

plot(time_steps, mag_energy_wave, label = "energy")
xlabel!("time (sec)")
ylabel!("energy")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_energy.png")

plot(time_steps, reluctivityperelement_wave, label = "mur")
xlabel!("time (sec)")
ylabel!("mur")
savefig((@__DIR__) * "/../img/3_ph_stedin_trafo_mur.png")

# printstyled(" ▸ Saving CSV files .... \r", color = :red)
# start = time_ns()
# Bx_mat = DataFrame(hcat(Bx...), :auto)
# CSV.write((@__DIR__) * "/../csv/3_ph_stedin_trafo_Bx.csv", Bx_mat)
# By_mat = DataFrame(hcat(By...), :auto)
# CSV.write((@__DIR__) * "/../csv/3_ph_stedin_trafo_By.csv", By_mat)
# B_mat = DataFrame(hcat(B...), :auto)
# CSV.write((@__DIR__) * "/../csv/3_ph_stedin_trafo_B.csv", B_mat)
# Hx_mat = DataFrame(hcat(Hx...), :auto)
# CSV.write((@__DIR__) * "/../csv/3_ph_stedin_trafo_Hx.csv", Hx_mat)
# Hy_mat = DataFrame(hcat(Hy...), :auto)
# CSV.write((@__DIR__) * "/../csv/3_ph_stedin_trafo_Hy.csv", Hy_mat)
# H_mat = DataFrame(hcat(H...), :auto)
# CSV.write((@__DIR__) * "/../csv/3_ph_stedin_trafo_H.csv", H_mat)
# mag_energy_mat = DataFrame(hcat(mag_energy...), :auto)
# CSV.write((@__DIR__) * "/../csv/3_ph_stedin_trafo_energy.csv", mag_energy_mat)
# elapsed = round((time_ns() - start)/10^9, digits=2)
# printstyled(" ✓ CSV files saved ("*string(elapsed)*" seconds)                               \n", color = :green)

printstyled(" ✓ Code execution completed !!!!", color = :magenta)