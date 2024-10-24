# #######################################################################################
# Precompile.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

try
    using Gmsh: gmsh
catch
    using gmsh
end

function precompile_tri_1e_frequency()

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
end

function precompile_voltagefed_tri_1e_frequency()

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
        z_length = [z_len_p; z_len_p; z_len_p;
                    z_len_s; z_len_s; z_len_s]
    
        ## Source voltage V
        ## One term for each of the windings, with a positive and negative part
        coil_voltage = [Vp*exp(1im*-2pi/3); Vp; Vp*exp(1im*2pi/3);
                        0; 0; 0]
    
        ## External resistance Rext
        ## One term for each of the windings, with a positive and negative part
        ext_resistance = [Rp/CFFp; Rp/CFFp; Rp/CFFp; 
                          Rs/CFFs; Rs/CFFs; Rs/CFFs]
    
        ## External inductance Lext
        ## One term for each of the windings, with a positive and negative part
        ext_inductance = Lext .* [1; 1; 1; 
                                  1; 1; 1]
        
        ## Source turn density T
        ## One term for each of the windings, with a positive and negative part
        sourcefunction(group_id, Tp, Ts) = [Tp*(1*(group_id==3) - 1*(group_id==4)),
                                            Tp*(1*(group_id==5) - 1*(group_id==6)),
                                            Tp*(1*(group_id==7) - 1*(group_id==8)),
                                            Ts*(-1*(group_id==9) + 1*(group_id==10)),
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
    
        return sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, core_elements, Rp, Rs
    end
    
    ## Load and generate mesh data
    mesh = fem_future_distribution_grids.Mesh_Data_STEDIN.get_mesh_data_tri_1e()
    
    ## Initialise parameters
    sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, core_elements, Rp, Rs = initialise_params(mesh)
    
    ## Calculate the vector potential
    u = fem_future_distribution_grids.FEM_VoltageFed_Tri_1e.fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length)
    
    ## Calculate the vector potential
    max_itr = 2
    u = fem_future_distribution_grids.FEM_VoltageFed_Tri_1e.fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, core_elements, max_itr)
    
    ## Post-process for magnetic field and current density
    Bx, By, B, Hx, Hy, H, mag_energy = fem_future_distribution_grids.Post_Process_Frequency.post_process(mesh, u, reluctivityperelement)
    
    # Calculate source current density
    Jel = fem_future_distribution_grids.Post_Process_Frequency.source_current_density(mesh, u, sourceperelement)
    
    # Calculate winding loss
    Pwindingp, Pwindings = fem_future_distribution_grids.Post_Process_Frequency.winding_loss(mesh, u, Rp, Rs)

    # Calculate core loss
    Pv, Pcore = fem_future_distribution_grids.Post_Process_Frequency.core_loss(mesh, B, z_length)
end

function precompile_transient_tri_1e_time()

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
        n_cycles = 1
        final_time = n_cycles * (2 * pi / omega)
        dt = (final_time - init_time) / (10 * n_cycles)
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
    
        # Find the linear indices of elements that are part of the core (e_group == 2)
        core_elements = findall(x -> x == 2, mesh.e_group)
        
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Initial parameters set ("*string(elapsed)*" seconds)                               \n", color = :green)
    
        return sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, har_content, har_phase, core_elements
    end
    
    ## Load and generate mesh data
    mesh = fem_future_distribution_grids.Mesh_Data_STEDIN.get_mesh_data_tri_1e()
    
    ## Initialise parameters
    sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, har_content, har_phase, core_elements = initialise_params(mesh)

    ## Calculate the vector potential
    u = fem_future_distribution_grids.FEM_Transient_Tri_1e.fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, har_content, har_phase)
    
    ## Calculate the vector potential
    max_itr = 2
    u = fem_future_distribution_grids.FEM_Transient_Tri_1e.fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, har_content, har_phase, core_elements, max_itr)
    
    ## Post-process for magnetic field and current density
    Bx, By, B, Hx, Hy, H, mag_energy = fem_future_distribution_grids.Post_Process_Time.post_process(mesh, u, reluctivityperelement, time_steps)
end

function precompile_transient_voltagefed_tri_1e_time()
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
        
        ## Conductivity
        conductivityfunction(group_id) = 0
        conductivityperelement = map(conductivityfunction, mesh.e_group)
    
        ## Specify time start, end and step
        init_time = 0
        n_cycles = 1
        final_time = n_cycles * (2 * pi / omega)
        dt = (final_time - init_time) / (10 * n_cycles)
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
    
        # Find the linear indices of elements that are part of the core (e_group == 2)
        core_elements = findall(x -> x == 2, mesh.e_group)
        
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Initial parameters set ("*string(elapsed)*" seconds)                               \n", color = :green)
    
        return sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, har_content, har_phase, core_elements, Rp, Rs
    end
    
    ## Load and generate mesh data
    mesh = fem_future_distribution_grids.Mesh_Data_STEDIN.get_mesh_data_tri_1e()
    
    ## Initialise parameters
    sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, har_content, har_phase, core_elements, Rp, Rs = initialise_params(mesh)

    ## Calculate the vector potential
    u = fem_future_distribution_grids.FEM_Transient_VoltageFed_Tri_1e.fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, har_content, har_phase)
    
    ## Calculate the vector potential
    max_itr = 2
    u = fem_future_distribution_grids.FEM_Transient_VoltageFed_Tri_1e.fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, har_content, har_phase, core_elements, max_itr)
    
    ## Post-process for magnetic field and current density
    Bx, By, B, Hx, Hy, H, mag_energy = fem_future_distribution_grids.Post_Process_Time.post_process(mesh, u, reluctivityperelement, time_steps)
    
    ## Calculate source current density
    Jel = fem_future_distribution_grids.Post_Process_Time.source_current_density(mesh, u, sourceperelement, time_steps)
    
    ## Calculate core loss
    Pv, Pcore = fem_future_distribution_grids.Post_Process_Time.core_loss(mesh, B, z_length, time_steps)

    ## Calculate winding loss
    Pwindingp, Pwindings = fem_future_distribution_grids.Post_Process_Time.winding_loss(mesh, u, Rp, Rs, time_steps)
end

## Suppress output
original_stdout = stdout
redirect_stdout(devnull)

## Run all above functions
precompile_tri_1e_frequency()
precompile_voltagefed_tri_1e_frequency()
precompile_transient_tri_1e_time()
precompile_transient_voltagefed_tri_1e_time()

## Release output
redirect_stdout(original_stdout)