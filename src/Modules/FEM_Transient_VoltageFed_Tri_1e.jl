# #######################################################################################
# FEM_Transient_VoltageFed_Tri_1e.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module FEM_Transient_VoltageFed_Tri_1e

    using SparseArrays
    using LinearAlgebra
    using fem_future_distribution_grids

    function fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, num_harmonic, phase_diff)
        ## Initialise vectors
        u = Vector{Vector{Float64}}(undef, length(time_steps))
        S = zeros(Float64, length(coil_voltage))
    
        ## Time step
        dt = time_steps[2] - time_steps[1]

        ## Assemble A, B, and f Matrices
        A, B, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement[1], conductivityperelement, omega)

        ## Time
        printstyled(" ▸ Modelling circuit equations .... \r", color = :red)
        start = time_ns()

        ## Circuit Equations
        K11 = B .+ dt .* A
        K12 = - dt .* f
        K21 = transpose(f) .* z_length
        K22 = (dt .* Diagonal(ext_resistance)) .+ Diagonal(ext_inductance)
        K = [K11 K12; K21 K22]
        
        M11 = B
        M12 = zeros(Float64, mesh.nnodes,length(coil_voltage))
        M21 = transpose(f) .* z_length
        M22 = Diagonal(ext_inductance)
        M = [M11 M12; M21 M22]

        ## RHS vector
        T = zeros(Float64, mesh.nnodes+length(coil_voltage))

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Circuit equations modelled ("*string(elapsed)*" seconds)                               \n", color = :green)
    
        ## Factorize the matrix for efficiency
        K_factored = factorize(K)

        ## Initial condition
        u[1] = zeros(Float64, mesh.nnodes + length(coil_voltage))

        ## Time
        printstyled(" ▸ Computing solution .... \r", color = :red)
        start = time_ns()
        
        for k = firstindex(time_steps) + 1 : lastindex(time_steps)
            ## Initialise voltage source
            fill!(S, 0.0)

            ## Sum harmonic contributions
            for e in eachindex(num_harmonic)
                @. S += imag(coil_voltage / num_harmonic[e] * exp(num_harmonic[e] * 1im * omega * time_steps[k] + phase_diff[e]))
            end

            ## RHS vector
            T[mesh.nnodes+1 : mesh.nnodes+length(coil_voltage)] = dt .* S

            ## Compute the solution at the next time step
            u[k] = copy(u[k-1])
            mul!(u[k], M, u[k-1])
            @. u[k] += T
            u[k] = K_factored \ u[k]
            
            ## Time
            if mod(k, 10) == 0
                progress = round(k / length(time_steps) * 100, digits=1)
                elapsed = round((time_ns() - start) / 1e9, digits=2)
                estimated = round(length(time_steps) * elapsed / k, digits=2)
                print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            end
        end

        ## Time 
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Solution computed ("*string(elapsed)*" seconds)                               \n", color = :green)
    
        return u
    end

    function mu_func(B)
        k1 = 3.8
        k2 = 2.17
        k3 = 396.2

        v = k1 * exp(k2 * B^2) + k3
        
        return (1 ./ v)
    end
        
    function fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, num_harmonic, phase_diff, core_elements, max_itr)
        ## Threshold value for the error
        threshold = 1e-2
        alpha = 0.9
    
        ## Initialise vectors
        u = Vector{Vector{Float64}}(undef, length(time_steps))
        S = zeros(Float64, length(coil_voltage))
    
        ## Time step
        dt = time_steps[2] - time_steps[1]

        ## Initial solution
        A, C, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement[1], conductivityperelement, omega)

        ## Time
        printstyled(" ▸ Modelling circuit equations .... \r", color = :red)
        start = time_ns()

        ## Circuit Equations
        K11 = C .+ dt .* A
        K12 = - dt .* f
        K21 = transpose(f) .* z_length
        K22 = (dt .* Diagonal(ext_resistance)) .+ Diagonal(ext_inductance)
        K = [K11 K12; K21 K22]
        
        M11 = C
        M12 = zeros(Float64, mesh.nnodes,length(coil_voltage))
        M21 = transpose(f) .* z_length
        M22 = Diagonal(ext_inductance)
        M = [M11 M12; M21 M22]

        ## RHS vector
        T = zeros(Float64, mesh.nnodes+length(coil_voltage))

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Circuit equations modelled ("*string(elapsed)*" seconds)                               \n", color = :green)

        ## Initial condition
        u[1] = zeros(Float64, mesh.nnodes + length(coil_voltage))

        ## Preallocate buffers for reuse
        u_temp = zeros(Float64, mesh.nnodes + length(coil_voltage))
        u_prev = zeros(Float64, mesh.nnodes + length(coil_voltage))
        u_hist = zeros(Float64, mesh.nnodes + length(coil_voltage))

        ## Relucitvity at initial time step
        reluctivity = copy(reluctivityperelement[1])

        ## Time
        printstyled(" ▸ Computing solution .... \r", color = :red)
        start = time_ns()
        
        for k = firstindex(time_steps) + 1 : lastindex(time_steps)
            ## Save a copy of the previous time step solution
            copy!(u_hist, u_temp)

            itr_count = max_itr
            for loop = 1:max_itr
                ## Save a copy of the previous solution
                copy!(u_prev, u_temp)
                
                ## Update history with a blend of previous and current solution
                @. u_hist = alpha * u_hist + (1 - alpha) * u_temp
                
                ## Compute B norm
                B = fem_future_distribution_grids.Post_Process_Frequency.B_norm(mesh, u_hist)

                ## Assign new value of mur in core elements
                @. reluctivity[core_elements] = 1.0 / mu_func(abs(B[core_elements]))

                ## Update solution
                A = fem_future_distribution_grids.Assemble_Matrices.assemble_A(mesh, reluctivity)
                @. K11 = C + dt * A
                K[1:mesh.nnodes, 1:mesh.nnodes] = K11
    
                ## Factorize the matrix for efficiency
                K_factored = factorize(K)

                ## Initialise voltage source
                fill!(S, 0.0)

                ## Sum harmonic contributions
                for e in eachindex(num_harmonic)
                    @. S += imag(coil_voltage / num_harmonic[e] * exp(num_harmonic[e] * 1im * omega * time_steps[k] + phase_diff[e]))
                end

                ## RHS vector
                T[mesh.nnodes+1 : mesh.nnodes+length(coil_voltage)] = dt .* S

                ## Compute the solution at the next time step
                mul!(u_temp, M, u_hist)
                @. u_temp += T
                u_temp = K_factored \ u_temp
            
                ## Check the error against the threshold
                if norm(u_temp - u_prev) <= threshold
                    itr_count = loop
                    break
                end
            end

            ## Save final values at current time step
            u[k] = copy(u_temp)
            reluctivityperelement[k] = copy(reluctivity)
            
            ## Time
            if mod(k, 1) == 0
                progress = round(k / length(time_steps) * 100, digits=1)
                elapsed = round((time_ns() - start) / 1e9, digits=2)
                estimated = round(length(time_steps) * elapsed / k, digits=2)
                print(" ▸ Computing nonlinear solution (Iteration "*string(Int(itr_count))*") .... " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            end
        end
                                                
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")

        return u
    end

end