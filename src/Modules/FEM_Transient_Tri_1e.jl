# #######################################################################################
# FEM_Transient_Tri_1e.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module FEM_Transient_Tri_1e

    using SparseArrays
    using LinearAlgebra
    using Statistics
    using fem_future_distribution_grids

    function fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, num_harmonic, phase_diff)
        ## Initialise vectors
        u = Vector{Vector{Float64}}(undef, length(time_steps))
        S = zeros(Float64, mesh.nnodes)
    
        ## Time step
        dt = time_steps[2] - time_steps[1]

        ## Assemble A, B, and f Matrices
        A, B, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement[1], conductivityperelement, omega)
        C = B .+ dt .* A

        ## Factorize the matrix for efficiency
        C_factored = factorize(C)

        ## Initial condition
        u[1] = zeros(Float64, mesh.nnodes)

        ## Time
        printstyled(" ▸ Computing solution .... \r", color = :red)
        start = time_ns()

        ## Backward Euler
        for k = firstindex(time_steps) + 1 : lastindex(time_steps)
            ## Initialise voltage source
            fill!(S, 0.0)

            ## Sum harmonic contributions
            for e in eachindex(num_harmonic)
                @. S += imag(f / num_harmonic[e] * exp(num_harmonic[e] * 1im * omega * time_steps[k] + phase_diff[e]))
            end
            
            ## Compute the solution at the next time step
            u[k] = copy(u[k-1])
            mul!(u[k], B, u[k-1])
            @. u[k] += dt * S
            u[k] = C_factored \ u[k]
            
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
    
    function fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, time_steps, num_harmonic, phase_diff, core_elements, max_itr)
        ## Threshold value for the error
        threshold = 1e-2
        alpha = 0.9
    
        ## Initialise vectors
        u = Vector{Vector{Float64}}(undef, length(time_steps))
        S = zeros(Float64, mesh.nnodes)
    
        ## Time step
        dt = time_steps[2] - time_steps[1]

        ## Initial solution
        A, C, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement[1], conductivityperelement, omega)
        K = C .+ dt .* A
    
        ## Initial condition
        u[1] = zeros(Float64, mesh.nnodes)

        ## Preallocate buffers for reuse
        u_temp = zeros(Float64, mesh.nnodes)
        u_prev = zeros(Float64, mesh.nnodes)
        u_hist = zeros(Float64, mesh.nnodes)

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
                @. K = C + dt * A

                ## Factorize the matrix for efficiency
                K_factored = factorize(K)

                ## Initialise voltage source
                fill!(S, 0.0)

                ## Sum harmonic contributions
                for e in eachindex(num_harmonic)
                    @. S += imag(f / num_harmonic[e] * exp(num_harmonic[e] * 1im * omega * time_steps[k] + phase_diff[e]))
                end

                ## Compute the solution at the next time step
                mul!(u_temp, C, u_hist)
                @. u_temp += dt * S
                u_temp = K_factored \ u_temp
            
                ## Check the error against the threshold
                epsilon = 1e-12
                avg_percentage_error = mean(abs.((u_temp .- u_prev) ./ (u_prev .+ epsilon)))
                if avg_percentage_error <= threshold
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

        ## Time 
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Solution computed ("*string(elapsed)*" seconds)                                                              \n", color = :green)

        return u
    end

end