# #######################################################################################
# FEM_VoltageFed_Tri_1e.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module FEM_VoltageFed_Tri_1e

    using LinearAlgebra
    using Statistics
    using fem_future_distribution_grids

    function fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length)
        ## Assemble A, B, and f Matrices
        A, B, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega)

        ## Time
        printstyled(" ▸ Modelling circuit equations .... \r", color = :red)
        start = time_ns()

        ## Circuit Equations
        K11 = A .+ 1im.*B
        K12 = -f
        K21 = (1im * omega * z_length) .* transpose(f)
        K22 = Diagonal(ext_resistance) .+ (1im * omega) .* Diagonal(ext_inductance)
        K = [K11 K12; K21 K22]

        ## RHS vector
        T = [zeros(Float64, mesh.nnodes); coil_voltage]

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Circuit equations modelled ("*string(elapsed)*" seconds)                               \n", color = :green)

        ## Time
        printstyled(" ▸ Computing solution .... \r", color = :red)
        start = time_ns()

        ## Factorize the matrix for efficiency
        K_factored = factorize(K)

        ## Compute the numerical solution
        u = K_factored \ T
        
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

    function fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, coil_voltage, ext_resistance, ext_inductance, z_length, core_elements, max_itr)        
        ## Threshold value for the error
        threshold = 1e-2
        alpha = 0.9

        ## Initial solution
        ## Assemble A, B, and f Matrices
        A, C, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega)

        ## Time
        printstyled(" ▸ Modelling circuit equations .... \r", color = :red)
        start = time_ns()

        ## Circuit Equations
        K11 = A .+ 1im.*C
        K12 = -f
        K21 = (1im * omega * z_length) .* transpose(f)
        K22 = Diagonal(ext_resistance) .+ (1im * omega) .* Diagonal(ext_inductance)
        K = [K11 K12; K21 K22]

        ## RHS vector
        T = [zeros(Float64, mesh.nnodes); coil_voltage]

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Circuit equations modelled ("*string(elapsed)*" seconds)                               \n", color = :green)

        ## Factorize the matrix for efficiency
        K_factored = factorize(K)

        ## Compute the numerical solution
        u = K_factored \ T

        ## Time
        printstyled(" ▸ Computing nonlinear solution .... \r", color = :red)
        start = time_ns()

        ## Preallocate buffers for reuse
        u_prev = copy(u)

        ## Save a copy of the initial solution
        u_hist = copy(u)

        for loop = 1:max_itr
            print(" ▸ Computing nonlinear solution (Iteration "*string(Int(loop))*") .... \r")
            
            ## Save a copy of the previous solution
            copy!(u_prev, u)

            ## Update history with a blend of previous and current solution
            @. u_hist = alpha * u_hist + (1 - alpha) * u
            
            ## Compute B norm
            B = fem_future_distribution_grids.Post_Process_Frequency.B_norm(mesh, u_hist)    

            ## Assign new value of mur in core elements
            @. reluctivityperelement[core_elements] = 1.0 / mu_func(abs(B[core_elements]))
            
            ## Update solution
            A = fem_future_distribution_grids.Assemble_Matrices.assemble_A(mesh, reluctivityperelement)
            @. K11 = A + 1im * C
            K[1:mesh.nnodes, 1:mesh.nnodes] = K11

            ## Factorize the matrix for efficiency
            K_factored = factorize(K)
    
            ## Compute the numerical solution
            u = K_factored \ T
            
            ## Check the error against the threshold
            epsilon = 1e-12
            avg_percentage_error = mean(abs.((u .- u_prev) ./ (u_prev .+ epsilon)))
            if avg_percentage_error <= threshold
                break
            end
        end

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Solution computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return u
    end

end