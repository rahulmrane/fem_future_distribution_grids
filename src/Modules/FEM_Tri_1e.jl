# #######################################################################################
# FEM_Tri_1e.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module FEM_Tri_1e

    using LinearAlgebra
    using Statistics
    using fem_future_distribution_grids

    function fem(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega)
        ## Assemble A, B, and f Matrices
        A, B, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega)
        C = A .+ 1im.*B

        ## Time
        printstyled(" ▸ Computing solution .... \r", color = :red)
        start = time_ns()

        ## Factorize the matrix for efficiency
        C_factored = factorize(C)
        
        ## Compute the numerical solution
        u = C_factored \ f
        
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Solution computed ("*string(elapsed)*" seconds)                               \n", color = :green)
    
        return vec(u)
    end

    function mu_func(B)
        k1 = 3.8
        k2 = 2.17
        k3 = 396.2

        v = k1 * exp(k2 * B^2) + k3
        
        return (1 ./ v)
    end

    function fem_nonlinear(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega, core_elements, max_itr)        
        ## Threshold value for the error
        threshold = 1e-4
        alpha = 0.9

        ## Initial solution
        ## Assemble A, B, and f Matrices
        A, C, f = fem_future_distribution_grids.Assemble_Matrices.assemble_matrices(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega)
        D = A .+ 1im.*C

        ## Factorize the matrix for efficiency
        D_factored = factorize(D)
        
        ## Compute the numerical solution
        u = D_factored \ f

        ## Time
        printstyled(" ▸ Computing nonlinear solution .... \r", color = :red)
        start = time_ns()

        ## Save a copy of the initial solution
        u_hist = copy(u)

        for loop = 1:max_itr
            printstyled(" ▸ Computing nonlinear solution (Iteration "*string(Int(loop))*") .... \r", color = :red)
            
            ## Save a copy of the previous solution
            u_prev = copy(u)

            ## Update history with a blend of previous and current solution
            u_hist .= (alpha .* u_hist) .+ ((1 - alpha) .* u)
            
            ## Compute B norm
            B = fem_future_distribution_grids.Post_Process_Frequency.B_norm(mesh, u_hist)

            ## Assign new value of mur in core elements
            reluctivityperelement[core_elements] .= 1 ./ mu_func.(abs.(B[core_elements]))
            
            ## Update solution
            A = fem_future_distribution_grids.Assemble_Matrices.assemble_A(mesh, reluctivityperelement)
            D = A .+ 1im.*C
    
            ## Factorize the matrix for efficiency
            D_factored = factorize(D)
            
            ## Compute the numerical solution
            u = D_factored \ f
            
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