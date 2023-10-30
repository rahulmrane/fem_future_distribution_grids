module FEM_Transient_Tri_1e

    using SparseArrays
    using LinearAlgebra
    include("Assemble_Matrices.jl");
    using .Assemble_Matrices
    include("FastSparse.jl");
    using .FastSparse
    include("Post_Process_Time.jl");
    using .Post_Process_Time

    export fem, fem_third_harmonic, fem_nonlinear, fem_nonlinear_third_harmonic

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, time_steps)
        ## Assemble A, B, and f Matrices
        A, B, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)

        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{Float64,1}}(undef, length(time_steps))

        K = factorize(B + dt * A)

        ## Initial condition
        u[1] = zeros(mesh_data.nnodes)

        print(" ▸ Computing solution .... \r")
        start = time_ns()
        ## Backward Euler 
        for k = 2:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            # Compute the solution at the next time step
            u[k] = vec(K \ (B * u[k-1] + dt .* imag(f .* exp(1im*omega*time_steps[k]))))
        end
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")

        return u
    end

    function fem_third_harmonic(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, time_steps, phi_diff)
        ## Assemble A, B, and f Matrices
        A, B, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
    
        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{Float64,1}}(undef, length(time_steps))

        K = factorize(B + dt * A)
    
        ## Initial condition
        u[1] = zeros(mesh_data.nnodes)
    
        print(" ▸ Computing solution .... \r")
        start = time_ns()
        ## Backward Euler 
        for k = 2:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            # Compute the solution at the next time step
            u[k] = vec(K \ (B * u[k-1] + dt .* (imag(f .* exp(1im*omega*time_steps[k])) + imag(f/3 .* exp(1im*3*omega*time_steps[k] + phi_diff)))))
        end
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")

        return u
    end

    
    function mu_func(B)
        k1 = 3.8;
        k2 = 2.17;
        k3 = 396.2;
        mu0 = 4e-7 * pi;
        v = k1 * exp(k2*B^2) + k3;
        return (1 ./ v)
    end
    
    function fem_nonlinear(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, time_steps)
        ## Assemble A, B, and f Matrices
        A, C, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
    
        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{Float64,1}}(undef, length(time_steps))
    
        ## Initial condition
        u[1] = zeros(mesh_data.nnodes)
    
        core_elements = findall(x->x==2, mesh_data.e_group)
        L = LinearIndices(mesh_data.e_group)
        core_elements = L[core_elements]
    
        ## Threshold value for the error
        threshold = 1e-3
        alpha = 0.9
    
        print(" ▸ Computing solution .... \r")
        start = time_ns()
        ## Backward Euler 
        for k = 2:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            u_hist = u[k-1]
            u_prev = u[k-1]
            u_temp = u[k-1]
            for loop = 1:10000
                u_prev = u_temp;
                u_hist = u_hist * alpha + u_temp * (1-alpha)
                
                ## Compute Bnorm
                B = Bnorm_per_timestep(mesh_data, u_hist)

                ## Assign new value of mur
                reluctivityperelement[core_elements] = 1 ./ mu_func.(B[core_elements]);
            
                ## Assemble A matrix again
                A = assemble_A(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)

                # Compute the solution at the next time step
                K = factorize(C + dt * A)
                u_temp = vec(K \ (C * u_hist + dt .* imag(f .* exp(1im*omega*time_steps[k]))))
            
                ## Check the error with the threshold values
                if norm(u_temp-u_prev) <= threshold
                    break;
                end
            end
            u[k] = u_temp
        end
        
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")

        return u
    end
    
    function fem_nonlinear_third_harmonic(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, time_steps, phi_diff)
        ## Assemble A, B, and f Matrices
        A, C, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
    
        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{Float64,1}}(undef, length(time_steps))
    
        ## Initial condition
        u[1] = zeros(mesh_data.nnodes)
    
        core_elements = findall(x->x==2, mesh_data.e_group)
        L = LinearIndices(mesh_data.e_group)
        core_elements = L[core_elements]
    
        ## Threshold value for the error
        threshold = 1e-3
        alpha = 0.9
    
        print(" ▸ Computing solution .... \r")
        start = time_ns()
        ## Backward Euler 
        for k = 2:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            u_hist = u[k-1]
            u_prev = u[k-1]
            u_temp = u[k-1]
            for loop = 1:10000
                u_prev = u_temp;
                u_hist = u_hist * alpha + u_temp * (1-alpha)
                
                ## Compute Bnorm
                B = Bnorm_per_timestep(mesh_data, u_hist)

                ## Assign new value of mur
                reluctivityperelement[core_elements] = 1 ./ mu_func.(B[core_elements]);
            
                ## Assemble A matrix again
                A = assemble_A(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)

                # Compute the solution at the next time step
                K = factorize(C + dt * A)
                u_temp = vec(K \ (C * u_hist + dt .* (imag(f .* exp(1im*omega*time_steps[k])) + imag(f/3 .* exp(1im*3*omega*time_steps[k] + phi_diff)))))
            
                ## Check the error with the threshold values
                if norm(u_temp-u_prev) <= threshold
                    break;
                end
            end
            u[k] = u_temp
        end
        
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")

        return u
    end

end