module FEM_Transient_VoltageFed_Tri_1e

    using SparseArrays
    using LinearAlgebra
    include("Assemble_Matrices.jl");
    using .Assemble_Matrices
    include("FastSparse.jl");
    using .FastSparse
    include("Post_Process_Time.jl");
    using .Post_Process_Time

    export fem, fem_nonlinear 

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, num_harmonic=[1], phase_diff=[0])
        ## Assemble A, B, and f Matrices
        A, B, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement[1], conductivityperelement, omega, bnd_node_ids)

        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{Float64,1}}(undef, length(time_steps))

        ## Circuit Equations
        K = [B+dt.*A -dt.*f; f'.*z_length (dt.*Diagonal(ext_resistance))+Diagonal(ext_inductance)]
        M = [B zeros(mesh_data.nnodes,size(coil_voltage)[1]); f'.*z_length Diagonal(ext_inductance)]
    
        K = factorize(K)

        ## Initial condition
        u[1] = zeros(mesh_data.nnodes + size(coil_voltage)[1])

        print(" ▸ Computing solution .... \r")
        start = time_ns()
        ## Backward Euler 
        for k = 2:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            S = zeros(size(coil_voltage)[1])
            for e = 1:size(num_harmonic)[1]
                S = S + imag(coil_voltage/num_harmonic[e] .* exp(num_harmonic[e]*1im*omega*time_steps[k] + phase_diff[e]))
            end
            T = [zeros(mesh_data.nnodes); dt.*S]
            u[k] = vec(K \ (M * u[k-1] + T))
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
        
    function fem_nonlinear(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, coil_voltage, ext_resistance, ext_inductance, z_length, time_steps, num_harmonic=[1], phase_diff=[0])
        ## Assemble A, B, and f Matrices
        A, C, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement[1], conductivityperelement, omega, bnd_node_ids)

        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{Float64,1}}(undef, length(time_steps))
    
        ## Initial condition
        u[1] = zeros(mesh_data.nnodes + size(coil_voltage)[1])

        ## Circuit Equations
        K = [C+dt.*A -dt.*f; f'.*z_length (dt.*Diagonal(ext_resistance))+Diagonal(ext_inductance)]
        M = [C zeros(mesh_data.nnodes,size(coil_voltage)[1]); f'.*z_length Diagonal(ext_inductance)]
    
        core_elements = findall(x->x==2, mesh_data.e_group)
        L = LinearIndices(mesh_data.e_group)
        core_elements = L[core_elements]
    
        ## Threshold value for the error
        threshold = 1e-4
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
            reluctivity = reluctivityperelement[k-1]
            for loop = 1:200
                u_prev = u_temp;
                u_hist[1:mesh_data.nnodes] = u_hist[1:mesh_data.nnodes] * alpha + u_temp[1:mesh_data.nnodes] * (1-alpha)
                
                ## Compute Bnorm
                B = Bnorm_per_timestep(mesh_data, u_hist)

                ## Assign new value of mur
                reluctivity[core_elements] = 1 ./ mu_func.(B[core_elements]);
            
                ## Assemble A matrix again
                A = assemble_A(mesh_data, sourceperelement, reluctivity, conductivityperelement, omega, bnd_node_ids)
                
                ## Circuit Equations
                K[1:mesh_data.nnodes,1:mesh_data.nnodes] = C+dt.*A

                # Compute the solution at the next time step
                S = zeros(size(coil_voltage)[1])
                for e = 1:size(num_harmonic)[1]
                    S = S + imag(coil_voltage/num_harmonic[e] .* exp(num_harmonic[e]*1im*omega*time_steps[k] + phase_diff[e]))
                end
                T = [zeros(mesh_data.nnodes); dt.*S]
                u_temp = vec(K \ (M * u_hist + T))
            
                ## Check the error with the threshold values
                if norm(u_temp[1:mesh_data.nnodes]-u_prev[1:mesh_data.nnodes]) <= threshold
                    break;
                end
            end
            u[k] = u_temp
            reluctivityperelement[k] = reluctivity
        end
                                                
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")

        return u
    end

end