module FEM_Transient_Tri_1e

    using SparseArrays
    using LinearAlgebra
    include("FastSparse.jl");
    using .FastSparse
    include("Post_Process_Time.jl");
    using .Post_Process_Time

    export fem, fem_third_harmonic, fem_nonlinear

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, time_steps)
        ## initialize global matrix A & B and global vector f
        f = zeros(Complex{Float64}, mesh_data.nnodes, 1)
        Asp = FastSparseMatrix(mesh_data.nelements)
        Bsp = FastSparseMatrix(mesh_data.nelements)

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            #....compute local vector contribution floc of the current element
            floc = mesh_data.area[element_id]/3*sourceperelement[element_id]*[1; 1; 1]

            #....compute local matrix contribution Aloc of the current element     
            Aloc = mesh_data.area[element_id]*reluctivityperelement[element_id]*(transpose(mesh_data.Eloc[element_id])*mesh_data.Eloc[element_id]);
            Bloc = mesh_data.area[element_id] / 3 * conductivityperelement[element_id] * omega * Matrix{Float64}(I, 3, 3);

            # Add local contribution to A & B
            add!(Asp, element_id, nodes, Aloc);
            add!(Bsp, element_id, nodes, Bloc);

            #....and add local contribution to global matrices
            f[nodes] += floc;
        end

        A = sparse(Asp.i_row, Asp.i_col, Asp.value, mesh_data.nnodes, mesh_data.nnodes, +);
        B = sparse(Bsp.i_row, Bsp.i_col, Bsp.value, mesh_data.nnodes, mesh_data.nnodes, +);

        ## Handle the boundary conditions
        A[bnd_node_ids,:] .= 0;
        A[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
        f[bnd_node_ids] .= 0;

        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{ComplexF64,1}}(undef, length(time_steps))

        K = factorize(B + dt * A)

        ## Initial condition
        u[1] = zeros(mesh_data.nnodes)

        ## Backward Euler 
        for k = 2:length(time_steps)
            # Compute the solution at the next time step
            u[k] = vec((K) \ (B * u[k-1] + dt .* f .* exp(1im*omega*time_steps[k])))
        end

        return u
    end

    function fem_third_harmonic(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, time_steps, phi_diff)
        ## initialize global matrix A & B and global vector f
        f = zeros(Complex{Float64}, mesh_data.nnodes, 1)
        Asp = FastSparseMatrix(mesh_data.nelements)
        Bsp = FastSparseMatrix(mesh_data.nelements)

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            #....compute local vector contribution floc of the current element
            floc = mesh_data.area[element_id]/3*sourceperelement[element_id]*[1; 1; 1]

            #....compute local matrix contribution Aloc of the current element     
            Aloc = mesh_data.area[element_id]*reluctivityperelement[element_id]*(transpose(mesh_data.Eloc[element_id])*mesh_data.Eloc[element_id]);
            Bloc = mesh_data.area[element_id] / 3 * conductivityperelement[element_id] * omega * Matrix{Float64}(I, 3, 3);

            # Add local contribution to A & B
            add!(Asp, element_id, nodes, Aloc);
            add!(Bsp, element_id, nodes, Bloc);

            #....and add local contribution to global matrices
            f[nodes] += floc;
        end

        A = sparse(Asp.i_row, Asp.i_col, Asp.value, mesh_data.nnodes, mesh_data.nnodes, +);
        B = sparse(Bsp.i_row, Bsp.i_col, Bsp.value, mesh_data.nnodes, mesh_data.nnodes, +);

        ## Handle the boundary conditions
        A[bnd_node_ids,:] .= 0;
        A[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
        f[bnd_node_ids] .= 0;
    
        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{ComplexF64,1}}(undef, length(time_steps))

        K = factorize(B + dt * A)
    
        ## Initial condition
        u[1] = zeros(mesh_data.nnodes)
    
        ## Backward Euler 
        for k = 2:length(time_steps)
            # Compute the solution at the next time step
            u[k] = vec((K) \ (B * u[k-1] + dt .* (f .* exp(1im*omega*time_steps[k]) + f/3 .* exp(1im*3*omega*time_steps[k] + phi_diff))))
        end

        return u
    end

    
    function mu_func(B)
        k1 = 3.8;
        k2 = 2.17;
        k3 = 396.2;
        mu0 = 4e-7 * pi;
        v = k1 * exp(k2*B^2) + k3;
        if (1 ./ v)/mu0 < 10
            v = 1 / 10 / mu0;  
        end
        return (1 ./ v)
    end
    
    function fem_nonlinear(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, time_steps)
        ## initialize global matrix A & C and global vector f
        f = zeros(Complex{Float64}, mesh_data.nnodes, 1)
        Csp = FastSparseMatrix(mesh_data.nelements)

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            #....compute local vector contribution floc of the current element
            floc = mesh_data.area[element_id]/3*sourceperelement[element_id]*[1; 1; 1]

            #....compute local matrix contribution Bloc of the current element
            Cloc = mesh_data.area[element_id] / 3 * conductivityperelement[element_id] * omega * Matrix{Float64}(I, 3, 3);

            # Add local contribution to B
            add!(Csp, element_id, nodes, Cloc);

            #....and add local contribution to global matrices
            f[nodes] += floc;
        end
    
        C = sparse(Csp.i_row, Csp.i_col, Csp.value, mesh_data.nnodes, mesh_data.nnodes, +);
    
        ## Handle the boundary conditions
        f[bnd_node_ids] .= 0;
    
        ## Specify time start, end and step
        dt = time_steps[2] - time_steps[1]
        u = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
    
        ## Initial condition
        u[1] = zeros(mesh_data.nnodes)
    
        mur_pts = findall(x->x==2, mesh_data.e_group)
        L = LinearIndices(mesh_data.e_group)
        mur_pts = L[mur_pts]
    
        ## Threshold value for the error
        threshold = 1e-4 .* ones(length(mur_pts))
    
        ## Backward Euler 
        for k = 2:length(time_steps)
            for loop = 1:100
                Asp = FastSparseMatrix(mesh_data.nelements)
            
                ## Perform a loop over the elements
                for (element_id, nodes) in enumerate(mesh_data.elements)
                    #....compute local matrix contribution Aloc of the current element
                    Aloc = mesh_data.area[element_id]*reluctivityperelement[element_id]*(transpose(mesh_data.Eloc[element_id])*mesh_data.Eloc[element_id]);

                    # Add local contribution to A
                    add!(Asp, element_id, nodes, Aloc);
                end

                A = sparse(Asp.i_row, Asp.i_col, Asp.value, mesh_data.nnodes, mesh_data.nnodes, +);

                ## Handle the boundary conditions
                A[bnd_node_ids,:] .= 0;
                A[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))

                K = factorize(C + dt * A)

                # Compute the solution at the next time step
                u[k] = vec((K) \ (C * u[k-1] + dt .* f .* exp(1im*omega*time_steps[k])))
                
                Bx, By, B, Hx, Hy, H, mag_energy = post_process_per_timestep(mesh_data, u[k], reluctivityperelement)
            
                ## Check the error with the threshold values
                if abs.((1 ./ mu_func.(abs.(B[mur_pts]))) - reluctivityperelement[mur_pts]) <= threshold
                    break;
                end

                ## Assign new value of mur
                reluctivityperelement[mur_pts] = 1 ./ mu_func.(abs.(B[mur_pts]));
            end
        end

        return u
    end

end