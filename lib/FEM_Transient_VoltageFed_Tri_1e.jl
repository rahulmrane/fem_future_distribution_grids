module FEM_Transient_VoltageFed_Tri_1e

    using SparseArrays
    using LinearAlgebra
    include("FastSparse.jl");
    using .FastSparse

    export fem, post_process

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, coil_voltage, ext_resistance, ext_inductance, time_steps)
        ## initialize global matrix A & B and global vector f
        f = zeros(Complex{Float64}, mesh_data.nnodes, size(sourceperelement[1])[2])
        Asp = FastSparseMatrix(mesh_data.nelements)
        Bsp = FastSparseMatrix(mesh_data.nelements)

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            #....retrieve global numbering of the local nodes of the current element
            node1_id = nodes[1]; node2_id = nodes[2]; node3_id = nodes[3];

            #....retrieve the x and y coordinates of the local nodes of the current element
            xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
            ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];

            #....compute surface area of the current element
            x12 = xnode2 - xnode1; x13 = xnode3-xnode1;
            y12 = ynode2 - ynode1; y13 = ynode3-ynode1;
            area_id = x12*y13 - x13*y12; area_id = abs(area_id)/2

            #....compute local vector contribution floc of the current element
            floc = area_id/3*[1; 1; 1]*sourceperelement[element_id]

            #....compute local matrix contribution Aloc of the current element
            Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]] \ UniformScaling(1.);
            Emat[3,:] .= 0;
            Aloc = area_id*reluctivityperelement[element_id]*(transpose(Emat)*Emat);
            Bloc = area_id / 3 * conductivityperelement[element_id] * omega * Matrix{Float64}(I, 3, 3);

            # Add local contribution to A & B
            add!(Asp, element_id, nodes, Aloc);
            add!(Bsp, element_id, nodes, Bloc);

            #....and add local contribution to global matrices
            f[nodes,:] += floc;
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

        ## Circuit Equations
        K = [B+dt.*A -dt.*f; f' (dt.*Diagonal(ext_resistance))+Diagonal(ext_inductance)]
        M = [B zeros(mesh_data.nnodes,2); f' Diagonal(ext_inductance)]

        ## Initial condition
        u[1] = zeros(mesh_data.nnodes + size(coil_voltage)[1])

        ## Backward Euler 
        for k = 2:length(time_steps)
            T = [zeros(mesh_data.nnodes); dt.*coil_voltage.*sin(omega*time_steps[k])]
            u[k] = vec(K \ (M * u[k-1] + T))
        end
    
        return u
    end

    function post_process(mesh_data, u, reluctivityperelement, time_steps)
        Bx = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
        By = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
        B = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
        Hx = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
        Hy = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
        H = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
        mag_energy = Vector{Array{ComplexF64,1}}(undef, length(time_steps))
    
        for k = 1:length(time_steps)
            Bx[k], By[k], B[k], Hx[k], Hy[k], H[k], mag_energy[k] = post_process_per_timestep(mesh_data, u[k], reluctivityperelement)
        end
    
        return Bx, By, B, Hx, Hy, H, mag_energy
    end

    function post_process_per_timestep(mesh_data, u, reluctivityperelement)
        Bx = zeros(Complex{Float64}, mesh_data.nelements);
        By = zeros(Complex{Float64}, mesh_data.nelements);
    
        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;
        
        u = u[1:mesh_data.nnodes]
    
        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            #....retrieve global numbering of the local nodes of the current element
            node1_id = nodes[1]; node2_id = nodes[2]; node3_id = nodes[3];
        
            #....retrieve the x and y coordinates of the local nodes of the current element
            xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
            ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];
        
            #....compute local matrix contribution Aloc of the current element
            Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]]\UniformScaling(1.);
            Emat[3,:] .= 0;
            c = u[[node1_id, node2_id, node3_id]];
            Bx[element_id] = sum(c .* Emat[2,:]);
            By[element_id] = -sum(c .* Emat[1,:]);
        end
        
        Hx = reluctivityperelement' .* Bx;
        Hy = reluctivityperelement' .* By;

        B = (Bx.^2 + By.^2).^0.5
        H = (Hx.^2 + Hy.^2).^0.5

        mag_energy = 0.5 .* B .* H
    
        return vec(Bx), vec(By), vec(B), vec(Hx), vec(Hy), vec(H), vec(mag_energy)
    end

end