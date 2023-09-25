module Post_Process_Time

    using LinearAlgebra

    export post_process

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