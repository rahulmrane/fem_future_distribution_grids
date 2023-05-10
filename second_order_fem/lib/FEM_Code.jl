module FEM_Code

    using SparseArrays
    using LinearAlgebra

    export fem, post_process

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
        ## initialize global matrix A and global vector f
        A = spzeros(Complex{Float64}, mesh_data.nnodes, mesh_data.nnodes)
        f = zeros(Complex{Float64}, mesh_data.nnodes, 1)
    
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
            floc = area_id/3*sourceperelement[element_id]*[1; 1; 1]
        
            #....compute local matrix contribution Aloc of the current element
            Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]] \ UniformScaling(1.);
            Emat[3,:] .= 0;
            Aloc = area_id*reluctivityperelement[element_id]*(transpose(Emat)*Emat);
            Bloc = 1im * area_id / 3 * conductivityperelement[element_id] * omega * Diagonal(ones(3));
        
            #....and add local contribution to global matrices
            f[nodes]       += floc;
            A[nodes,nodes] += (Aloc + Bloc);
        end

        ## Handle the boundary conditions
        A[bnd_node_ids,:] .= 0;
        A[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
        f[bnd_node_ids] .= 0;

        ## Compute the numerical solution
        u = A\f
    
        return u
    end

    function post_process(mesh_data, u, sourceperelement, reluctivityperelement, conductivityperelement, omega)
        Bx = zeros(Complex{Float64}, mesh_data.nelements);
        By = zeros(Complex{Float64}, mesh_data.nelements);
        Jel = zeros(mesh_data.nelements);
    
        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;
    
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
        
            # Calculate eddy current loss
            sigma = conductivityperelement[element_id];
            Jel[element_id] = norm(sourceperelement[element_id] + omega * sigma * 1/3 * sum(c));
        end
        
        Hx = reluctivityperelement' .* Bx;
        Hy = reluctivityperelement' .* By;

        B = (Bx.^2 + By.^2).^0.5
        H = (Hx.^2 + Hy.^2).^0.5

        mag_energy = 0.5 .* B .* H
    
        return Bx, By, B, Hx, Hy, H, mag_energy, Jel
    end

end