module FEM_Tri_1e

    using SparseArrays
    using LinearAlgebra
    include("FastSparse.jl");
    using .FastSparse

    export fem

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
        ## initialize global matrix A & B and global vector f
        f = zeros(Complex{Float64}, mesh_data.nnodes, 1)
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
            floc = area_id/3*sourceperelement[element_id]*[1; 1; 1]

            #....compute local matrix contribution Aloc of the current element
            Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]] \ UniformScaling(1.);
            Emat[3,:] .= 0;
            Aloc = area_id*reluctivityperelement[element_id]*(transpose(Emat)*Emat);
            Bloc = area_id / 3 * conductivityperelement[element_id] * omega * Matrix{Float64}(I, 3, 3);

            # Add local contribution to A & B
            add!(Asp, element_id, nodes, Aloc);
            add!(Bsp, element_id, nodes, Bloc);

            #....and add local contribution to global matrices
            f[nodes] += floc;
        end

        A = sparse(Asp.i_row, Asp.i_col, Asp.value, mesh_data.nnodes, mesh_data.nnodes, +);
        B = sparse(Bsp.i_row, Bsp.i_col, Bsp.value, mesh_data.nnodes, mesh_data.nnodes, +);

        A = A + 1im*B;

        ## Handle the boundary conditions
        A[bnd_node_ids,:] .= 0;
        A[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
        f[bnd_node_ids] .= 0;

        ## Compute the numerical solution
        u = A\f

        return u
    end

end