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