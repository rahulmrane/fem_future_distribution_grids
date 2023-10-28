module Assemble_Matrices
    
    using SparseArrays
    using LinearAlgebra
    include("FastSparse.jl");
    using .FastSparse

    export assemble_matrices, assemble_A

    function assemble_matrices(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
        print(" ▸ Constructing Matrices .... \r")
        start = time_ns()
    
        ## initialize global matrix A & B and global vector f
        if string(typeof(sourceperelement)) == "Matrix{Float64}"
            f = zeros(Float64, mesh_data.nnodes, 1)
        elseif string(typeof(sourceperelement)) == "Matrix{ComplexF64}"
            f = zeros(complex(Float64), mesh_data.nnodes, 1)
        else
            f = zeros(Float64, mesh_data.nnodes, size(sourceperelement[1])[2])
        end
        Asp = FastSparseMatrix(mesh_data.nelements)
        Bsp = FastSparseMatrix(mesh_data.nelements)

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            #....compute local vector contribution floc of the current element
            floc = mesh_data.area[element_id]/3*[1; 1; 1]*sourceperelement[element_id]

            #....compute local matrix contribution Aloc of the current element     
            Aloc = mesh_data.area[element_id]*reluctivityperelement[element_id]*(transpose(mesh_data.Eloc[element_id])*mesh_data.Eloc[element_id]);
            Bloc = mesh_data.area[element_id] / 3 * conductivityperelement[element_id] * omega * Matrix{Float64}(I, 3, 3);

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
        B[bnd_node_ids,:] .= 0;
        B[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
        f[bnd_node_ids] .= 0;
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Matrices constructed ("*string(elapsed)*" seconds)                               ")
    
        return A, B, f
    end

    function assemble_A(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
        ## initialize global matrix A
        Asp = FastSparseMatrix(mesh_data.nelements)

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            #....compute local vector contribution floc of the current element
            floc = mesh_data.area[element_id]/3*[1; 1; 1]*sourceperelement[element_id]

            #....compute local matrix contribution Aloc of the current element     
            Aloc = mesh_data.area[element_id]*reluctivityperelement[element_id]*(transpose(mesh_data.Eloc[element_id])*mesh_data.Eloc[element_id]);

            # Add local contribution to A
            add!(Asp, element_id, nodes, Aloc);
        end

        A = sparse(Asp.i_row, Asp.i_col, Asp.value, mesh_data.nnodes, mesh_data.nnodes, +);

        ## Handle the boundary conditions
        A[bnd_node_ids,:] .= 0;
        A[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
    
        return A
    end

end