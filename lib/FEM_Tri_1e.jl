module FEM_Tri_1e

    using LinearAlgebra
    include("Assemble_Matrices.jl");
    using .Assemble_Matrices

    export fem

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
        ## Assemble A, B, and f Matrices
        A, B, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
        A = A + 1im*B;
    
        print(" ▸ Computing solution .... \r")
        start = time_ns()
        ## Compute the numerical solution
        A = factorize(A)
        u = A\f
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")
    
        return u
    end

end