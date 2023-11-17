module FEM_VoltageFed_Tri_1e

    using LinearAlgebra
    include("Assemble_Matrices.jl");
    using .Assemble_Matrices

    export fem

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, coil_voltage, ext_resistance, ext_inductance, z_length)
        ## Assemble A, B, and f Matrices
        A, B, f = assemble_matrices(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)

        A = A + 1im*B;

        ## Circuit Equations
        K = [A -f; 1im*omega*f'.*z_length Diagonal(ext_resistance)+1im*omega*Diagonal(ext_inductance)]
        T = [zeros(mesh_data.nnodes); coil_voltage]
        
        print(" ▸ Computing solution .... \r")
        start = time_ns()
        ## Compute the numerical solution
        K = factorize(K)
        u = K\T
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Solution computed ("*string(elapsed)*" seconds)                               ")
    
        return u
    end

end