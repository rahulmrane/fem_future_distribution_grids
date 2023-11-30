module Post_Process_Frequency

    using LinearAlgebra
    using StaticArrays

    export post_process, B_norm, core_loss, source_current_density, winding_loss

    function post_process(mesh_data, u, reluctivityperelement)
        print(" ▸ Computing post processing variables .... \r")
        start = time_ns()
    
        Bx = zeros(ComplexF64, mesh_data.nelements);
        By = zeros(ComplexF64, mesh_data.nelements);

        u = u[1:mesh_data.nnodes]
    
        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Emat = SMatrix{3,3}(mesh_data.Eloc[element_id]);
            c = SVector{3}(u[[nodes[1], nodes[2], nodes[3]]]);
            Bx[element_id] = sum(c .* Emat[2,:]);
            By[element_id] = -sum(c .* Emat[1,:]);
        end
        
        Hx = reluctivityperelement .* Bx;
        Hy = reluctivityperelement .* By;

        B = (Bx.^2 + By.^2).^0.5
        H = (Hx.^2 + Hy.^2).^0.5

        mag_energy = 0.5 .* B .* H
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Post processing variables computed ("*string(elapsed)*" seconds)                               ")
    
        return Bx, By, B, Hx, Hy, H, mag_energy
    end

    function B_norm(mesh_data, u, reluctivityperelement)    
        Bx = zeros(ComplexF64, mesh_data.nelements);
        By = zeros(ComplexF64, mesh_data.nelements);

        u = u[1:mesh_data.nnodes]
    
        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Emat = SMatrix{3,3}(mesh_data.Eloc[element_id]);
            c = SVector{3}(u[[nodes[1], nodes[2], nodes[3]]]);
            Bx[element_id] = sum(c .* Emat[2,:]);
            By[element_id] = -sum(c .* Emat[1,:]);
        end

        B = (Bx.^2 + By.^2).^0.5
    
        return B
    end

    function core_loss(mesh_data, B, z_length) 
        print(" ▸ Computing core loss .... \r")
        start = time_ns()
        Pv = zeros(ComplexF64, mesh_data.nelements);
    
        a = 1.53
        b = 1.6
        kv = 2.25
        K = kv * 50^a
        kl = 0.955 # Lamination fill factor
        L = z_length[1]
        Pcore = 0
    
        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Pv[element_id] = K * B[element_id]^b;
            Pcore = Pcore + norm(Pv[element_id]) * mesh_data.area[element_id] * kl * L;
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Core loss computed ("*string(elapsed)*" seconds)                               ")

        return Pv, Pcore
    end

    function source_current_density(mesh_data, u, sourceperelement)
        print(" ▸ Computing current density .... \r")
        start = time_ns()
        Jel = zeros(ComplexF64, mesh_data.nelements);        
    
        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Jel[element_id] = sum(SVector{size(sourceperelement[1])[1]}(u[mesh_data.nnodes+1:mesh_data.nnodes+size(sourceperelement[1])[1]]) .* SVector{size(sourceperelement[1])[1]}(sourceperelement[element_id]))
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Current density computed ("*string(elapsed)*" seconds)                               ")

        return Jel
    end

    function winding_loss(mesh_data, u, ext_resistance)
        print(" ▸ Computing winding loss .... \r")
        start = time_ns()
        Ip = SVector{Int(size(ext_resistance)[1]/2)}(u[mesh_data.nnodes+1:mesh_data.nnodes+Int(size(ext_resistance)[1]/2)])
        Rp = SVector{Int(size(ext_resistance)[1]/2)}(ext_resistance[1:Int(size(ext_resistance)[1]/2)]/2)
        Pwindingp = sum(norm.(Ip) .^2 .* Rp/2)

        Is = SVector{Int(size(ext_resistance)[1]/2)}(u[mesh_data.nnodes+Int(size(ext_resistance)[1]/2)+1:mesh_data.nnodes+size(ext_resistance)[1]])
        Rs = SVector{Int(size(ext_resistance)[1]/2)}(ext_resistance[Int(size(ext_resistance)[1]/2)+1:size(ext_resistance)[1]])
        Pwindings = sum(norm.(Is) .^2 .* Rs/2)

        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Winding loss computed ("*string(elapsed)*" seconds)                               ")

        return Pwindingp, Pwindings
    end

end