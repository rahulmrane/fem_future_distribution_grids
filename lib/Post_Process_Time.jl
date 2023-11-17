module Post_Process_Time

    using LinearAlgebra
    using StaticArrays

    export post_process, Bnorm_per_timestep, core_loss, source_current_density, winding_loss

    function post_process(mesh_data, u, reluctivityperelement, time_steps)
        print(" ▸ Computing post processing variables .... \r")
        start = time_ns()
    
        Bx = Vector{Array{Float64,1}}(undef, length(time_steps))
        By = Vector{Array{Float64,1}}(undef, length(time_steps))
        B = Vector{Array{Float64,1}}(undef, length(time_steps))
        Hx = Vector{Array{Float64,1}}(undef, length(time_steps))
        Hy = Vector{Array{Float64,1}}(undef, length(time_steps))
        H = Vector{Array{Float64,1}}(undef, length(time_steps))
        mag_energy = Vector{Array{Float64,1}}(undef, length(time_steps))

        for k = 1:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            Bx[k], By[k], B[k], Hx[k], Hy[k], H[k], mag_energy[k] = post_process_per_timestep(mesh_data, u[k], reluctivityperelement[k])
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Post processing variables computed ("*string(elapsed)*" seconds)                               ")

        return Bx, By, B, Hx, Hy, H, mag_energy
    end

    function post_process_per_timestep(mesh_data, u, reluctivityperelement)
        Bx = zeros(Float64, mesh_data.nelements);
        By = zeros(Float64, mesh_data.nelements);

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

        return vec(Bx), vec(By), vec(B), vec(Hx), vec(Hy), vec(H), vec(mag_energy)
    end

    function Bnorm_per_timestep(mesh_data, u)
        Bx = zeros(Float64, mesh_data.nelements);
        By = zeros(Float64, mesh_data.nelements);

        u = u[1:mesh_data.nnodes]

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Emat = SMatrix{3,3}(mesh_data.Eloc[element_id]);
            c = SVector{3}(u[[nodes[1], nodes[2], nodes[3]]]);
            Bx[element_id] = sum(c .* Emat[2,:]);
            By[element_id] = -sum(c .* Emat[1,:]);
        end

        B = (Bx.^2 + By.^2).^0.5

        return vec(B)
    end

    function core_loss(mesh_data, B, z_length, time_steps) 
        print(" ▸ Computing core loss .... \r")
        start = time_ns()
    
        Pv = Vector{Array{Float64,1}}(undef, length(time_steps))
        Pcore = Vector{Float64}(undef, length(time_steps))

        for k = 1:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            Pv[k], Pcore[k] = core_loss_per_timestep(mesh_data, B[k], z_length)
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Core loss computed ("*string(elapsed)*" seconds)                               ")

        return Pv, Pcore
    end

    function core_loss_per_timestep(mesh_data, B, z_length)
        Pv = zeros(Float64, mesh_data.nelements);
    
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

        return vec(Pv), Pcore
    end

    function source_current_density(mesh_data, u, sourceperelement, time_steps)
        print(" ▸ Computing current density .... \r")
        start = time_ns()
    
        Jel = Vector{Array{Float64,1}}(undef, length(time_steps))

        for k = 1:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            Jel[k] = source_current_density_per_timestep(mesh_data, u[k], sourceperelement)
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Current density computed ("*string(elapsed)*" seconds)                               ")

        return Jel
    end

    function source_current_density_per_timestep(mesh_data, u, sourceperelement)
        Jel = zeros(Float64, mesh_data.nelements);        
    
        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Jel[element_id] = sum(SVector{size(sourceperelement[1])[1]}(u[mesh_data.nnodes+1:mesh_data.nnodes+size(sourceperelement[1])[1]]) .* SVector{size(sourceperelement[1])[1]}(sourceperelement[element_id]))
        end

        return vec(Jel)
    end

    function winding_loss(mesh_data, u, ext_resistance, time_steps)
        print(" ▸ Computing winding loss .... \r")
        start = time_ns()
    
        Pwindingp = Vector{Float64}(undef, length(time_steps))
        Pwindings = Vector{Float64}(undef, length(time_steps))

        for k = 1:length(time_steps)
            progress = round(k/length(time_steps)*100, digits=1)
            elapsed = round((time_ns() - start)/10^9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            print(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r")
            Pwindingp[k], Pwindings[k] = winding_loss_per_timestep(mesh_data, u[k], ext_resistance)
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Winding loss computed ("*string(elapsed)*" seconds)                               ")

        return Pwindingp, Pwindings
    end

    function winding_loss_per_timestep(mesh_data, u, ext_resistance)
        Ip = SVector{Int(size(ext_resistance)[1]/2)}(u[mesh_data.nnodes+1:mesh_data.nnodes+Int(size(ext_resistance)[1]/2)])
        Rp = SVector{Int(size(ext_resistance)[1]/2)}(ext_resistance[1:Int(size(ext_resistance)[1]/2)]/2)
        Pwindingp = sum(norm.(Ip) .^2 .* Rp/2)

        Is = SVector{Int(size(ext_resistance)[1]/2)}(u[mesh_data.nnodes+Int(size(ext_resistance)[1]/2)+1:mesh_data.nnodes+size(ext_resistance)[1]])
        Rs = SVector{Int(size(ext_resistance)[1]/2)}(ext_resistance[Int(size(ext_resistance)[1]/2)+1:size(ext_resistance)[1]])
        Pwindings = sum(norm.(Is) .^2 .* Rs/2)

        return Pwindingp, Pwindings
    end

end