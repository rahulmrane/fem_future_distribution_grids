module Post_Process_Time

    using LinearAlgebra

    export post_process, post_process_per_timestep, Bnorm_per_timestep

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
            Bx[k], By[k], B[k], Hx[k], Hy[k], H[k], mag_energy[k] = post_process_per_timestep(mesh_data, u[k], reluctivityperelement)
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        println(" ✓ Post processing variables computed ("*string(elapsed)*" seconds)                               ")

        return Bx, By, B, Hx, Hy, H, mag_energy
    end

    function post_process_per_timestep(mesh_data, u, reluctivityperelement)
        Bx = zeros(Float64, mesh_data.nelements);
        By = zeros(Float64, mesh_data.nelements);

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        u = u[1:mesh_data.nnodes]

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Emat = mesh_data.Eloc[element_id]
            c = u[[nodes[1], nodes[2], nodes[3]]];
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

    function Bnorm_per_timestep(mesh_data, u)
        Bx = zeros(Float64, mesh_data.nelements);
        By = zeros(Float64, mesh_data.nelements);

        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;

        u = u[1:mesh_data.nnodes]

        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Emat = mesh_data.Eloc[element_id]
            c = u[[nodes[1], nodes[2], nodes[3]]];
            Bx[element_id] = sum(c .* Emat[2,:]);
            By[element_id] = -sum(c .* Emat[1,:]);
        end

        B = (Bx.^2 + By.^2).^0.5

        return vec(B)
    end

end