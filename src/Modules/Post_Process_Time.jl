# #######################################################################################
# Post_Process_Time.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module Post_Process_Time

    using LinearAlgebra
    using StaticArrays
    using fem_future_distribution_grids

    function post_process(mesh, u, reluctivityperelement, time_steps)
        ## Time
        printstyled(" ▸ Computing post processing variables .... \r", color = :red)
        start = time_ns()
    
        ## Initialise vectors
        Bx = Vector{eltype(u)}(undef, length(time_steps))
        By = Vector{eltype(u)}(undef, length(time_steps))
        B = Vector{eltype(u)}(undef, length(time_steps))
        Hx = Vector{eltype(u)}(undef, length(time_steps))
        Hy = Vector{eltype(u)}(undef, length(time_steps))
        H = Vector{eltype(u)}(undef, length(time_steps))
        mag_energy = Vector{eltype(u)}(undef, length(time_steps))

        for k in eachindex(time_steps)
            ## Suppress output
            original_stdout = stdout
            redirect_stdout(devnull)

            ## Compute post processing variables at each time step
            Bx[k], By[k], B[k], Hx[k], Hy[k], H[k], mag_energy[k] = fem_future_distribution_grids.Post_Process_Frequency.post_process(mesh, u[k], reluctivityperelement[k])

            ## Release output
            redirect_stdout(original_stdout)
            
            ## Time 
            progress = round(k / length(time_steps) * 100, digits=1)
            elapsed = round((time_ns() - start) / 1e9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            printstyled(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r", color = :red)
        end
    
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Post processing variables computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Bx, By, B, Hx, Hy, H, mag_energy
    end

    function source_current_density(mesh, u, sourceperelement, time_steps)
        ## Time
        printstyled(" ▸ Computing current density .... \r", color = :red)
        start = time_ns()
    
        ## Initialise vectors
        Jel = Vector{eltype(u)}(undef, length(time_steps))

        for k in eachindex(time_steps)
            ## Suppress output
            original_stdout = stdout
            redirect_stdout(devnull)

            ## Compute post processing variables at each time step
            Jel[k] = fem_future_distribution_grids.Post_Process_Frequency.source_current_density(mesh, u[k], sourceperelement)

            ## Release output
            redirect_stdout(original_stdout)
            
            ## Time 
            progress = round(k / length(time_steps) * 100, digits=1)
            elapsed = round((time_ns() - start) / 1e9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            printstyled(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r", color = :red)
        end
    
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Current density computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Jel
    end

    function core_loss(mesh, B, z_length, time_steps)
        ## Time 
        printstyled(" ▸ Computing core loss .... \r", color = :red)
        start = time_ns()
    
        ## Initialise vectors
        Pv = Vector{eltype(B)}(undef, length(time_steps))
        Pcore = eltype(B)(undef, length(time_steps))

        for k in eachindex(time_steps)
            ## Suppress output
            original_stdout = stdout
            redirect_stdout(devnull)

            ## Compute post processing variables at each time step
            Pv[k], Pcore[k] = fem_future_distribution_grids.Post_Process_Frequency.core_loss(mesh, B[k], z_length)

            ## Release output
            redirect_stdout(original_stdout)
            
            ## Time 
            progress = round(k / length(time_steps) * 100, digits=1)
            elapsed = round((time_ns() - start) / 1e9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            printstyled(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r", color = :red)
        end
    
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Core loss computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Pv, Pcore
    end

    function winding_loss(mesh, u, Rp, Rs, time_steps)
        ## Time
        printstyled(" ▸ Computing winding loss .... \r", color = :red)
        start = time_ns()
    
        ## Initialise vectors
        Pwindingp = eltype(u)(undef, length(time_steps))
        Pwindings = eltype(u)(undef, length(time_steps))

        for k = 1:length(time_steps)
            ## Suppress output
            original_stdout = stdout
            redirect_stdout(devnull)

            ## Compute post processing variables at each time step
            Pwindingp[k], Pwindings[k] = fem_future_distribution_grids.Post_Process_Frequency.winding_loss(mesh, u[k], Rp, Rs)

            ## Release output
            redirect_stdout(original_stdout)
            
            ## Time 
            progress = round(k / length(time_steps) * 100, digits=1)
            elapsed = round((time_ns() - start) / 1e9, digits=2)
            estimated = round(length(time_steps) * elapsed / k, digits=2)
            printstyled(" ▸ " * string(progress) * "% (" * string(elapsed) * " of est. " * string(estimated) *" s)                \r", color = :red)
        end
    
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Winding loss computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Pwindingp, Pwindings
    end

end