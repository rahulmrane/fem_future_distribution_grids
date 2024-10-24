# #######################################################################################
# Post_Process_Frequency.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module Post_Process_Frequency

    using LinearAlgebra
    using StaticArrays

    function post_process(mesh, u, reluctivityperelement)
        ## Time
        printstyled(" ▸ Computing post processing variables .... \r", color = :red)
        start = time_ns()

        ## Initialise vectors
        Bx = Vector{eltype(u)}(undef, mesh.nelements)
        By = Vector{eltype(u)}(undef, mesh.nelements)
        Hx = Vector{eltype(u)}(undef, mesh.nelements)
        Hy = Vector{eltype(u)}(undef, mesh.nelements)
        B = Vector{eltype(u)}(undef, mesh.nelements)
        H = Vector{eltype(u)}(undef, mesh.nelements)
        mag_energy = Vector{eltype(u)}(undef, mesh.nelements)

        @views u_post = u[1:mesh.nnodes]
    
        ## Perform a loop over the elements
        for element_id = 1:mesh.nelements
            element = mesh.Elements[element_id]
        
            uloc = SVector{3}(u_post[element.nodes])
            
            Bx[element_id] = dot(uloc, element.Emat[2, :])
            By[element_id] = -dot(uloc, element.Emat[1, :])
        end
        
        @. Hx = reluctivityperelement * Bx;
        @. Hy = reluctivityperelement * By;

        @. B = (Bx^2 + By^2)^0.5
        @. H = (Hx^2 + Hy^2)^0.5

        @. mag_energy = 0.5 * B * H

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Post processing variables computed ("*string(elapsed)*" seconds)                               \n", color = :green)
    
        return Bx, By, B, Hx, Hy, H, mag_energy
    end

    function B_norm(mesh, u)   
        ## Initialise vectors 
        Bx = Vector{eltype(u)}(undef, mesh.nelements)
        By = Vector{eltype(u)}(undef, mesh.nelements)
        B = Vector{eltype(u)}(undef, mesh.nelements)

        @views u_post = u[1:mesh.nnodes]
    
        ## Perform a loop over the elements
        for element_id = 1:mesh.nelements
            element = mesh.Elements[element_id]
        
            uloc = SVector{3}(u_post[element.nodes])
            
            Bx[element_id] = dot(uloc, element.Emat[2, :])
            By[element_id] = -dot(uloc, element.Emat[1, :])
        end

        @. B = (Bx^2 + By^2)^0.5
    
        return B
    end

    function source_current_density(mesh, u, sourceperelement)
        ## Time
        printstyled(" ▸ Computing current density .... \r", color = :red)
        start = time_ns()

        ## Initialise vectors
        Jel = zeros(eltype(u), mesh.nelements)
        
        ## number of coils
        ncoils = size(sourceperelement[1])[1]
        
        ## Coil currents
        @views u_coils =  u[mesh.nnodes+1:mesh.nnodes+ncoils]
    
        ## Perform a loop over the elements
        for element_id = 1:mesh.nelements
            ## Compute Jel
            Jel[element_id] = dot(u_coils, sourceperelement[element_id])
        end
    
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Current density computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Jel
    end

    function core_loss(mesh, B, z_length)
        ## Time
        printstyled(" ▸ Computing core loss .... \r", color = :red)
        start = time_ns()

        ## Initialise vectors
        Pv = zeros(eltype(B), mesh.nelements)

        ## Initialize total core loss
        Pcore = 0.
    
        ## Coefficients
        a = 1.53
        b = 1.6
        kv = 2.25
        K = kv * 50^a
        kl = 0.955                      # Lamination fill factor
        L = z_length[1]
    
        ## Perform a loop over the elements
        for element_id = 1:mesh.nelements
            element = mesh.Elements[element_id]

            ## Compute Pv
            @views Pv[element_id] = K * B[element_id]^b

            ## Compute Pcore
            @views Pcore += norm(Pv[element_id]) * element.area * kl * L
        end
    
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Core loss computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Pv, Pcore
    end

    function winding_loss(mesh, u, Rp, Rs)
        ## Time
        printstyled(" ▸ Computing winding loss .... \r", color = :red)
        start = time_ns()

        # Number of windings (assuming primary and secondary have the same number)
        nwindings = div(length(u) - mesh.nnodes, 2)

        # Primary side winding current (Ip)
        @views Ip = u[mesh.nnodes+1:mesh.nnodes+nwindings]

        # Primary winding loss
        Pwindingp = sum(abs2.(Ip) .* (Rp / 2))

        # Secondary side winding current (Is)
        @views Is = u[mesh.nnodes+nwindings+1:mesh.nnodes+(2*nwindings)]

        # Secondary winding loss
        Pwindings = sum(abs2.(Is) .* (Rs / 2))

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Winding loss computed ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Pwindingp, Pwindings
    end

end