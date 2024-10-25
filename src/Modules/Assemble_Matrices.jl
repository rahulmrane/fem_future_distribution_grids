# #######################################################################################
# Assemble_Matrices.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module Assemble_Matrices
    
    using SparseArrays
    using LinearAlgebra
    using StaticArrays
    using fem_future_distribution_grids

    function assemble_matrices(mesh, sourceperelement, reluctivityperelement, conductivityperelement, omega)
        ## Time
        printstyled(" ▸ Constructing matrices .... \r", color = :red)
        start = time_ns()

        ## Compute global matrix A
        A = assemble_A(mesh, reluctivityperelement)

        ## Compute global matrix B
        B = assemble_B(mesh, conductivityperelement, omega)

        ## Compute global vector f
        f = assemble_f(mesh, sourceperelement)
        
        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Matrices constructed ("*string(elapsed)*" seconds)                               \n", color = :green)
    
        return A, B, f
    end

    function assemble_A(mesh, reluctivityperelement)
        ## Initialize global matrix A
        I = Int64[]
        J = Int64[]
        Avalues = Float64[]

        sizehint!(I, 9*mesh.nelements)
        sizehint!(J, 9*mesh.nelements)
        sizehint!(Avalues, 9*mesh.nelements)

        ## Perform a loop over the elements
        for element_id = 1:mesh.nelements
            element = mesh.Elements[element_id]

            Iloc = SVector{9}((element.nodes[i] for i in axes(element.Emat,1) for j in axes(element.Emat,2))...)
            Jloc = SVector{9}((element.nodes[j] for i in axes(element.Emat,1) for j in axes(element.Emat,2))...)
        
            ## Compute local matrix contribution Aloc of the current element
            Amat = SMatrix{3,3}(element.area .* reluctivityperelement[element_id] .* (transpose(element.Emat) * element.Emat))
            Aloc = vec(Amat)

            ## Add local contribution to A
            for j=1:9
                push!(I, Iloc[j])
                push!(J, Jloc[j])
                push!(Avalues, Aloc[j])
            end
        end

        A = sparse(I, J, Avalues, mesh.nnodes, mesh.nnodes)

        ## Handle the boundary conditions
        A[mesh.bnd_node_ids,:] .= 0
        A[mesh.bnd_node_ids,mesh.bnd_node_ids] .= Diagonal(ones(size(mesh.bnd_node_ids)))
    
        return A
    end

    function assemble_B(mesh, conductivityperelement, omega)
        ## Initialize global matrix B
        I = Int64[]
        J = Int64[]
        Bvalues = Float64[]

        sizehint!(I, 9*mesh.nelements)
        sizehint!(J, 9*mesh.nelements)
        sizehint!(Bvalues, 9*mesh.nelements)

        ## Perform a loop over the elements
        for element_id = 1:mesh.nelements
            element = mesh.Elements[element_id]

            Iloc = SVector{9}((element.nodes[i] for i in axes(element.Emat,1) for j in axes(element.Emat,2))...)
            Jloc = SVector{9}((element.nodes[j] for i in axes(element.Emat,1) for j in axes(element.Emat,2))...)
        
            ## Compute local matrix contribution Bloc of the current element     
            Bmat = (element.area ./ 3 .* conductivityperelement[element_id] .* omega) .* SMatrix{3,3}(LinearAlgebra.I)
            Bloc = vec(Bmat)

            ## Add local contribution to B
            for j=1:9
                push!(I, Iloc[j])
                push!(J, Jloc[j])
                push!(Bvalues, Bloc[j])
            end
        end

        B = sparse(I, J, Bvalues, mesh.nnodes, mesh.nnodes)

        ## Handle the boundary conditions
        B[mesh.bnd_node_ids,:] .= 0
        B[mesh.bnd_node_ids,mesh.bnd_node_ids] .= Diagonal(ones(size(mesh.bnd_node_ids)))
    
        return B
    end

    function assemble_f(mesh, sourceperelement::Vector{T}) where T <: Number
        ## Initialize global vector f
        f = zeros(eltype(sourceperelement), mesh.nnodes)
        ones = SVector(1, 1, 1)

        ## Indices of non-zero elements in sourceperelement
        indices = findall(!iszero, sourceperelement)

        ## Perform a loop over the elements
        for element_id in indices
            element = mesh.Elements[element_id]
            
            ## Assemble f matrix
            floc = (element.area/3 .* sourceperelement[element_id]) .* ones
            
            ## Add local contribution to global matrices
            @views f[element.nodes] += floc
        end

        ## Handle the boundary conditions
        f[mesh.bnd_node_ids] .= 0
    
        return f
    end

    function assemble_f(mesh, sourceperelement::Vector{Vector{T}}) where T <: Number
        ## Convert to matrix
        sourceperelementmat = stack(sourceperelement, dims=1)

        ## Number of coils
        ncoils = size(sourceperelementmat)[2]

        ## Initialize global vector f
        f = zeros(eltype(sourceperelementmat), mesh.nnodes, ncoils)

        for coil_id = 1:ncoils
            ## Compute floc for each coil
            floc = assemble_f(mesh, sourceperelementmat[:,coil_id])

            ## Add local contribution to global matrices
            @views f[:,coil_id] = floc
        end
    
        return f
    end

end