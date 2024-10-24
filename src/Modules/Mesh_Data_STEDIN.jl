# #######################################################################################
# Mesh_Data_STEDIN.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module Mesh_Data_STEDIN

    try
        using Gmsh: gmsh
    catch
        using gmsh
    end
    using LinearAlgebra
    using StaticArrays
    using GeometryBasics

    ## Struct to hold 2D point
    struct Point
        x::Float64                                                        # x coordinates
        y::Float64                                                        # y coordinates 
    end

    ## Struct to hold a single mesh elemen
    struct Element
        area::Float64                                                     # area of the element 
        Emat::StaticArraysCore.MMatrix{3, 3, Float64, 9}                  # Emat of each element
        nodes::StaticArraysCore.SVector{3, UInt64}                        # list of nodes
    end

    ## Struct to hold entire mesh
    struct Mesh
        nnodes::Int64                                                     # number of nodes
        nelements::Int64                                                  # number of elements
        Elements::Vector{Element}                                         # list of elements
        e_group::Vector{Int64}                                            # array containing the physical group number of each element
        xnode::Vector{Float64}                                            # array of x node values
        ynode::Vector{Float64}                                            # array of y node values
        bnd_node_ids::Vector{UInt64}                                      # indices of nodes where Dirichlet bc are applied
    end

    function get_mesh_data_tri_1e()
        gmsh.initialize()
        
        ## Time
        printstyled(" ▸ Reading mesh file .... \r", color = :red)
        start = time_ns()

        ## Suppress output
        original_stdout = stdout
        redirect_stdout(devnull)

        ## Read mesh from file
        gmsh.open("Geometry/mesh/stedin_transformer.msh")

        ## Release output
        redirect_stdout(original_stdout)

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Mesh file loaded ("*string(elapsed)*" seconds)                               \n", color = :green)

        ## Time
        printstyled(" ▸ Generating required mesh data .... \r", color = :red)
        start = time_ns()
    
        mesh = gmsh.model.mesh::Module
    
        ## Get and sort the mesh nodes
        node_ids, node_coord, _ = mesh.getNodes()::Tuple{Vector{UInt64}, Vector{Float64}, Vector{Float64}}
        nnodes = length(node_ids)
        
        ## Sort the node coordinates by ID, such that Node one sits at row 1
        tosort = Matrix{Float64}(undef, length(node_ids), 3)
        tosort[:, 1] = node_ids
        tosort[:, 2] = @view node_coord[1:3:end]
        tosort[:, 3] = @view node_coord[2:3:end]

        sorted = sortslices(tosort, dims = 1)
        xnode = sorted[:,2]
        ynode = sorted[:,3]
    
        ## Get the mesh elements
        _, element_ids, element_connectivity = mesh.getElements(2)::Tuple{Vector{Int32}, Vector{Vector{UInt64}}, Vector{Vector{UInt64}}}
        nelements = length(element_ids[1])

        ngroup = groups(mesh::Module)
    
        e_group = Vector{Int64}(undef, nelements)
        Elements = Vector{Element}(undef, nelements)
    
        for element_id in 1:nelements
            e1 = element_connectivity[1][3*(element_id-1)+1]
            e2 = element_connectivity[1][3*(element_id-1)+2]
            e3 = element_connectivity[1][3*(element_id-1)+3]
        
            p1 = Point(sorted[e1,2], sorted[e1,3])
            p2 = Point(sorted[e2,2], sorted[e2,3])
            p3 = Point(sorted[e3,2], sorted[e3,3])
        
            area = area_triangle(p1, p2, p3)
        
            nodes = SVector(e1, e2, e3)
            
            ## Determine which physical group the element belongs to
            ## 1  = Oil
            ## 2  = Core
            ## 3  = HV winding phase 1 left
            ## 4  = HV winding phase 1 right
            ## 5  = HV winding phase 2 left
            ## 6  = HV winding phase 2 right
            ## 7  = HV winding phase 3 left
            ## 8  = HV winding phase 3 right
            ## 9  = LV winding phase 1 left
            ## 10 = LV winding phase 1 right
            ## 11 = LV winding phase 2 left
            ## 12 = LV winding phase 2 right
            ## 13 = LV winding phase 3 left
            ## 14 = LV winding phase 3 right
            for i in eachindex(ngroup)
                G = 0
                for j in (e1, e2, e3)
                    G += in(j, ngroup[i][1])
                    if G == 3
                        e_group[element_id] = i
                        break
                    end
                end
            end
            
            ## Compute local matrix contribution Aloc of the current element
            Xmat = SMatrix{3,3}(p1.x, p2.x, p3.x, p1.y, p2.y, p3.y, 1, 1, 1) 
            rhs = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 1.)
            Emat = MMatrix{3,3}(Xmat\rhs)
            Emat[3,:] .= 0
        
            Elements[element_id] = Element(area, Emat, nodes)
        end
    
        ## Handle the boundary conditions
        bnd_node_ids, _ = gmsh.model.mesh.getNodesForPhysicalGroup(1, 1)::Tuple{Vector{UInt64}, Vector{Float64}}
        
        gmsh.finalize()

        ## Time
        elapsed = round((time_ns() - start)/10^9, digits=2)
        printstyled(" ✓ Mesh data generated ("*string(elapsed)*" seconds)                               \n", color = :green)

        return Mesh(nnodes, nelements, Elements, e_group, xnode, ynode, bnd_node_ids)::Mesh_Data_STEDIN.Mesh
    end

    function groups(mesh::Module)
        ngroup = Vector{Tuple{Vector{UInt64}, Vector{Float64}}}(undef, 14)
    
        ## Create groups of elements for the subdomains
        for i in 1:14
            ngroup[i] = mesh.getNodesForPhysicalGroup(2, i)::Tuple{Vector{UInt64}, Vector{Float64}}
        end

        return ngroup::Vector{Tuple{Vector{UInt64}, Vector{Float64}}}
    end

    ## Compute the area of the triangle with vertices p1, p2 and p3
    function area_triangle(p1::Point, p2::Point, p3::Point)
        x12 = p2.x - p1.x
        x13 = p3.x-p1.x
        y12 = p2.y - p1.y
        y13 = p3.y-p1.y
    
        area_id = x12*y13 - x13*y12
        area_id = abs(area_id)/2.
    
        return area_id::Float64
    end

    # Obtain a GeometryBasics.Mesh object suitable for plotting with Makie
    function get_cell_mesh_tri_1e(cell_val)
        gmsh.initialize()

        ## Suppress output
        original_stdout = stdout
        redirect_stdout(devnull)

        ## Read mesh from file
        gmsh.open("Geometry/mesh/stedin_transformer.msh")

        ## Release output
        redirect_stdout(original_stdout)

        mesh = gmsh.model.mesh::Module

        ## Get and sort the mesh nodes
        _, node_coord, _ = mesh.getNodes()::Tuple{Vector{UInt64}, Vector{Float64}, Vector{Float64}}
        
        ## Get the mesh elements
        _, element_ids, element_connectivity = mesh.getElements(2)::Tuple{Vector{Int32}, Vector{Vector{UInt64}}, Vector{Vector{UInt64}}}
        nelements = length(element_ids[1])

        points = zeros(GeometryBasics.Point{2, Float64}, nelements * 3)            # Array of vertex coordinates (x,y)
        trif = zeros(GeometryBasics.TriangleFace{Int}, nelements)                  # Array of triangular faces (n1, n2, n3)
        node_val = zeros(Float64, nelements * 3)

        for element_id in 1:nelements
            n1idx = 3 * (element_id - 1) + 1
            n2idx = 3 * (element_id - 1) + 2
            n3idx = 3 * (element_id - 1) + 3

            n1 = element_connectivity[1][n1idx]
            n2 = element_connectivity[1][n2idx]
            n3 = element_connectivity[1][n3idx]

            points[n1idx] = GeometryBasics.Point{2}(node_coord[3*(n1-1) + 1], node_coord[3*(n1-1) + 2])
            points[n2idx] = GeometryBasics.Point{2}(node_coord[3*(n2-1) + 1], node_coord[3*(n2-1) + 2])
            points[n3idx] = GeometryBasics.Point{2}(node_coord[3*(n3-1) + 1], node_coord[3*(n3-1) + 2])

            node_val[n1idx] = cell_val[element_id]
            node_val[n2idx] = cell_val[element_id]
            node_val[n3idx] = cell_val[element_id]

            trif[element_id] = (n1idx, n2idx, n3idx)
        end

        msh = GeometryBasics.Mesh(points, trif)

        gmsh.finalize()

        return msh, node_val
    end

end