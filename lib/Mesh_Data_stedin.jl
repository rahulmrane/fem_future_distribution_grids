module Mesh_Data_stedin

    using LinearAlgebra
    using StaticArrays

    export get_mesh_data_tri_1e

    struct meshdata
        nnodes::Int64                              # number of nodes
        nelements::Int64                           # number of elements
        xnode::Vector{Float64}                     # array of x coordinates
        ynode::Vector{Float64}                     # array of y coordinates
        elements::Vector{Vector{Int64}}            # more conveniently structured connectivity array
        e_group::Vector{Int64}                     # array containing the physical group number of each element
        area::Vector{Float64}                      # area of each element
        Eloc::Vector{MMatrix{3, 3, Float64, 9}}    # Emat of each element
    end

    function get_mesh_data_tri_1e(gmsh::Module)
        ## Get and sort the mesh nodes
        node_ids, node_coord, _ = gmsh.model.mesh.getNodes()
        nnodes = length(node_ids)
        
        #..sort the node coordinates by ID, such that Node one sits at row 1
        tosort = [node_ids node_coord[1:3:end] node_coord[2:3:end]];
        sorted = sortslices(tosort , dims = 1);
        node_ids = sorted[:,1]
        xnode = sorted[:,2]
        ynode = sorted[:,3]
    
        ## Get the mesh elements
        element_types, element_ids, element_connectivity = gmsh.model.mesh.getElements(2)
        nelements = length(element_ids[1])
    
        ## Create groups of elements for the subdomains
        ngroup1 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 1)
        ngroup2 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 2)
        ngroup3 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 3)
        ngroup4 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 4)
        ngroup5 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 5)
        ngroup6 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 6)
        ngroup7 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 7)
        ngroup8 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 8)
        ngroup9 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 9)
        ngroup10 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 10)
        ngroup11 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 11)
        ngroup12 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 12)
        ngroup13 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 13)
        ngroup14 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 14)
        ngroup15 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 15)
        ngroup16 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 16)
    
        e_group = zeros(nelements)
        area = zeros(nelements)
        Eloc = []
        elements = [zeros(Int, 3) for i in 1:nelements];
        for element_id in 1:nelements
            node1_id = element_connectivity[1][3*(element_id-1)+1]
            node2_id = element_connectivity[1][3*(element_id-1)+2]
            node3_id = element_connectivity[1][3*(element_id-1)+3]
            
            # Determine which physical group the element belongs to
            G1  = sum(node1_id.== ngroup1[1])+sum(node2_id.== ngroup1[1])+sum(node3_id.== ngroup1[1]) # Oil
            G2  = sum(node1_id.== ngroup2[1])+sum(node2_id.== ngroup2[1])+sum(node3_id.== ngroup2[1]) # Core
            G3  = sum(node1_id.== ngroup3[1])+sum(node2_id.== ngroup3[1])+sum(node3_id.== ngroup3[1]) # HV winding phase 1 left
            G4  = sum(node1_id.== ngroup4[1])+sum(node2_id.== ngroup4[1])+sum(node3_id.== ngroup4[1]) # HV winding phase 1 right
            G5  = sum(node1_id.== ngroup5[1])+sum(node2_id.== ngroup5[1])+sum(node3_id.== ngroup5[1]) # HV winding phase 2 left
            G6  = sum(node1_id.== ngroup6[1])+sum(node2_id.== ngroup6[1])+sum(node3_id.== ngroup6[1]) # HV winding phase 2 right
            G7  = sum(node1_id.== ngroup7[1])+sum(node2_id.== ngroup7[1])+sum(node3_id.== ngroup7[1]) # HV winding phase 3 left
            G8  = sum(node1_id.== ngroup8[1])+sum(node2_id.== ngroup8[1])+sum(node3_id.== ngroup8[1]) # HV winding phase 3 right
            G9  = sum(node1_id.== ngroup9[1])+sum(node2_id.== ngroup9[1])+sum(node3_id.== ngroup9[1]) # LV winding phase 1 left
            G10 = sum(node1_id.== ngroup10[1])+sum(node2_id.== ngroup10[1])+sum(node3_id.== ngroup10[1]) # LV winding phase 1 right
            G11 = sum(node1_id.== ngroup11[1])+sum(node2_id.== ngroup11[1])+sum(node3_id.== ngroup11[1]) # LV winding phase 2 left
            G12 = sum(node1_id.== ngroup12[1])+sum(node2_id.== ngroup12[1])+sum(node3_id.== ngroup12[1]) # LV winding phase 2 right
            G13 = sum(node1_id.== ngroup13[1])+sum(node2_id.== ngroup13[1])+sum(node3_id.== ngroup13[1]) # LV winding phase 3 left
            G14 = sum(node1_id.== ngroup14[1])+sum(node2_id.== ngroup14[1])+sum(node3_id.== ngroup14[1]) # LV winding phase 3 right
            
            if G1 == 3
                e_group[element_id] = 1;
            elseif G2 == 3
                e_group[element_id] = 2;
            elseif G3 == 3
                e_group[element_id] = 3;
            elseif G4 == 3
                e_group[element_id] = 4;
            elseif G5 == 3
                e_group[element_id] = 5;
            elseif G6 == 3
                e_group[element_id] = 6;
            elseif G7 == 3
                e_group[element_id] = 7;
            elseif G8 == 3
                e_group[element_id] = 8;
            elseif G9 == 3
                e_group[element_id] = 9;
            elseif G10 == 3
                e_group[element_id] = 10;
            elseif G11 == 3
                e_group[element_id] = 11;
            elseif G12 == 3
                e_group[element_id] = 12;
            elseif G13 == 3
                e_group[element_id] = 13;
            elseif G14 == 3
                e_group[element_id] = 14;
            end
    
            # Store connectivity in a convenient format
            elements[element_id] = [node1_id, node2_id, node3_id];
            
            #....retrieve the x and y coordinates of the local nodes of the current element
            xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
            ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];
    
            #....compute surface area of the current element
            area[element_id] = ((xnode2 - xnode1)*(ynode3-ynode1) - (xnode3-xnode1)*(ynode2 - ynode1))/2; 
            
            #....compute local matrix contribution Aloc of the current element
            Xmat = SMatrix{3,3}(xnode1, xnode2, xnode3, ynode1, ynode2, ynode3, 1, 1, 1) 
            rhs  = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 1.) 
            Emat = MMatrix{3,3}(Xmat\rhs)
            Emat[3,:] .= 0;
            push!(Eloc,Emat)
        end
    
        return meshdata(nnodes, nelements, xnode, ynode, elements, e_group, area, Eloc)
    end

end