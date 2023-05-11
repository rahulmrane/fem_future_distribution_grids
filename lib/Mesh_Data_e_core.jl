module Mesh_Data_e_core

    export get_mesh_data_tri_1e

    struct mesh_data
        nnodes  # number of nodes
        xnode   # array of x coordinates
        ynode   # array of y coordinates
        nelements   # number of elements
        e_group     # array containing the physical group number of each element
        elements    #  more conveniently structured connectivity array
    end

    function get_mesh_data_tri_1e(gmsh)
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
        ngroup1 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 401)
        ngroup2 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 402)
        ngroup3 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 403)
        ngroup4 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 404)
        ngroup5 = gmsh.model.mesh.getNodesForPhysicalGroup(2, 405)
    
        e_group = zeros(1,nelements)
        elements = [zeros(Int, 3) for i in 1:nelements];
        for element_id in 1:nelements
            node1_id = element_connectivity[1][3*(element_id-1)+1]
            node2_id = element_connectivity[1][3*(element_id-1)+2]
            node3_id = element_connectivity[1][3*(element_id-1)+3]
            G1 = sum(node1_id.== ngroup1[1])+sum(node2_id.== ngroup1[1])+sum(node3_id.== ngroup1[1]) #Air
            G2 = sum(node1_id.== ngroup2[1])+sum(node2_id.== ngroup2[1])+sum(node3_id.== ngroup2[1]) #upper core
            G3 = sum(node1_id.== ngroup3[1])+sum(node2_id.== ngroup3[1])+sum(node3_id.== ngroup3[1]) #lower core
            G4 = sum(node1_id.== ngroup4[1])+sum(node2_id.== ngroup4[1])+sum(node3_id.== ngroup4[1]) #winding right
            G5 = sum(node1_id.== ngroup5[1])+sum(node2_id.== ngroup5[1])+sum(node3_id.== ngroup5[1]) #winding left
            if G1 == 3
                e_group[element_id] = 1;
            elseif G2 == 3 || G3 == 3
                e_group[element_id] = 2;
            elseif G4 == 3
                e_group[element_id] = 4;
            elseif G5 == 3
                e_group[element_id] = 5;
            end

            # Store connectivity in a convenient format
            elements[element_id] = [node1_id, node2_id, node3_id];
        end
    
        return mesh_data(nnodes, xnode, ynode, nelements, e_group, elements)
    end

end