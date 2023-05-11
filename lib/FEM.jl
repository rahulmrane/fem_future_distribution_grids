module FEM

    include("FEM_Tri_1e.jl");
    using .FEM_Tri_1e

    export fem, post_process

    function fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids, order, mesh_type)
        if cmp(order, "first")
            if cmp(mesh_type, "tri")
                FEM_Tri_1e.fem(mesh_data, sourceperelement, reluctivityperelement, conductivityperelement, omega, bnd_node_ids)
            end
        end
    end

    function post_process(mesh_data, u, sourceperelement, reluctivityperelement, conductivityperelement, omega, order, mesh_type)
        if isequal(order,1)
            if cmp(mesh_type, "tri")
                FEM_Tri_1e.post_process(mesh_data, u, sourceperelement, reluctivityperelement, conductivityperelement, omega)
            end
        end
    end
    
end