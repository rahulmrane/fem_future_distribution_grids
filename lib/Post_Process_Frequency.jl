module Post_Process_Frequency

    using LinearAlgebra

    export post_process

    function post_process(mesh_data, u, reluctivityperelement)
        Bx = zeros(Complex{Float64}, mesh_data.nelements);
        By = zeros(Complex{Float64}, mesh_data.nelements);
    
        xnode = mesh_data.xnode;
        ynode = mesh_data.ynode;
    
        ## Perform a loop over the elements
        for (element_id, nodes) in enumerate(mesh_data.elements)
            Emat = mesh_data.Eloc[element_id];
            c = u[[nodes[1], nodes[2], nodes[3]]];
            Bx[element_id] = sum(c .* Emat[2,:]);
            By[element_id] = -sum(c .* Emat[1,:]);
        end
        
        Hx = reluctivityperelement' .* Bx;
        Hy = reluctivityperelement' .* By;

        B = (Bx.^2 + By.^2).^0.5
        H = (Hx.^2 + Hy.^2).^0.5

        mag_energy = 0.5 .* B .* H
    
        return Bx, By, B, Hx, Hy, H, mag_energy
    end

end