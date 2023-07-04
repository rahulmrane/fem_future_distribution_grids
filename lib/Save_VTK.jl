module Save_VTK

    using WriteVTK

    export save_vtk

    function save_vtk(mesh_data, u, B, H, mag_energy, reluctivityperelement, name)
        mu0 = 4e-7 * pi;
    
        # Define nodes (points) and elements (cells)
        points = [mesh_data.xnode mesh_data.ynode]';
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el) for el in mesh_data.elements];

        # Create VTK file structure using nodes and elements
        vtkfile = vtk_grid("../vtk/"*name*".vtu", points, cells);

        # Store data in the VTK file
        vtkfile["Voltage", VTKPointData()] = u;
        vtkfile["B", VTKCellData()] = B;
        vtkfile["H", VTKCellData()] = H;
        vtkfile["Magnetic Energy", VTKCellData()] = mag_energy;
        vtkfile["Permeability", VTKCellData()]  = 1 ./ (mu0 * reluctivityperelement);

        # Save the file
        outfiles = vtk_save(vtkfile);
    end

end