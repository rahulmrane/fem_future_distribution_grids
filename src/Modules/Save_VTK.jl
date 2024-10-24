# #######################################################################################
# Save_VTK.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module Save_VTK

    using WriteVTK

    function save_vtk(mesh, u, B, H, mag_energy, reluctivityperelement, name, current_dir)
        mu0 = 4e-7 * pi
    
        ## Define nodes (points) and elements (cells)
        points = [mesh.xnode mesh.ynode]'
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el.nodes) for el in mesh.Elements]

        ## Create VTK file structure using nodes and elements
        vtkfile = vtk_grid(current_dir * "/../vtk/"*name*".vtu", points, cells)

        ## Store data in the VTK file
        vtkfile["Vector Potential", VTKPointData()] = u
        vtkfile["B", VTKCellData()] = B
        vtkfile["H", VTKCellData()] = H
        vtkfile["Magnetic Energy", VTKCellData()] = mag_energy
        vtkfile["Permeability", VTKCellData()]  = 1 ./ (mu0 * reluctivityperelement)

        ## Save the file
        vtk_save(vtkfile)
    end

    function save_vtk(mesh, u, B, H, mag_energy, reluctivityperelement, Jel, name, current_dir)
        mu0 = 4e-7 * pi
    
        ## Define nodes (points) and elements (cells)
        points = [mesh.xnode mesh.ynode]'
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el.nodes) for el in mesh.Elements]

        ## Create VTK file structure using nodes and elements
        vtkfile = vtk_grid(current_dir * "/../vtk/"*name*".vtu", points, cells)

        ## Store data in the VTK file
        vtkfile["Vector Potential", VTKPointData()] = u
        vtkfile["B", VTKCellData()] = B
        vtkfile["H", VTKCellData()] = H
        vtkfile["Magnetic Energy", VTKCellData()] = mag_energy
        vtkfile["Permeability", VTKCellData()]  = 1 ./ (mu0 * reluctivityperelement)
        vtkfile["Source Current Density", VTKCellData()] = Jel

        ## Save the file
        vtk_save(vtkfile)
    end

    function save_vtk(mesh, u, B, H, mag_energy, reluctivityperelement, Jel, Pv, name, current_dir)
        mu0 = 4e-7 * pi
    
        ## Define nodes (points) and elements (cells)
        points = [mesh.xnode mesh.ynode]'
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el.nodes) for el in mesh.Elements]

        ## Create VTK file structure using nodes and elements
        vtkfile = vtk_grid(current_dir * "/../vtk/"*name*".vtu", points, cells)

        ## Store data in the VTK file
        vtkfile["Vector Potential", VTKPointData()] = u
        vtkfile["B", VTKCellData()] = B
        vtkfile["H", VTKCellData()] = H
        vtkfile["Magnetic Energy", VTKCellData()] = mag_energy
        vtkfile["Permeability", VTKCellData()]  = 1 ./ (mu0 * reluctivityperelement)
        vtkfile["Source Current Density", VTKCellData()] = Jel
        vtkfile["Pv", VTKCellData()] = Pv

        ## Save the file
        vtk_save(vtkfile)
    end

end