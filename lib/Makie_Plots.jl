module Makie_Plots

    using CairoMakie
    CairoMakie.activate!(type = "svg")
    using GeometryBasics
    using Colors
    using ColorSchemes
    export get_cell_mesh_tri_1e, plot_surface_tri_1e, plot_mesh_tri_1e

    # Obtain a GeometryBasics.Mesh object suitable for plotting with Makie
    function get_cell_mesh_tri_1e(cell_val, gmsh)
        node_ids, node_coord, _ = gmsh.model.mesh.getNodes()
        nNode = length(node_ids)

        eType, eTag, eConn = gmsh.model.mesh.getElements(2)
        nEl = length(eTag[1])

        points   = zeros(Point{2, Float64}, nEl * 3)   # Array of vertex coordinates (x,y)
        trif     = zeros(TriangleFace{Int}, nEl)      # Array of triangular faces (n1, n2, n3)
        node_val = zeros(Float64, nEl * 3)

        for e = 1:nEl
            n1idx = 3 * (e - 1) + 1; n1 = eConn[1][n1idx]
            n2idx = 3 * (e - 1) + 2; n2 = eConn[1][n2idx]
            n3idx = 3 * (e - 1) + 3; n3 = eConn[1][n3idx]

            points[n1idx] = Point{2}(node_coord[3*(n1-1) + 1], node_coord[3*(n1-1) + 2])
            points[n2idx] = Point{2}(node_coord[3*(n2-1) + 1], node_coord[3*(n2-1) + 2])
            points[n3idx] = Point{2}(node_coord[3*(n3-1) + 1], node_coord[3*(n3-1) + 2])

            node_val[n1idx] = cell_val[e]
            node_val[n2idx] = cell_val[e]
            node_val[n3idx] = cell_val[e]

            trif[e] = (n1idx, n2idx, n3idx)
        end

        msh = GeometryBasics.Mesh(points, trif)

        return msh, node_val
    end

    # Plot a surface plot using Makie
    function plot_surface_tri_1e(cell_val, name, gmsh)
        msh, colors = get_cell_mesh_tri_1e(cell_val, gmsh)
        cmap = ColorScheme([RGB(0.788, 0.788, 0.784), RGB(0.592, 0.180, 0.996), RGB(1.0, 0.631, 0.216), RGB(0.671, 1.0, 0.224)])
        f, ax, pl = mesh(msh, color = colors, colormap = cmap, shading = false)
        ax.aspect = DataAspect()
        save("../img/"*name*".png", f)
        return current_figure()
    end

    # Plot a meshed plot using Makie
    function plot_mesh_tri_1e(cell_val, name, gmsh)
        msh, colors = get_cell_mesh_tri_1e(cell_val, gmsh)
        cmap = ColorScheme([RGB(0.788, 0.788, 0.784), RGB(0.592, 0.180, 0.996), RGB(1.0, 0.631, 0.216), RGB(0.671, 1.0, 0.224)])
        f, ax, pl = mesh(msh, color = colors, colormap = cmap, shading = false)
        wireframe!(ax, msh, color=(:black), linewidth=0.5, transparency=false)
        ax.aspect = DataAspect()
        save("../img/"*name*".png", f)
        return current_figure()
    end

end