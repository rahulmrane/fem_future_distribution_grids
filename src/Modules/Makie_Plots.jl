# #######################################################################################
# Makie_Plots.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module Makie_Plots

    using CairoMakie
    CairoMakie.activate!(type = "svg")
    using Colors
    using ColorSchemes

    # Plot a surface plot using Makie
    function plot_surface_tri_1e(msh, colors, name, current_dir)
        cmap = ColorScheme([RGB(0.788, 0.788, 0.784), RGB(0.592, 0.180, 0.996), RGB(1.0, 0.631, 0.216), RGB(0.671, 1.0, 0.224)])
        f, ax, pl = mesh(msh, color = colors, colormap = cmap, shading = Makie.automatic)
        ax.aspect = DataAspect()
        save(current_dir * "/../img/"*name*".png", f)
        return current_figure()
    end

    # Plot a meshed plot using Makie
    function plot_mesh_tri_1e(msh, colors, name, current_dir)
        cmap = ColorScheme([RGB(0.788, 0.788, 0.784), RGB(0.592, 0.180, 0.996), RGB(1.0, 0.631, 0.216), RGB(0.671, 1.0, 0.224)])
        f, ax, pl = mesh(msh, color = colors, colormap = cmap, shading = Makie.automatic)
        wireframe!(ax, msh, color=(:black), linewidth=0.5, transparency=false)
        ax.aspect = DataAspect()
        save(current_dir * "/../img/"*name*".png", f)
        return current_figure()
    end

end