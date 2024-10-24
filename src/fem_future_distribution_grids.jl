# #######################################################################################
# fem_future_distribution_grids.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

module fem_future_distribution_grids

    for file in readdir("src/Modules/")
        include("Modules/" * file)
        sym = Symbol(chop(file, tail=3))
        @eval using .$sym
    end

    include("Precompile.jl")

end # module fem_future_distribution_grids