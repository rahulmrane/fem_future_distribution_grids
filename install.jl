# #######################################################################################
# install.jl for fem_future_distribution_grids
# Rahul Rane
# #######################################################################################

using Pkg

print("\n\n\tActivating environment in $(pwd())...\n")
Pkg.activate(@__DIR__)
print("\n\n\tInstantiating environment... (i.e. downloading + precompiling packages)\n");
Pkg.instantiate()
Pkg.precompile()