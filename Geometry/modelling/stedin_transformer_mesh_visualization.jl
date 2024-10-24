# #######################################################################################
# STEDIN Transformer Mesh Visualization
# Rahul Rane
# #######################################################################################

try
    using Gmsh: gmsh
catch
    using gmsh
end

gmsh.initialize()

gmsh.open("Geometry/mesh/stedin_transformer.msh")

if true 
    gmsh.fltk.run() 
end

gmsh.finalize()

using Plots

A = Vector(0:0.001:0.06)
omega = 2*pi*50
X = imag(1*exp.(1im*omega*A))
Y = imag(1*exp.(1im*omega*A)*exp(1im*2*pi/3))
Z = imag(1*exp.(1im*omega*A)*exp(-1im*2*pi/3))
display(plot(A,[X,Y,Z]))

print(" ✓ Code execution completed !!!!")