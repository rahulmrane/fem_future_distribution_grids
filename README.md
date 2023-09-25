# Future Distribution Grids

The goal of this project is to enhance the EE4375: Finite Element Modelling for Electrical Energy Applications course. It mostly entails applying finite element method (FEM) knowledge to the modeling of a power transformer's magnetic field, thermal field, and electrical circuit parameters.

The code is mostly written in Julia, a programming language, and is used to construct FEM models to calculate the iron and copper losses in power transformers and, in addition, to analyze the substation's temperature.

## Contents
- `Navigation Pane`: Use this file to navigate around the repository
- `Documentation`: Contains documentation of the required material
- `lib`: Contains files for the required user defined functions
- `General FEM`: General code for FEM models
  - `first_order`: Use of first order elements
    - `triangle`: Use of triangular elements
      - E-shaped Ferrite Core Transformer Geometry and FEM Modelling
      - STEDIN Transformer Geometry and FEM Modelling
    - `quadrilateral`: Use of quadrilateral elements
    - `hybrid`: Use of hybrid meshing
      - STEDIN Transformer Hybrid Geometry
  - `second_order`: Use of second order elements
    - `triangle`: Use of triangular elements
      - E-shaped Ferrite Core Transformer Geometry and FEM Modelling
- `Magnetic Field Simulations`
  - `Frequency Domain`
    - `Current-fed Analysis : Without Eddy Currents`: FEM Analysis with no conductivity of the core
    - `Current-fed Analysis : With Eddy Currents`: FEM Analysis with presence of conductivity of the core
    - `Voltage-fed Coupled Circuit Analysis`: FEM Analysis for Voltage-fed Couple Circuit Analysis
    - `Voltage-fed Analysis : Non-Linear BH Curve`: FEM Analysis incorporating non-linearity of the core material
  - `Time Domain`
    - `Current-fed Analysis`: FEM Analysis with presence of Current-fed coils
    - `Voltage-fed Coupled Circuit Analysis`: FEM Analysis for Voltage-fed Couple Circuit Analysis
    - `Voltage-fed Analysis : Non-Linear BH Curve`: FEM Analysis incorporating non-linearity of the core material

- General structure of subfolders :
  - `img`: Contains images obtained
  - `mesh`: Contains GMSH output files
  - `modelling`: Contains .ipynb files for the required code
  - `vtk`: Contains .vtu files for Paraview visualization

## Useful Links
- Description of the project and the required work allotment : https://github.com/ziolai/finite_element_electrical_engineering/blob/main/project-based-assignment/modeling_distribution_transformer/modeling_distribution_transformer.ipynb
- Main repository of the EE4375 course : https://github.com/ziolai/finite_element_electrical_engineering
- Repository of Gijs Lagerweij : https://github.com/gijswl/ee4375_fem_ta
- Repository of Auke Schaap and Philip Soliman : https://github.com/aukeschaap/am-transformers
- Thesis of Max van Dijk : https://repository.tudelft.nl/islandora/object/uuid%3A15b25b42-e04b-4ff2-a187-773bc170f061?collection=education
- Electrical Equipment and Machines: Finite Element Analysis course (MOOC available on NPTEL, IIT Bombay) PDF: https://drive.google.com/file/d/1wiyJuqohQMM8lPlCGI2hVKD1cLqO6cTP/view
- Electrical Equipment and Machines: Finite Element Analysis course (MOOC available on NPTEL, IIT Bombay) Videos: https://www.youtube.com/playlist?list=PLOzRYVm0a65evnes6W1TFYThU8-4hfqwc
- Finite Element Methods interesting github page: http://hplgit.github.io/INF5620/doc/pub/sphinx-fem/index.html

## Refer below for Interactive Navigation in the Repository
### Documentation
<a href="Documentation.pdf">Documentation</a>

### STEDIN Transformer Geometry
- <a href="Geometry/modelling/stedin_transformer_geometry_definition.ipynb">Modelling Code</a>
- <a href="Geometry/img/stedin_transformer_mesh.png">Meshed Image</a>
- <a href="Geometry/img/stedin_transformer.png">Normal Image</a>

### Frequency Domain Simulations

- Current-fed Analysis : Without Eddy Current Effects
  - <a href="Magnetic Field Simulations/Frequency Domain/Current-fed Analysis - Without Eddy Currents/modelling/stedin_transformer_fem_modelling.ipynb">Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Frequency Domain/Current-fed Analysis - Without Eddy Currents/modelling/single_phase_stedin_transformer_fem_modelling.ipynb">Single-Phase STEDIN Transformer Modelling Code</a>
  
- Current-fed Analysis : With Eddy Current Effects
  - <a href="Magnetic Field Simulations/Frequency Domain/Current-fed Analysis - With Eddy Currents/modelling/stedin_transformer_fem_modelling.ipynb">Three-Phase STEDIN Transformer Modelling Code</a>
  
- Voltage-fed Coupled Circuit Analysis
  - <a href="Magnetic Field Simulations/Frequency Domain/Voltage-fed Coupled Circuit Analysis/modelling/stedin_transformer_fem_modelling.ipynb">Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Frequency Domain/Voltage-fed Coupled Circuit Analysis/modelling/single_phase_stedin_transformer_fem_modelling.ipynb">Single-Phase STEDIN Transformer Modelling Code</a>
  
- Voltage-fed Analysis - Non-Linear BH Curve
  - <a href="Magnetic Field Simulations/Frequency Domain/Voltage-fed Analysis - Non-Linear BH Curve/modelling/stedin_transformer_fem_modelling.ipynb">Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Frequency Domain/Voltage-fed Analysis - Non-Linear BH Curve/modelling/single_phase_stedin_transformer_fem_modelling.ipynb">Single-Phase STEDIN Transformer Modelling Code</a>

#### Required Code Files - Frequency Domain
- <a href="lib/Post_Process_Frequency.jl">Post Processing</a>
- <a href="lib/FEM_Tri_1e.jl">FEM Code for Current-fed Analysis</a>
- <a href="lib/FEM_VoltageFed_Tri_1e.jl">FEM Code for Voltage-fed Analysis</a>

### Time Domain Simulations

- Current-fed Analysis
  - <a href="Magnetic Field Simulations/Time Domain/Current-fed Analysis/modelling/stedin_transformer_fem_modelling.ipynb">Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Time Domain/Current-fed Analysis/modelling/single_phase_stedin_transformer_fem_modelling.ipynb">Single-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Time Domain/Current-fed Analysis/modelling/harmonic_stedin_transformer_fem_modelling.ipynb">Third Harmonic Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Time Domain/Current-fed Analysis/modelling/harmonic_single_phase_stedin_transformer_fem_modelling.ipynb">Third Harmonic Single-Phase STEDIN Transformer Modelling Code</a>
  
- Voltage-fed Coupled Circuit Analysis
  - <a href="Magnetic Field Simulations/Time Domain/Voltage-fed Coupled Circuit Analysis/modelling/stedin_transformer_fem_modelling.ipynb">Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Time Domain/Voltage-fed Coupled Circuit Analysis/modelling/single_phase_stedin_transformer_fem_modelling.ipynb">Single-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Time Domain/Voltage-fed Coupled Circuit Analysis/modelling/harmonic_stedin_transformer_fem_modelling.ipynb">Third Harmonic Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Time Domain/Voltage-fed Coupled Circuit Analysis/modelling/harmonic_single_phase_stedin_transformer_fem_modelling.ipynb">Third Harmonic Single-Phase STEDIN Transformer Modelling Code</a>
  
- Voltage-fed Analysis - Non-Linear BH Curve
  - <a href="Magnetic Field Simulations/Time Domain/Voltage-fed Analysis - Non-Linear BH Curve/modelling/stedin_transformer_fem_modelling.ipynb">Three-Phase STEDIN Transformer Modelling Code</a>
  - <a href="Magnetic Field Simulations/Time Domain/Voltage-fed Analysis - Non-Linear BH Curve/modelling/single_phase_stedin_transformer_fem_modelling.ipynb">Single-Phase STEDIN Transformer Modelling Code</a>

#### Required Code Files - Time Domain
- <a href="lib/Post_Process_Time.jl">Post Processing</a>
- <a href="lib/FEM_Transient_Tri_1e.jl">FEM Code for Current-fed Analysis</a>
- <a href="lib/FEM_Transient_VoltageFed_Tri_1e.jl">FEM Code for Voltage-fed Analysis</a>

### Other Library Files
- <a href="lib/Save_VTK.jl">Saving VTK Files Module</a>
- <a href="lib/Makie_Plots.jl">Makie Plots Module</a>
- <a href="lib/FastSparse.jl">Fast Sparse Module</a>
- Reading Mesh Data
  - <a href="lib/Mesh_Data_stedin.jl">STEDIN Transformer</a>
  - <a href="lib/Mesh_Data_e_core.jl">E-Core Transformer</a>
