# Future Distribution Grids

The goal of this project is to enhance the EE4375: Finite Element Modelling for Electrical Energy Applications course. It mostly entails applying finite element method (FEM) knowledge to the modeling of a power transformer's magnetic field, thermal field, and electrical circuit parameters.

The code is mostly written in Julia, a programming language, and is used to construct FEM models to calculate the iron and copper losses in power transformers and, in addition, to analyze the substation's temperature.

## Contents
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
  - `Without Eddy Currents`: FEM Analysis with no conductivity of the core
  - `With Eddy Currents`: FEM Analysis with presence of conductivity of the core
  - `Non-Linear BH Curve`: FEM Analysis incorporating non-linearity of the core material
  - `Voltage-fed Couple Circuit Analysis`: FEM Analysis for Voltage fed Couple Circuit Analysis
  - `Transient Analysis`: FEM Analysis for Transient Analysis
  - `Trasient + Voltage-fed Couple Circuit Analysis`: FEM Analysis for Transient + Voltage fed Couple Circuit Analysis

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
