# Future Distribution Grids

The goal of this project is to enhance the EE4375: Finite Element Modelling for Electrical Energy Applications course. It mostly entails applying finite element method (FEM) knowledge to the modeling of a power transformer's magnetic field, thermal field, and electrical circuit parameters.

The code is mostly written in Julia, a programming language, and is used to construct FEM models to calculate the iron and copper losses in power transformers and, in addition, to analyze the substation's temperature.

## Contents
- `lib`: Contains files for the required user defined functions
- `General FEM`: General code for FEM models
  - `first_order`: Use of first order elements
    - `triangle`: Use of triangular elements
      - Contains of the folder :
        - `img`: Contains images obtained
        - `mesh`: Contains GMSH output files
        - `modelling`: Contains .ipynb files for the required code
      - Different cases covered :
        - E-shaped Ferrite Core Transformer Geometry and FEM Modelling
        - STEDIN Transformer Geometry and FEM Modelling
    - `quadrilateral`: Use of quadrilateral elements
    - `hybrid`: Use of hybrid meshing
      - Contains of the folder :
        - `img`: Contains images obtained
        - `mesh`: Contains GMSH output files
        - `modelling`: Contains .ipynb files for the required code
      - Different cases covered :
        - STEDIN Transformer Hybrid Geometry
  - `second_order`: Use of second order elements
    - `triangle`: Use of triangular elements
      - Contains of the folder :
        - `mesh`: Contains GMSH output files
        - `modelling`: Contains .ipynb files for the required code
      - Different cases covered :
        - E-shaped Ferrite Core Transformer Geometry and FEM Modelling
- `Magnetic Field Simulations`
- `Thermal Field Simulations`

## Useful Links
- Description of the project and the required work allotment : https://github.com/ziolai/finite_element_electrical_engineering/blob/main/project-based-assignment/modeling_distribution_transformer.ipynb
- Main repository of the EE4375 course : https://github.com/ziolai/finite_element_electrical_engineering
- Repository of Gijs Lagerweij : https://github.com/gijswl/ee4375_fem_ta
- Repository of Auke Schaap and Philip Soliman : https://github.com/aukeschaap/am-transformers
- Thesis of Max van Dijk : https://repository.tudelft.nl/islandora/object/uuid%3A15b25b42-e04b-4ff2-a187-773bc170f061?collection=education
