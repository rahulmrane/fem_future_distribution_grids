# Future Distribution Grids

The goal of this project is to enhance the EE4375: Finite Element Modelling for Electrical Energy Applications course. It mostly entails applying finite element method (FEM) knowledge to the modeling of a power transformer's magnetic field, thermal field, and electrical circuit parameters.

The code is mostly written in Julia, a programming language, and is used to construct FEM models to calculate the iron and copper losses in power transformers and, in addition, to analyze the substation's temperature.

## Contents
- `first_order_fem`: Use of First order FEM code
  - 'Contains of the folder :'
    - `lib`: Contains files for the required user defined functions
    - `mesh`: Contains GMSH output files
    - `modelling`: Contains .ipynb files for the required FEM code
  - 'Different cases covered :'
    - E-shaped Ferrite Core Transformer Geometry and FEM Modelling
    - STEDIN Transformer Geometry and FEM Modelling
- `second_order_fem`: Use of First order FEM code
  - 'Contains of the folder :'
    - `mesh`: Contains GMSH output files
    - `modelling`: Contains .ipynb files for the required FEM code
  - 'Different cases covered :'
    - E-shaped Ferrite Core Transformer Geometry and FEM Modelling
- `lib`: Contains files for the required user defined functions

## Useful Links
- Description of the project and the required work allotment : https://github.com/ziolai/finite_element_electrical_engineering/blob/main/project-based-assignment/modeling_distribution_transformer.ipynb
- Main repository of the EE4375 course : https://github.com/ziolai/finite_element_electrical_engineering
- Repository of Gijs Lagerweij : https://github.com/gijswl/ee4375_fem_ta
- Thesis of Max van Dijk : https://repository.tudelft.nl/islandora/object/uuid%3A15b25b42-e04b-4ff2-a187-773bc170f061?collection=education