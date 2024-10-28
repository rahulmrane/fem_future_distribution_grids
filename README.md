# Future Distribution Grids

The project aims to enhance the EE4375 course, Finite Element Modelling for Electrical Energy Applications, at Delft University of Technology. This involves applying finite element method (FEM) techniques to model a power transformer's magnetic and thermal fields, as well as its electrical circuit parameters.

The code, primarily written in Julia, is used to develop FEM models that calculate iron and copper losses in power transformers and analyze substation temperature distribution.

## Contents
- `Documentation`: Contains documentation of the required material.
- `src`: Contains files for the required user defined functions.
- `Magnetic Field Simulations`
  - `Frequency Domain`
    - `Current-fed Analysis`: FEM analysis with current-fed coils.
    - `Voltage-fed Coupled Circuit Analysis`: FEM analysis of a voltage-fed coupled circuit.
    - `Non-Linear Voltage-fed Analysis`: FEM analysis with non-linear core material.
  - `Time Domain`
    - `Current-fed Analysis`: FEM analysis with current-fed coils.
    - `Non-Linear Current-fed Analysis`: FEM analysis with non-linear core material.
    - `Voltage-fed Coupled Circuit Analysis`: FEM analysis of a voltage-fed coupled circuit.
    - `Non-Linear Voltage-fed Analysis`: FEM analysis with non-linear core material.

- General structure of subfolders :
  - `img`: Contains obtained images.
  - `mesh`: Contains GMSH output files.
  - `modelling`: Contains .jl files with the code.
  - `vtk`: Contains .vtu files for Paraview visualization.

## Useful Links
Here’s a refined list of references for easy access and to help maintain clarity.
- Project Description and Work Allotment : <a href="https://github.com/ziolai/finite_element_electrical_engineering/blob/main/project-based-assignment/modeling_distribution_transformer/modeling_distribution_transformer.ipynb">Modeling Distribution Transformer Project Assignment</a>
  - Contains a detailed description of the project goals, structure, and specific work requirements for modeling a distribution transformer.
- Main Repository for EE4375 Course : <a href="https://github.com/ziolai/finite_element_electrical_engineering">EE4375 Finite Element Electrical Engineering Repository</a>
  - Central repository for the course, containing core resources, assignments, and FEM models related to electrical engineering.
- Gijs Lagerweij's Repository : <a href="https://github.com/gijswl/ee4375_fem_ta">EE4375 FEM Teaching Assistant Repository</a>
  - Repository by TA Gijs Lagerweij, containing additional resources, code examples, and guidance for the EE4375 course.
- Auke Schaap and Philip Soliman's Repository : <a href="https://github.com/aukeschaap/am-transformers">AM Transformers Project Repository</a>
  - Contains project materials by Auke Schaap and Philip Soliman focusing on advanced transformer modeling techniques.
- Thesis by Max van Dijk : <a href="https://repository.tudelft.nl/islandora/object/uuid%3A15b25b42-e04b-4ff2-a187-773bc170f061?collection=education">Finite Element Analysis Thesis, TU Delft Repository</a>
  - Max van Dijk’s research thesis providing in-depth insights into finite element methods applied to electrical engineering problems.
- NPTEL Finite Element Analysis Course (PDF) : <a href="https://drive.google.com/file/d/1wiyJuqohQMM8lPlCGI2hVKD1cLqO6cTP/view">Finite Element Analysis for Electrical Equipment and Machines (IIT Bombay)</a>
  - Downloadable course material from the NPTEL MOOC covering FEM applications in electrical equipment and machines.
- NPTEL Finite Element Analysis Course (Video Series) : <a href="https://www.youtube.com/playlist?list=PLOzRYVm0a65evnes6W1TFYThU8-4hfqwc">YouTube Playlist - Finite Element Analysis for Electrical Equipment and Machines</a>
  - Complete video lecture series of the IIT Bombay NPTEL MOOC on FEM in electrical engineering.
- Finite Element Methods - HPLGit Resources : <a href="http://hplgit.github.io/INF5620/doc/pub/sphinx-fem/index.html">HPLGit Sphinx-FEM Documentation</a>
  - A comprehensive online resource detailing finite element methods, with code examples, tutorials, and theoretical explanations.
