---
title: 'GIBBON: The Geometry and Image-Based Bioengineering add-On'
tags:
  - Image-Based modeling
  - Visualization
  - Meshing
  - Finite Element Analysis
authors:
 - name: Kevin M Moerman
   orcid: 0000-0003-3768-4269
   affiliation: 1
affiliations:
 - name: Massachusetts Institute of Technology
   index: 1
date: 1 June 2017
bibliography: paper.bib
---

# Summary
[GIBBON](https://github.com/gibbonCode/GIBBON) [@Moerman2017a], which loosely stands for Geometry and Image-Based Bioengineering add-ON, is an open-source MATLAB® toolbox providing a single open-source framework for many aspects of computational (bio)mechanics such as: image segmentation, meshing, boundary conditions specification, and finite element analysis (FEA). Visualizations for the core functionality are shown in the figure below.

![A Graphical summary of the GIBBON toobox](GIBBON_overview.png)

A more detailed discussion of the core functionality and example implementations with links to documentation are next discussed.

* **Image segmentation**: In patient-, or subject-specific biomechanics, the geometry information is often derived from image data (e.g. Magnetic Resonace Imaging (MRI)). GIBBON offers image filtering and smoothening methods, and has a graphical user interface for 3D image segmentation (see `HELP_imx.m`). The segmented image data can be converted to 3D surface models which can be meshed for FEA.
* **Computer Aided Design (CAD) tools**: Sometimes geometry is instead imported or designed. Using GIBBON geometry can be imported from common mesh based CAD files (such as STL, see `HELP_import_STL`). However, several CAD style commands have also been implemented in GIBBON such as polygon rounding, revolution (see `HELP_polyRevolve`), extrusion (see `HELP_polyExtrude`), and sweeping and lofting (see `HELP_polyLoftLinear` and `HELP_sweepLoft`).
* **Surface meshing tools**: Two-dimensional multi-region triangular meshing (See functions like `HELP_regionTriMesh2D` and `HELP_multiRegionTriMeshUneven2D`), resampling meshes geodesically (see `DEMO_geodesic_remeshing`), smoothening (`DEMO_surface_smooth_methods`), and surface mesh refinement (see `HELP_subtri`, `HELP_subTriDual` and `HELP_subQuad`), and mesh type conversions (see `HELP_tri2quad`, `HELP_quad2tri`), and mesh dual computation (see `HELP_patch_dual`).
* **Volumetric meshing**: Tetrahedral volumetric meshing meshing (and constrained Delaunay tessellation) is enabled through an interface with the [TetGen](http://wias-berlin.de/software/tetgen/) [@Si2015] package (`HELP_runTetGen`). Hexahedral meshes for some geometry types can be directly coded (e.g. sphere `HELP_hexMeshSphere`, boxes `HELP_hexMeshBox` and lattices `HELP_element2HexLattice`). For general input surfaces mixed tetrahedral-hexahedral meshing is also available (`DEMO_MixedTetHexMeshing`).
* **Lattice structures**: Volumetric and surface conforming lattices
* **Finite element analysis (FEA)**: FEA is enabled by acting as a wrapper for the [FEBio](http://febio.org/) project [@Maas2012]. GIBBON can be used as a pre- and post- processor for FEA as it enables the code based development of meshes, boundary conditions, and input files. FEBio files can be directly exported based on dedicated MATLAB® structures. Furthermore, GIBBON can be used to start and control FEBio simulations. As such, iterative and inverse FEA (e.g. based on MATLAB® optimization routines) is also enabled.
* **Visualization**  
GIBBON expands the standard MATLAB® visualization capabilities by adding 3D image and voxel visualization (see `HELP_im2patch`), meshed geometries (`HELP_XXXX`), finite element models (`HELP_element2patch`), and all visualization methods enable multiple colormaps to be used in each figure or axis window. Furthermore GIBBON offers a custom figure window `cFigure` containing 3D rotation options (`HELP_vcw`) that mimic CAD behavior of 3D scene rendering, and high quality figure exporting options (`HELP_efw`).

To date GIBBON has been used for image analysis and visualization [@Moerman2012a], continuum mechanics[@Moerman2017a], soft tissue biomechanics ([@Cooney2015] [@Takaza2013]), and subject-specific and inverse FEA of soft biological soft tissue ([@Moerman2017a] [@Sengeh2016]). Current research with GIBBON is also focused on subject-specific computational modeling for the automated generation of 3D printable prosthetic devices with optimized and spatially varying mechanical behavior [@Moerman2016c].

# References
