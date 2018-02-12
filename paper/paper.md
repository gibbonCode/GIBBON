---
title: 'GIBBON: The Geometry and Image-Based Bioengineering add-On'
tags:
  - Computational mechanics
  - Biomechanics
  - Finite Element Analysis
  - Image-Based modeling
  - Visualization
  - Meshing
authors:
 - name: Kevin M Moerman
   orcid: 0000-0003-3768-4269
   affiliation: 1
affiliations:
 - name: Massachusetts Institute of Technology
   index: 1
date: 9 February 2018
bibliography: paper.bib
---

# Summary
[GIBBON](https://www.gibboncode.org), which loosely stands for Geometry and Image-Based Bioengineering add-ON, is a MATLAB® toolbox providing a single open-source framework for many aspects of computational (bio)mechanics such as: image segmentation, meshing, boundary conditions specification, finite element analysis (FEA), and visualization. A schematic of the core functionality of GIBBON is shown in the figure below.

![A Graphical summary of the GIBBON toolbox](GIBBON_overview.png)

Below is a more detailed discussion of the core functionality where reference is made to implementations in the documentation.

* **Image segmentation**: In patient-, or subject-specific biomechanics, the geometry information is often derived from image data (e.g. Magnetic Resonance Imaging (MRI)). GIBBON offers image filtering and smoothing methods, and has a graphical user interface for 3D image segmentation (`HELP_imx.m`). The segmented image data can be converted to 3D surface models (`DEMO_imx_levelset_surface_compare`) which can be meshed for FEA (`HELP_runTetGen`).
* **Computer Aided Design (CAD) tools**: Sometimes geometry is instead imported or designed. Using GIBBON, geometry can be imported from common mesh based CAD files (such as STL, `HELP_import_STL`). For generating geometries within MATLAB®, GIBBON also provides several CAD-style commands such as polygon rounding (`HELP_filletCurve`), revolution (`HELP_polyRevolve`), extrusion (`HELP_polyExtrude`), and sweeping and lofting (`HELP_polyLoftLinear` and `HELP_sweepLoft`). Simple geometries such as spheres (`HELP_geoSphere`), boxes (`HELP_quadBox`), platonic solids (`HELP_platonic_solid`), and rhombic dodecahedra (`HELP_rhombicDodecahedron`) can also be directly created using GIBBON.
* **Surface meshing tools**: 2D multi-region triangular meshing (e.g. `HELP_regionTriMesh2D` and `HELP_multiRegionTriMeshUneven2D`), resampling meshes geodesically (`DEMO_geodesic_remeshing`), smoothing (`DEMO_surface_smooth_methods`), and surface mesh refinement (e.g. `HELP_subtri`, `HELP_subTriDual` and `HELP_subQuad`), mesh type conversions (e.g. `HELP_tri2quad`, `HELP_quad2tri`), and mesh dual computation (`HELP_patch_dual`). Geometries can also be exported to the STL format e.g. for computer aided manufacture and 3D printing.
* **Volumetric meshing**: Tetrahedral meshing (and constrained Delaunay tessellation) of multi-region domains is enabled through an interface with the [TetGen](http://wias-berlin.de/software/tetgen/) [@Si2015] package (`HELP_runTetGen` and `HELP_constrainedDelaunayTetGen`). Hexahedral meshes for some geometry types can be directly coded (e.g. spheres `HELP_hexMeshSphere`, boxes `HELP_hexMeshBox` and lattices `HELP_element2HexLattice`). For general input surfaces multi-region mixed tetrahedral-hexahedral meshing is also available (e.g. `DEMO_MixedTetHexMeshing`).
* **Lattice structures**: One method to generate surface geometry for lattices is the use of triply-periodic functions (`HELP_triplyPeriodicMinimal`). Functions to convert element descriptions, such as tetrahedral and hexahedral elements, to lattice structures have also been implemented (`HELP_element2lattice` and `HELP_element2HexLattice`). These allow for the creation of 3D boundary conforming lattice structures on arbitrary input geometry. Exporting of hexahedral elements is also supported allowing for FEA on the created lattice structures (`DEMO_FEBio_hexLattice_compression`).
* **Finite element analysis (FEA)**: GIBBON interfaces with the free software [FEBio](http://febio.org/) [@Maas2012] for FEA (source code available on FEBio website). GIBBON can be used as a pre- and post- processor for FEBio as it enables code-based development of meshes, boundary conditions, and input files. FEBio files can be directly exported based on dedicated MATLAB® structures (`HELP_febioStruct2xml`). Furthermore, GIBBON can be used to start and control FEBio simulations. As such, iterative and inverse FEA (e.g. based on MATLAB® optimization routines) is also enabled. All `DEMO_febio_...` files are FEBio demos, e.g. `DEMO_febio_0001_cube_uniaxial` is a simple uniaxial loading example, and `DEMO_FEBio_iFEA_uniaxial_01` is an example of inverse FEA. Other demos cover tension, compression, shear, applied forces, applied pressures, applied displacements, bending, poroelasticity, dynamic and viscoelastic analysis, contact and indentation problems, multi-generational materials for pre-load analysis.
* **Visualization**: GIBBON expands the standard MATLAB® visualization capabilities by adding 3D image and voxel visualization (`HELP_im2patch` and `HELP_sliceViewer`), meshed geometries (`HELP_gpatch` and `HELP_meshView`), finite element models (`HELP_element2patch`), and colormapped vector data (`HELP_quiverVec`), and all visualization methods enable multiple colormaps to be used in each figure or axis window. Furthermore GIBBON offers a custom figure window `cFigure` containing 3D rotation options (`HELP_vcw`) that mimic CAD behavior of 3D scene rendering, and high quality figure exporting options (`HELP_efw`). Advanced graphics animation creation and exporting capabilities through a figure window based GUI are also enabled (`HELP_anim8`).

To date GIBBON has been used for image analysis and visualization [@Moerman2012a], continuum mechanics [@Moerman2016], soft tissue biomechanics [@Takaza2013] [@Cooney2015], subject-specific and inverse FEA of soft biological soft tissue _in-vivo_ [@Moerman2017] [@Sengeh2016], and automated generation of parametric scalable models of the lumbar spine [@Lavecchia2018]. The author's personal research with GIBBON is currently focused on subject-specific computational modeling for the automated generation of 3D printable prosthetic devices with optimized and spatially varying mechanical behavior [@Moerman2016c].

# Acknowledgements
GIBBON has not received dedicated funding but development of GIBBON has taken place as part of several other research projects. As such GIBBON has been partially financially supported by: Science Foundation Ireland (Research Frontiers Grant No. 06/RF/ENMO76) in the period 2006-2010, the Dutch Technology Foundation (project 12398) in the period 2012-2015, the Robert Wood Johnson Foundation (RWJF-ID72293) in the period 2015-2016, and National Institute of Health (R01EB024531-01) presently.

# References
