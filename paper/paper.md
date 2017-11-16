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

Your submission should probably be somewhere between 250-1000 words.

## Summary
A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience
A clear statement of need that illustrates the purpose of the software
A list of key references including a link to the software archive
Mentions (if applicable) of any ongoing research projects using the software or recent scholarly publications enabled by it

[GIBBON](https://github.com/gibbonCode/GIBBON) [@Moerman2017] is an open-source MATLAB® toolbox for image-based computational (bio)mechanics. GIBBON allows for image-segmentation, surface meshing, volumetric meshing, and boundary conditions specification and finite element analysis.

The core functionality covers:
* **Image processing and segmentation**
* **Surface geometry creation, meshing, and mesh adjustment**
    * Remeshing, smoothening, refinement.
* **Volumetric meshing**
    * Tetrahedral meshing (and constrained Delaunay tessellation) is enabled through an interface with the [TetGen](http://wias-berlin.de/software/tetgen/) [@Si2015] package.
    * Basic mixed tetrahedral-hexahedral meshing through mesh pillowing
    * Hexahedral meshing of geometry primitives and latticed structures
* **Lattice structure creation**
* **Boundary conditions specification**
* **Finite element analysis (FEA)** is enabled by acting as a wrapper for the [FEBio](http://febio.org/) project [@Maas2012]. GIBBON can be used as a pre- and post- processor for FEA as it enables the code based development of meshes, boundary conditions, and input files. FEBio files can be directly exported based on dedicated MATLAB® structures. Furthermore, GIBBON can be used to start and control FEBio simulations. As such, iterative and inverse FEA (e.g. based on MATLAB® optimization routines) is also enabled.  

![A Graphical summary of the GIBBON toobox](/docs/img/GIBBON_overview.jpg)
_A Graphical summary of the GIBBON toobox_

## References
