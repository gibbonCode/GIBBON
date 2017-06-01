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
[GIBBON](www.gibboncode.org) (The Geometry and Image-Based Bioengineering add-On) is an open-source MATLAB toolbox for segmentation, image-based modeling, visualization, meshing, and finite element analysis. GIBBON includes an array of image and geometry visualization and processing tools and is interfaced with free open source software such as [TetGen](http://wias-berlin.de/software/tetgen/) [@Si2015], for robust tetrahedral meshing, and [FEBio](http://febio.org/) [@Maas2012] for finite element analysis. The combination provides a highly flexible image-based modelling environment and enables advanced inverse finite element analysis.

![Overview of GIBBON](docs/html/GIBBON_overview.jpg){:width="100%"}

# Application highlights
## Segmentation  
The `imx.m` function provides a graphical user interface for segmenting 3D image data.

## Meshing  
`runTetGen.m`

## Geometric tools  
Models, `geoSphere.m`
Mesh conversions `hex2tet.m`
Refinement, resampling `subTri.m`, `triSurfaceRemesh.m`  
Lattice structures  
CAD (computer aided design) tools `polyExtrude.m`, `polyLoftLinear.m`

## Visualization
`cFigure.m` `vcw.m` `im2patch.m` `gpatch.m` `quiver3DPatch.m` `fourthOrderTensorView.m`

## Finite Element Analysis
`runMonitorFEBio.m`

# Setting up GIBBON
## Installation
The GIBBON website (and the `readme.md` file) contains detailed installation instructions. For basic use only the GIBBON folder and sub-directories need to be added to the MATLAB path. For use with FEBio for finite element analysis FEBio needs to be installed and the patch to FEBio needs to be provided. If exporting of figures to images and animations is of interest the export_fig library needs to be downloaded and added to the MATLAB path also.

## Third party packages
* __FEBio__ [@Maas2012], a finite element solver, which enables finite element analysis with GIBBON. FEBio is available through [FEBio](http://febio.org/) website.
* __TetGen__, an open source tetrahedral meshing package [@Si2015]. TetGen is utilized in GIBBON for tetrahedral meshing, see also the `runTetGen.m` function. TetGen is included with GIBBON but alternative versions may be available through the [TetGen](http://wias-berlin.de/software/tetgen/) website.
* __export_fig__, a MATLAB library for exporting publication quality figures. export_fig is used in GIBBON with the `efw.m` (export figure widget) to export figures, and with `anim8.m` (animation function) to export animated .gif files. export_fig can be obtained via the [export_fig](https://github.com/altmany/export_fig) GitHub repository.

## Use of GIBBON in research


# References
