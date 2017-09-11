---
layout: page
title: "About"
logo: "img/home-bg.jpg"
description: "About GIBBON"
header-img: "img/home-bg.jpg"
---

# Summary
[GIBBON](www.gibboncode.org) (The Geometry and Image-Based Bioengineering add-On) is an open-source MATLAB toolbox for segmentation, image-based modeling, visualization, meshing, and finite element analysis. GIBBON includes an array of image and geometry visualization and processing tools and is interfaced with free open source software such as [TetGen](http://wias-berlin.de/software/tetgen/), for robust tetrahedral meshing, and [FEBio](http://febio.org/) for finite element analysis. The combination provides a highly flexible image-based modelling environment and enables advanced inverse finite element analysis.   
![Overview of GIBBON](html/GIBBON_overview.jpg){:width="100%"}

# Application highlights   
- [Segmentation](#Segmentation)  
- [Meshing](#Meshing)  
- [Finite element analysis](#FEA)
- [Geometric tools](#Geometric)
- [Visualization](#Visualization)    

### Segmentation  <a name="Segmentation"></a>    
The `imx.m` function provides a graphical user interface for segmenting 3D image data. The demo below stems from: `HELP_imx`.   
![Segmentation](/img/imx_demo.gif){:width="70%"}

### Meshing <a name="Meshing"></a>   
Multi-material tetrahedral meshing is enabled using [TetGen](http://wias-berlin.de/software/tetgen/). The TetGen interface is based on the `runTetGen.m` function. The demo below comes from the help file `HELP_runTetGen`.
![Tetrahedral meshing](/img/bunnyMesh.gif){:width="70%"}

### Finite Element Analysis <a name="FEA"></a>   
Finite element analysis is enabled through the FEBio interface (see also the `runMonitorFEBio` function.
The image below is for large strain analysis of a twisting bar and stems from the demo `DEMO_FEBio_bar_twist`. Other `DEMO_FEBio_...` files cover uni-axial tension/compression, bending, indentation, viscoelastic analysis, contact and indentation problems, multi-generational materials for pre-load analysis.   
![Large strain analysis](/img/barTwist.gif){:width="70%"}

### Geometric tools <a name="Geometric"></a>   
Models, `geoSphere`   
Mesh conversions `hex2tet`   
Refinement, resampling `subTri`, `triSurfaceRemesh`     
Lattice structures  `element2lattice`   
CAD (computer aided design) tools `polyExtrude`, `polyLoftLinear`, `sweepLoft`

### Visualization <a name="Visualization"></a>    
`im2patch` `gpatch` `quiver3DPatch` `fourthOrderTensorView`

GIBBON uses a custom figure type called `cFigure`. It contains a white background, is maximized by default, and includes a view control widget (`vcw`) and an export figure widget (`efw`). The `vcw` function enables 3D view control similar to CAD packages. Rotation, zooming, and panning can be performed using a 3 button "mouse" (or equivalent input device), e.g. middle-click rotates the view, left click drags/moves the view, right click allows zooming in and out.

| GIBBON | Classic MATLAB |   
|:--|--:|   
|  |  |   

 ![](/img/gibbonViewControl.gif){:width="100%"}

# Setting up GIBBON
See [installation](www.gibboncode.org/Installation/).

## Use of GIBBON in research
