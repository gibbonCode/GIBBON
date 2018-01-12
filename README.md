
<img src="docs/html/gibbonLogo.jpg" href="https://gibboncode.org" alt="GIBBON" width="100%">   

[![DOI](https://zenodo.org/badge/93002532.svg)](https://zenodo.org/badge/latestdoi/93002532)
[![Join the chat at https://gitter.im/GIBBONchat/Lobby](https://badges.gitter.im/GIBBONchat/Lobby.svg)](https://gitter.im/GIBBONchat/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![License](https://img.shields.io/badge/License-GNU_GPLv3-orange.svg)](https://github.com/gibbonCode/GIBBON/blob/master/LICENSE)

# Table of contents
- [Project Summary](#Summary)  
- [Application highlights](#Application)  
- [Installation](#Installation)  
- [Getting started](#Start)
- [Testing](#Test)
- [License](#License)  
- [Contributing](#Contributing)  
- [Code of conduct](#CodeOfConduct)  
- [Road Map](#RoadMap)  

# Project summary <a name="Summary"></a>
GIBBON (The Geometry and Image-Based Bioengineering add-On) is an open-source MATLAB toolbox by [Kevin M. Moerman](kevimoerman.org) and includes an array of image and geometry visualization and processing tools and is interfaced with free open source software such as [TetGen](http://wias-berlin.de/software/tetgen/), for robust tetrahedral meshing, and [FEBio](http://febio.org/) for finite element analysis. The combination provides a highly flexible image-based modelling environment and enables advanced inverse finite element analysis.

[![GIBBON overview](docs/html/GIBBON_overview.jpg)](https://gibboncode.org)

## Application highlights <a name="Application"></a>  
- [Segmentation](#Segmentation)  
- [Computer Aided Design (CAD) tools](#CAD)  
- [Surface meshing tools](#SurfaceMeshing)  
- [Volumetric meshing](#Meshing)  
- [Lattice structures](#Lattice)
- [Finite element analysis](#FEA)
- [Visualization](#Visualization)    

### Segmentation  <a name="Segmentation"></a>    
In patient-, or subject-specific biomechanics, the geometry information is often derived from image data (e.g. Magnetic Resonance Imaging (MRI)). GIBBON offers image filtering and smoothening methods, and has a graphical user interface for 3D image segmentation (see `HELP_imx.m`). The segmented image data can be converted to 3D surface models which can be meshed for FEA.  

<img src="docs/img/imx_demo.gif" width="70%">   
<!-- ![Segmentation](docs/img/imx_demo.gif) -->

### Computer Aided Design (CAD) tools <a name="CAD"></a>  
Sometimes geometry is instead imported or designed. Using GIBBON geometry can be imported from common mesh based CAD files (such as STL, see `HELP_import_STL`). However, several CAD style commands have also been implemented in GIBBON such as polygon rounding, revolution (see `HELP_polyRevolve`), extrusion (see `HELP_polyExtrude`), and sweeping and lofting (see `HELP_polyLoftLinear` and `HELP_sweepLoft`).   

### Surface meshing tools<a name="SurfaceMeshing"></a>   
2D multi-region triangular meshing (see for instance `HELP_regionTriMesh2D` and `HELP_multiRegionTriMeshUneven2D`), resampling meshes geodesically (see `DEMO_geodesic_remeshing`), smoothening (`DEMO_surface_smooth_methods`), and surface mesh refinement (see `HELP_subtri`, `HELP_subTriDual` and `HELP_subQuad`), mesh type conversions (see `HELP_tri2quad`, `HELP_quad2tri`), and mesh dual computation (see `HELP_patch_dual`). Geometries can also be exported to the STL format e.g. for computer aided manufacture and 3D printing.   

### Volumetric meshing <a name="Meshing"></a>   
Tetrahedral meshing (and constrained Delaunay tessellation) of multi-region domains is enabled through an interface with the [TetGen](http://wias-berlin.de/software/tetgen/) [@Si2015] package (`HELP_runTetGen`). Hexahedral meshes for some geometry types can be directly coded (e.g. sphere `HELP_hexMeshSphere`, boxes `HELP_hexMeshBox` and lattices `HELP_element2HexLattice`). For general input surfaces multi-region mixed tetrahedral-hexahedral meshing is also available (`DEMO_MixedTetHexMeshing`).   

<img src="docs/img/bunnyMesh.gif" width="70%">   
<!-- ![Tetrahedral meshing](docs/img/bunnyMesh.gif) -->

### Lattice structures <a name="Lattice"></a>
Various lattice structure tools have been implemented. One method to generate surface geometry for lattices is the use of triply-periodic functions (see `HELP_triplyPeriodicMinimal` and `DEMO_FEBio_trabeculae_compression`). Functions to convert element descriptions, such as tetrahedral and hexahedral elements, to lattice structures have also been implemented (`HELP_element2lattice` and `HELP_element2HexLattice`). These allow for the creation of 3D boundary conforming lattice structures on arbitrary input geometry. Exporting of hexahedral elements is also supported allowing for FEA on the created lattice structures.

### Finite Element Analysis <a name="FEA"></a>
GIBBON interfaces with the free software [FEBio](http://febio.org/) for FEA (source code available on FEBio website). GIBBON can be used as a pre- and post- processor for FEA as it enables the code based development of meshes, boundary conditions, and input files. FEBio files can be directly exported based on dedicated MATLAB® structures. Furthermore, GIBBON can be used to start and control FEBio simulations (see also the `runMonitorFEBio` function). As such, iterative and inverse FEA (e.g. based on MATLAB® optimization routines) is also enabled. All `DEMO_FEBio_...` files contain FEBio examples, `DEMO_FEBio_iFEA_uniaxial_01` is a simple inverse FEA example.
The image below is for large strain analysis of a twisting bar and stems from the demo `DEMO_FEBio_bar_twist`. Other `DEMO_FEBio_...` files cover uni-axial tension/compression, bending, indentation, viscoelastic analysis, contact and indentation problems, multi-generational materials for pre-load analysis.   

<img src="docs/img/barTwist.gif" width="70%">   
<!-- ![Large strain analysis](docs/img/barTwist.gif) -->

### Visualization <a name="Visualization"></a>    
GIBBON expands the standard MATLAB® visualization capabilities by adding 3D image and voxel visualization (see `HELP_im2patch` and `HELP_sliceViewer`), meshed geometries (`HELP_gpatch`), finite element models (`HELP_element2patch`), and colormapped vector data (`HELP_quiverVec`), and all visualization methods enable multiple colormaps to be used in each figure or axis window. Furthermore GIBBON offers a custom figure window `cFigure` containing 3D rotation options (`HELP_vcw`) that mimic CAD behavior of 3D scene rendering, and high quality figure exporting options (`HELP_efw`). Advanced graphics animation creation and exporting capabilities through a figure window based GUI are also enabled (see `HELP_anim8`).   

# Installation <a name="Installation"></a>  

### 1. Installing 3rd party packages
Skip this step if finite element analysis and figure exporting are not required.

| Package | Description | Included? | Download |
|:--|:--|:--:|--:|
|[__FEBio__](https://febio.org) <br/> <img src="docs/img/logos/febioLogo.png" href="https://febio.org" alt="FEBIO" width="100%">|FEBio is a finite element solver and is used in GIBBON for all finite element analysis. Use of FEBio is featured in the many `DEMO_FEBio...` files. FEBio version 2.5.0 or newer is recommended. |__No__|[__FEBio website__](https://febio.org) |
|[__export_fig__](https://github.com/altmany/export_fig) <br/> <img src="docs/img/logos/export_fig_logo.jpg" href="https://github.com/altmany/export_fig" alt="export_fig" width="100%">| <br/> `export_fig` helps to export publication quality images (e.g. .png, .jpg, .pdf, .eps), in GIBBON it is integrated in the export figure widget `efw` to export such images from the `cFigure` window directly. `export_fig` is also used for exporting images for creation of .gif animations with the GIBBON `anim8` function. |__No__|[__Get via GitHub__](https://github.com/altmany/export_fig) <br/> <br/> [__Download zip__](https://github.com/altmany/export_fig/archive/master.zip)|
|<br/> [__TetGen__]() <br/> <img src="docs/img/logos/tetgenLogo.gif" href="http://wias-berlin.de/software/tetgen/" alt="TetGen" width="100px">| <br/> Is used for tetrahedral meshing (and possibly constrained 3D Delaunay tessellation). See for instance `HELP_runTetGen.m`. |__Yes__| For other versions: [__TetGen website__](http://wias-berlin.de/software/tetgen/)|

### 2. Run `installGibbon.m`
By running `installGibbon.m` the GIBBON, FEBio, and export_fig path definitions will be added and saved to MATLAB. The help and documentation will also be integrated. Once finished you will be asked to __restart MATLAB__. `installGibbon.m` can be found in the main GIBBON folder.

# Getting started <a name="Start"></a>

### Access the integrated help
* To access the help documentation from MATLAB click on the HELP browser then click o `GIBBON toolbox` under `Supplemental Software` as shown below. This will open the toolbox help and documentation which is now searchable and integrated just like the rest of MATLAB's help and documentation.  

<img src="docs/gif_helpSearch.gif" alt="Help integration" width="100%">

### Where to find functions and the executable help and demo files
* The `lib` folder contains all GIBBON's functions and the `lib_ext` contains "external functions" i.e. functions developed by others included with GIBBON. The `docs` folder contains the help&documentation, and demo files which when "published" (using MATLAB's publish functionality) create the .html documentation files (found in `docs/html`) which are integrated in MATLAB.  

* The source for the help information for any function `functionName` is named `HELP_functionName`, and  the source for demos have `DEMO_` as part of the name. Therefore if one is interested in reproducing or starting off from codes in the help and documentation simply start typing code names starting in `HELP_` or `DEMO_` in the MATLAB command window, e.g. `HELP_ind2patch` can be used to generate the help information for the `ind2patch` function. Users can start editing the file by typing `open HELP_ind2patch` in the command window. By publishing (MATLAB publish functionality) the HELP_ or DEMO_ files .html files are created in the `docs\html` folder. As such if users alter/contribute code in the `lib` folder and generate associated `HELP_` or `DEMO_` files, new help and documentation is added. For new help and documentation to become known and visible to MATLAB run the `createHelpDemoDocumentation` function and restart MATLAB.  

* Many of the `DEMO_` files focus on the use of FEBio. The demo `DEMO_FEBio_block_uniaxial_compression` for instance features a simple cube that undergoes a 30% compression. Other demos focus on different load types, single versus multi-step analysis, different materials and inverse analysis (e.g. `DEMO_FEBio_iFEA_uniaxial_01`).

# Testing <a name="Test"></a>
GIBBON's core functionality can be tested by running `testGibbon('all','test');`. Use `testGibbon('demo','test');` or `testGibbon('help','test');` for running the demo or help files only. Use the `'publish'` option to test and publish the output to the integrated help and documentation`testGibbon('all','publish');`.   
GIBBON is currently developed and tested using the most recent version of MATLAB (or the latest pre-release) and has been tested on Windows 10, Ubuntu 14.10/16.04/17.10, and Mac OS. Most of GIBBON's functionality is compatible with older MATLAB versions, especially MATLAB R2014a and newer (Delaunay tessellation and toolbox help integration are amongst things that have undergone large change). Please inform the developers (or open an issue) if a particular function does not work for your MATLAB environment. It is likely that codes can be altered to work for your version.    
A large portion of GIBBON's functionality does not rely on special MATLAB toolboxes. However some functions do. Here is a list of toolboxes which appear to be used in GIBBON:
* Image Processing Toolbox
* Bioinformatics Toolbox
* Statistics and Machine Learning Toolbox
* Computer Vision System Toolbox
* Neural Network Toolbox
* Symbolic Math Toolbox
* Curve Fitting Toolbox
* Parallel Computing Toolbox
* Mapping Toolbox    

Geodesic remeshing (`DEMO_geodesic_remeshing`) currently only works for Windows OS, for other operational systems special mex file compilation is required or an alternative (experimental function) is used.

# License <a name="License"></a>
[![License](https://img.shields.io/badge/License-GNU_GPLv3-orange.svg)](https://github.com/gibbonCode/GIBBON/blob/master/LICENSE)

# Contributing <a name="Contributing"></a>
See [CONTRIBUTING](CONTRIBUTING.md)    
[![Join the chat at https://gitter.im/GIBBONchat/Lobby](https://badges.gitter.im/GIBBONchat/Lobby.svg)](https://gitter.im/GIBBONchat/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# Code of conduct <a name="CodeOfConduct"></a>
See [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md)

# Roadmap <a name="RoadMap"></a>
Coming soon
