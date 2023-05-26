# About GIBBON

# Project summary
GIBBON (The Geometry and Image-Based Bioengineering add-On) is an open-source MATLAB toolbox by [Kevin M. Moerman](https://kevinmoerman.org) and includes an array of image and geometry visualization and processing tools and is interfaced with free open source software such as [TetGen](http://wias-berlin.de/software/tetgen/), for robust tetrahedral meshing, and [FEBio](http://febio.org/) for finite element analysis. The combination provides a highly flexible image-based modelling environment and enables advanced inverse finite element analysis.

![GIBBON overview](/html/GIBBON_overview.jpg){:width="100%"}

## Application highlights <a name="Application"></a>  
- [Segmentation](#Segmentation)  
- [Computer Aided Design (CAD) tools](#CAD)  
- [Surface meshing tools](#SurfaceMeshing)  
- [Volumetric meshing](#Meshing)  
- [Lattice structures](#Lattice)
- [Finite element analysis](#FEA)
- [Visualization](#Visualization)    

### Segmentation  <a name="Segmentation"></a>    
GIBBON offers image filtering and smoothing methods, and has a graphical user interface for 3D image segmentation (`HELP_imx.m`). The segmented image data can be converted to 3D surface models (`DEMO_imx_levelset_surface_compare`) which can be meshed for FEA (`HELP_runTetGen`).   

<div>
<img src="/img/imx_demo.gif" width="50%">   
<img src="/img/footrevolve.gif" width="25%">   
<img src="/img/elefootSurfs.png" width="15%">   
</div>

### Computer Aided Design (CAD) tools <a name="CAD"></a>  
Using GIBBON, geometry can be imported from common mesh based CAD files (such as STL, `HELP_import_STL`). For generating geometries within MATLAB®, GIBBON also provides several CAD-style commands such as polygon rounding (`HELP_filletCurve`), revolution (`HELP_polyRevolve`), extrusion (`HELP_polyExtrude`), and sweeping and lofting (`HELP_polyLoftLinear` and `HELP_sweepLoft`). Simple geometries such as spheres (`HELP_geoSphere`), boxes (`HELP_quadBox`), platonic solids (`HELP_platonic_solid`), and rhombic dodecahedra (`HELP_rhombicDodecahedron`) can also be directly created using GIBBON.  

<div>
<img src="/img/gallery/sweepLoft1.gif" width="40%">
<img src="/img/gallery/stentHex.gif" width="40%">   
</div>

### Surface meshing tools<a name="SurfaceMeshing"></a>   
2D multi-region triangular meshing (e.g. `HELP_regionTriMesh2D` and `HELP_multiRegionTriMeshUneven2D`), resampling meshes geodesically (`DEMO_geodesic_remeshing`), smoothing (`DEMO_surface_smooth_methods`), and surface mesh refinement (e.g. `HELP_subtri`, `HELP_subTriDual` and `HELP_subQuad`), mesh type conversions (e.g. `HELP_tri2quad`, `HELP_quad2tri`), and mesh dual computation (`HELP_patch_dual`). Geometries can also be exported to the STL format e.g. for computer aided manufacture and 3D printing.

<div>
<img src="/img/gallery/im_surf_refine1.gif" width="20%">   
<img src="/img/gallery/elephantDance.gif" width="30%">   
<img src="/img/gallery/regionMesh1.gif" width="30%">   
</div>

### Volumetric meshing <a name="Meshing"></a>   
Tetrahedral meshing (and constrained Delaunay tessellation) of multi-region domains is enabled through an interface with the [TetGen](http://wias-berlin.de/software/tetgen/) package (`HELP_runTetGen` and `HELP_constrainedDelaunayTetGen`). Hexahedral meshes for some geometry types can be directly coded (e.g. spheres `HELP_hexMeshSphere`, boxes `HELP_hexMeshBox` and lattices `HELP_element2HexLattice`). For general input surfaces multi-region mixed tetrahedral-hexahedral meshing is also available (e.g. `DEMO_MixedTetHexMeshing`).

<div>
<img src="/img/bunnyMesh.gif" width="35%">
<img src="/img/mixedMesh.png" width="35%">
<img src="/img/gallery/pillowedHexMeshFemur.png" width="70%">
</div>

### Lattice structures <a name="Lattice"></a>
One method to generate surface geometry for lattices is the use of triply-periodic functions (`HELP_triplyPeriodicMinimal`). Functions to convert element descriptions, such as tetrahedral and hexahedral elements, to lattice structures have also been implemented (`HELP_element2lattice` and `HELP_element2HexLattice`). These allow for the creation of 3D boundary conforming lattice structures on arbitrary input geometry. Exporting of hexahedral elements is also supported allowing for FEA on the created lattice structures (`DEMO_febio_0026_hexlattice_compression`).

<div>
<img src="/img/latticeCompress.gif" width="40%">   
<img src="/img/dualClad.gif" width="40%">   
<img src="/img/octet.gif" width="40%">   
<img src="/img/gallery/latticePart.gif" width="40%">   
</div>

### Finite Element Analysis <a name="FEA"></a>
For finite element analysis GIBBON currently links with either the free and open source software [FEBio](http://febio.org/) or with Simulia ABAQUS. Both the FEBio and ABAQUS interface is based on MATLAB® structures. The image below shows the coding of a material section in a MATLAB® structure (top row) and how these components are represented in the input files for FEBio or ABAQUS (bottom row). Through this structure to input file conversion process **any FEBio or ABAQUS functionality can be directly coded in MATLAB®**.
<div>
<img src="/img/FEA_interface_syntax.jpg" width="75%">   
</div>

##### FEBio
GIBBON can be used as a pre- and post- processor for FEBio as it enables code-based development of meshes, boundary conditions, and input files. FEBio files can be directly exported based on dedicated MATLAB® structures (`HELP_febioStruct2xml`). Furthermore, GIBBON can be used to start and control FEBio simulations. As such, iterative and inverse FEA (e.g. based on MATLAB® optimization routines) is also enabled. All `DEMO_febio_...` files are FEBio demos, e.g. `DEMO_febio_0001_cube_uniaxial` is a simple uniaxial loading example, and `DEMO_febio_0042_inverse_FEA_cube_uniaxial` is an example of inverse FEA.    
The image below is for large strain analysis of a twisting bar and stems from the demo `DEMO_febio_0004_beam_twist`. Other demos cover tension, compression, shear, applied forces, applied pressures, applied displacements, bending, poroelasticity, dynamic and viscoelastic analysis, contact and indentation problems, multi-generational materials for pre-load analysis.     

<div>
<img src="/img/gallery/propFlap.gif" width="40%">     
<img src="/img/gallery/mammography.gif" width="40%">
<img src="/img/gallery/softRobot.gif" width="40%">
<img src="/img/gallery/clotSlide.gif" width="40%">
</div>

#### Abaqus
The interface for ABAQUS is a recent development. Users can look at `HELP_abaqusStruct2inp` to study how input files are coded. The demo `DEMO_abaqus_0001_cube_uniaxial` is for uniaxial loading of a cube and steps through geometry creation, setting up the ABAQUS structure, saving the .inp file, running the job, and importing the results for visualization. Data is imported into MATLAB® using `importAbaqusDat` which parses ABAQUS `.DAT` files.

### Visualization <a name="Visualization"></a>    
GIBBON expands the standard MATLAB® visualization capabilities by adding 3D image and voxel visualization (`HELP_im2patch` and `HELP_sv3`), meshed geometries (`HELP_gpatch` and `HELP_meshView`), finite element models (`HELP_element2patch`), and colormapped vector data (`HELP_quiverVec`), and all visualization methods enable multiple colormaps to be used in each figure or axis window. Furthermore GIBBON offers a custom figure window `cFigure` containing 3D rotation options (`HELP_vcw`) that mimic CAD behavior of 3D scene rendering, and high quality figure exporting options (`HELP_efw`). Advanced graphics animation creation and exporting capabilities through a figure window based GUI are also enabled (`HELP_anim8`).   

GIBBON uses a custom figure type called `cFigure`. It contains a white background, is maximized by default, and includes a view control widget (`vcw`) and an export figure widget (`efw`). The `vcw` function enables 3D view control similar to CAD packages. Rotation, zooming, and panning can be performed using a 3 button "mouse" (or equivalent input device), e.g. middle-click rotates the view, left click drags/moves the view, right click allows zooming in and out.

| GIBBON | Classic MATLAB |   
|:--|--:|   
|  |  |   

 ![](/img/gibbonViewControl.gif){:width="100%"}

## Setting up GIBBON
See [installation]({{ site.baseurl }}/Installation/).
