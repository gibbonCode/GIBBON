# Roadmap

This roadmap is a place to start to investigate how you can contribute to the GIBBON project. Check out the different milestones listed below, but feel free to explore the issues by label too.

Please check out our [contribution guidelines](https://github.com/gibbonCode/GIBBON/blob/master/CONTRIBUTING.md) and [code of conduct](https://github.com/gibbonCode/GIBBON/blob/master/CODE_OF_CONDUCT.md) to help you get started, and the [README](https://github.com/gibbonCode/GIBBON/blob/master/README.md) for an overview if you haven't read it yet!

## Short term milestones

### 1. Implement medical device design optimization functionality
* Enable iterative optimization and topology optimization 
* Add demos covering this funcitonality

### 2. Lattice structure geometry creation for FEA simulation and 3D printing
* Element to lattice conversion
* Dual based lattice creation
* Direct coding of diamond and octet-truss lattices
* Demos for FEA of lattices and STL creation

### 3. Functionality for branched blood vessel architecture meshing
* Avoid self-intersection for polyTube and related functions
* Finalize convex-hull based branch meshing

### 4. Implement GEOGRAM remeshing functionality
* The open source project GEOMGRAM (you can test functionality here: https://members.loria.fr/Bruno.Levy/GEOGRAM/geobox.html) contains very efficient remeshing algorithms to produce high quality triangulations in input geometry. Making these accessible from GIBBON would make it a powerful solution for remeshing which is currently rather inefficient.

## Long term milestones

### 1. Enable import of FEBio XML files

### 2. Enable import of FEBio .XPLT files

### 3. Add visualization GUI for FEBio based on XPLT output

### 4. Add more ABAQUS demos

## Implementations in other languages

### Julia implementation
A Julia implementation of GIBBON is planned for 2021. Developments will take place here: https://github.com/gibbonCode/GIBBON.jl

### Support Octave implementation
Preliminary investigations have shown implementation of GIBBON in Octave is in principle possible. However, it does come with reduced performance, especially in terms of mesh visualization which major drawback at the moment. If more efficient visualization solutions can be found the Octave implementation can be supported.
