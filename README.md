
<img src="docs/html/gibbonLogo.jpg" href="https://gibboncode.org" alt="GIBBON" width="100%">   

[![Citing](https://zenodo.org/badge/DOI/10.5281/zenodo.44404.svg)](https://doi.org/10.5281/zenodo.44404) [![Join the chat at https://gitter.im/GIBBONchat/Lobby](https://badges.gitter.im/GIBBONchat/Lobby.svg)](https://gitter.im/GIBBONchat/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![License](https://img.shields.io/badge/License-GNU_GPLv3-orange.svg)](https://github.com/gibbonCode/GIBBON/blob/master/LICENSE)

# Table of contents
- [Project Summary](#Summary)  
- [Installation](#Installation)  
- [Getting started](#Start)  
- [License](#License)  
- [Contributing](#Contributing)  
- [Code of conduct](#CodeOfConduct)  
- [Road Map](#RoadMap)  

## Project summary <a name="Summary"></a>
GIBBON (The Geometry and Image-Based Bioengineering add-On) is an open-source MATLAB toolbox by [Kevin M. Moerman](kevimoerman.org) and includes an array of image and geometry visualization and processing tools and is interfaced with free open source software such as [TetGen](http://wias-berlin.de/software/tetgen/), for robust tetrahedral meshing, and [FEBio](http://febio.org/) for finite element analysis. The combination provides a highly flexible image-based modelling environment and enables advanced inverse finite element analysis.

<img src="docs/html/GIBBON_overview.jpg" href="https://gibboncode.org" alt="GIBBON overview" width="100%">

# Installation <a name="Installation"></a>  
### 1. Installing 3rd party packages
Skip this step if finite element analysis and figure exporting are not required.

| Package | Purpose | Included? | Download |
|:--|:--|:--:|--:|
|[__FEBio__](https://febio.org) <br/> <img src="docs/img/logos/febioLogo.png" href="https://febio.org" alt="FEBIO" width="100%">|FEBio is a finite element solver and is used in GIBBON for all finite element analysis. Use of FEBio is featured in many of the `DEMO_FEBio...` files |__No__|[__FEBio website__](https://febio.org) |
|[__export_fig__](https://github.com/altmany/export_fig) <br/> <img src="docs/img/logos/export_fig_logo.jpg" href="https://github.com/altmany/export_fig" alt="export_fig" width="100%">| <br/> `export_fig` helps to export publication quality images (e.g. .png, .jpg, .pdf, .eps), in GIBBON it is integrated in the export figure widget `efw` to export such images from the `cFigure` window directly. `export_fig` is also used for exporting images for creation of .gif animations with the GIBBON `anim8` function |__No__|[__Get via GitHub__](https://github.com/altmany/export_fig) <br/> <br/> [__Download zip__](https://github.com/altmany/export_fig/archive/master.zip)|
|<br/> [__TetGen__]() <br/> <img src="docs/img/logos/tetgenLogo.gif" href="http://wias-berlin.de/software/tetgen/" alt="TetGen" width="100px">| <br/> Is used for tetrahedral meshing (and possibly constrained 3D Delaunay tessellation). See for instance `HELP_runTetGen.m`|__Yes__| For other versions: [__TetGen website__](http://wias-berlin.de/software/tetgen/)|

### 2. Run `installGibbon.m`
By running `installGibbon.m` the GIBBON, FEBio, and export_fig path definitions will be added and saved to MATLAB. The help and documentation will also be integrated. Once finished you will be asked to __restart MATLAB__. `installGibbon.m` can be found in the main GIBBON folder.

## Getting started <a name="Start"></a>

##### Access the integrated help
* To access the help documentation from MATLAB click on the HELP browser then click o `GIBBON toolbox` under `Supplemental Software` as shown below. This will open the toolbox help and documentation which is now searchable and integrated just like the rest of MATLAB's help and documentation.  

<img src="docs/gif_helpSearch.gif" alt="Help integration" width="100%">

##### Where to find functions and the executable help and demo files
* The `lib` folder contains all GIBBON's functions and the `lib_ext` contains "external functions" i.e. functions developed by others included with GIBBON. The `docs` folder contains the help&documentation, and demo files which when "published" (using MATLAB's publish functionality) create the .html documentation files (found in `docs/html`) which are integrated in MATLAB.  

* The source for the help information for any function `functionName` is named `HELP_functionName`, and  the source for demos have `DEMO_` as part of the name. Therefore if one is interested in reproducing or starting off from codes in the help and documentation simply start typing code names starting in `HELP_` or `DEMO_` in the MATLAB command window, e.g. `HELP_ind2patch` can be used to generate the help information for the `ind2patch` function. Users can start editing the file by typing `open HELP_ind2patch` in the command window. By publishing (MATLAB publish functionality) the HELP_ or DEMO_ files .html files are created in the `doc\html` folder. As such if users alter/contribute code in the `lib` folder and generate associated `HELP_` or `DEMO_` files, new help and documentation is added. For new help and documentation to become known and visible to MATLAB run the `createHelpDemoDocumentation` function and restart MATLAB.  

* Many of the `DEMO_` files focus on the use of FEBio. The demo `DEMO_FEBio_block_uniaxial_compression` for instance features a simple cube that undergoes a 30% compression. Other demos focus on different load types, single versus multi-step analysis, different materials and inverse analysis (e.g. `DEMO_FEBio_iFEA_uniaxial_01`).

## License <a name="License"></a>
[![License](https://img.shields.io/badge/License-GNU_GPLv3-orange.svg)](https://github.com/gibbonCode/GIBBON/blob/master/LICENSE)

## Contributing <a name="Contributing"></a>
See [CONTRIBUTING](CONTRIBUTING.md)    
[![Join the chat at https://gitter.im/GIBBONchat/Lobby](https://badges.gitter.im/GIBBONchat/Lobby.svg)](https://gitter.im/GIBBONchat/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## Code of conduct <a name="CodeOfConduct"></a>
See [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md)

## Roadmap <a name="RoadMap"></a>
Coming soon
