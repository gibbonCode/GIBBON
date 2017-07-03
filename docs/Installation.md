---
layout: page
title: "Installation"
description: "How to install GIBBON"
header-img: "img/home-bg.jpg"
---

#### Table of content
* [Download](#Download)  
* [Installation](#Installation)

***

# Download <a name="Download"></a>

|[Setup GIT/GitHub (recommended)](https://github.com/gibbonCode/GIBBON)|[Download .zip file directly](https://github.com/gibbonCode/GIBBON/archive/master.zip)|
|:-|-:|
|[![Github repository](img/gibbon_github.png){:height="125px"}](https://github.com/gibbonCode/GIBBON)|[![Download](img/gibbonDownLoadBox.png){:height="125px"}](https://github.com/gibbonCode/GIBBON/archive/master.zip)|  

```
git init   
git clone https://github.com/gibbonCode/GIBBON.git   
```  
(__New to GIT?__ see these [learning resources](https://help.github.com/articles/git-and-github-learning-resources/) and this [10 min. GIT tutorial](https://try.github.io/levels/1/challenges/1))

# Installation <a name="Installation"></a>  
The steps below guide you through a streamlined installation procedure using the `installGibbon.m` function<sup>\*</sup>.   
\*<sub>If you prefer manual installation do the following: 1) Add the GIBBON folder (with subfolders) to the path and save the path definitions, 2) Run `createHelpDemoDocumentation.m` to integrate the help and documentation, 3) For the 3rd party packages: 3a) Add the `export_fig` folder to the path and save the path definitions, 3b) Go to the config folder in _../GIBBON/config_ and edit the _FEBioPath.txt_ file to contain the full path to the FEBio executable </sub>

### 1. Installing 3rd party packages
Skip this step if finite element analysis and figure exporting are not required.

| Package | Purpose | Included? | Download |
|:--|:--|:--:|--:|
|[__FEBio__](https://febio.org) <br/> [![FEBio](/img/logos/febioLogo.png){:height="100px"}](https://febio.org)|FEBio is a finite element solver and is used in GIBBON for all finite element analysis. Use of FEBio is featured in many of the `DEMO_FEBio...` files |__No__|[__FEBio website__](https://febio.org) |
|[__export_fig__](https://github.com/altmany/export_fig) <br/> [![export_fig](/img/logos/export_fig_logo.jpg){:height="100px"}](https://github.com/altmany/export_fig)| <br/> `export_fig` helps to export publication quality images (e.g. .png, .jpg, .pdf, .eps), in GIBBON it is integrated in the export figure widget `efw` to export such images from the `cFigure` window directly. `export_fig` is also used for exporting images for creation of .gif animations with the GIBBON `anim8` function |__No__|[__Get via GitHub__](https://github.com/altmany/export_fig) <br/> <br/> [__Download zip__](https://github.com/altmany/export_fig/archive/master.zip)|
|<br/> [__TetGen__]() <br/> [![tetGen](/img/logos/tetgenLogo.gif){:height="100px"}](https://wias-berlin.de/software/tetgen/)| <br/> Is used for tetrahedral meshing (and possibly constrained 3D Delaunay tessellation). See for instance `HELP_runTetGen.m`|__Yes__| For other versions: [__TetGen website__](https://wias-berlin.de/software/tetgen/)|

### 2. Run `installGibbon.m`
By running `installGibbon.m` the GIBBON, FEBio, and export_fig path definitions will be added and saved to MATLAB. The help and documentation will also be integrated. Once finished you will be asked to __restart MATLAB__. `installGibbon.m` can be found in the main GIBBON folder.
