---
layout: page
title: "Installation"
description: "How to install GIBBON"
header-img: "img/home-bg.jpg"
---

#### Table of content
* [Getting GIBBON](#GettingGIBBON)
  * [GitHub and GIT](#git)
  * [Direct download](#direct)
* [Installation](#Installation)
  * [Setting up GIBBON](#setup)
  * [Setting up third party packages](#3rdparties)

***

# Getting GIBBON <a name="GettingGIBBON"></a>
To stay up to date with new developments it is recommended to keep up with the latest version on the GitHub repository by setting up git. Alternatively users may download a zipfile directly.

### GitHub and GIT <a name="git"></a>
GIBBON is hosted on [this](https://github.com/gibbonCode/GIBBON) GitHub repository.  
[![Download](img/gibbon_github.png){:height="150px"}](https://github.com/gibbonCode/GIBBON)

```
git init
git clone https://github.com/gibbonCode/GIBBON.git
```

### Direct download <a name="direct"></a>
Alternatively download the [.zip file](https://github.com/gibbonCode/GIBBON/archive/master.zip) directly

[![Download](img/gibbonDownLoadBox.png){:height="150px"}](https://github.com/gibbonCode/GIBBON/archive/master.zip)

# Installation <a name="Installation"></a>  

#### 1. __Installing 3rd party packages__
* **FEBio** is the finite element solver used by GIBBON. FEBio is not provided with GIBBON user need to download a desired release from the [FEBio website](http://febio.org/) and install it.    
* **export_fig**  is a MATLAB library which GIBBON uses for exporting figures. The export_fig library is available through its [gitHub](https://github.com/altmany/export_fig) page (or the [MathWorks file exchange](http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)).
Already bundled with GIBBON:
* **TetGen** is a tetrahedral meshing package which is provided with GIBBON. If an alternative release is required visit the [TetGen website](http://wias-berlin.de/software/tetgen/), and replace the existing files (in `...\GIBBON\lib_ext\tetGen)` as desired.

#### 2. __Run `installGibbon.m`__
The file can be found here `...\GIBBON\installGibbon.m` and will add the GIBBON, FEBio, and export_fig path definitions to MATLAB, and will also integrate the help and documentation.

#### 3. __Restart MATLAB__
