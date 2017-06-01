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
GIBBON is hosted on [this](https://github.com/Kevin-Mattheus-Moerman/GIBBON) GitHub repository.  
[![Download](/img/gibbon_github.png){:height="150px"}](https://github.com/Kevin-Mattheus-Moerman/GIBBON)

```
git init
git clone https://github.com/Kevin-Mattheus-Moerman/GIBBON.git
```

### Direct download <a name="direct"></a>
Alternatively download the [.zip file](https://github.com/Kevin-Mattheus-Moerman/GIBBON/archive/master.zip) directly

[![Download](/img/gibbonDownLoadBox.png){:height="150px"}](https://github.com/Kevin-Mattheus-Moerman/GIBBON/archive/master.zip)

# Installation
### Setting up GIBBON <a name="setup"></a>
Follow the steps below.

1. __Add the toolbox__ (e.g. `..\gibbon`) including subdirectories to the MATLAB path (see MATLAB home tab and the `Add path` button, or see MATLAB help). To permanently add GIBBON to the path save the added path definitions.
![](/img/gif_addToPath.gif){:width="100%"}
2. __Integrate the GIBBON help and demos__ in MATLAB by running the function `createHelpDemoDocumentation` (found in the createHelpDoc folder, after adding GIBBON to the path it can be run directly from the command window).
![](/img/gif_createHelpDoc.gif){:width="100%"}
3. __Restart MATLAB__ to allow it to update help and documentation definitions.
4. To access the help documentation from MATLAB click on the HELP browser click the GIBBON link under `Supplemental Software` as shown below. This will open the toolbox help and documentation which is now searchable and integrated into MATLAB.
![](/img/gif_helpSearch.gif){:width="100%"}

### Setting up third party packages <a name="3rdparties"></a>
* **FEBio** FEBio is a finite element package. FEBio is not provided with GIBBON. Install a desired release (see: http://febio.org/) and change the path name to FEBio in the configuration file `FEBioPath.txt`found in the `config` folder. However the config path can also be ignored by always directly specifying the location in the code (see FEBio related demos).
* **TetGen** TetGen is a tetrahedral meshing package. TetGen is provided with GIBBON and can be found in `...\GIBBON\lib_ext\tetGen`. No special treatment should be required to run and use TetGen. If an alternative release is required visit the [TetGen website](http://wias-berlin.de/software/tetgen/), and replace the existing files as desired.
* **export_fig** export_fig is not included with GIBBON. These figure exporting MATLAB tools (used by the GIBBON `efw` function for exporting figures and by the GIBBON `anim8` to export .gif animations) can be found on [GitHub](https://github.com/altmany/export_fig) or on [the MathWorks file exchange](http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig).
