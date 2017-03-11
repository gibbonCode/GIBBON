
![](http://www.gibboncode.org/wp-content/uploads/2014/05/header1.png)

## Project summary
GIBBON (The Geometry and Image-Based Bioengineering add-On) is an open-source MATLAB toolbox by [Kevin M. Moerman](kevimoerman.org) and includes an array of image and geometry visualization and processing tools and is interfaced with free open source software such as [TetGen](http://wias-berlin.de/software/tetgen/), for robust tetrahedral meshing, and [FEBio](http://febio.org/) for finite element analysis. The combination provides a highly flexible image-based modelling environment and enables advanced inverse finite element analysis.

![](http://www.gibboncode.org/wp-content/uploads/2014/05/Picture1.png)

## Installation
See also: http://www.gibboncode.org/installation/

Follow the steps below. Note: A mltbx file is no longer provided (file changed with each release and so it is not very compatible with GIT/SVN based version control).

1. Add the toolbox folder (e.g. for SVN the patch typically ends in `..\gibbon\trunk`) including subdirectories to the MATLAB path (see MATLAB home tab and the `Add path` button, or see MATLAB help). To permanently add GIBBON to the path save the added path definitions. 
2. To integrate the GIBBON help and demonstrations in MATLAB run the function `createHelpDemoDocumentation` (found in the createHelpDoc folder, after adding GIBBON to the path it can be run from the command window).
3. Restart MATLAB to allow it to update help and documentation definitions.
3. To access the help documentation from MATLAB click on the HELP browser then click `Supplemental Software` as shown below.
![](http://www.gibboncode.org/wp-content/uploads/2014/05/helpFront1.png)
This will open the following window showing the toolbox help and documentation integrated into MATLAB:
![](http://www.gibboncode.org/wp-content/uploads/2014/05/helpFront2.png)

GIBBON is now installed and the help and documentation is integrated allong side your default MATLAB help and documentation. The help and documentation is also searchable. 

## Setting up third party packages
* **FEBio** FEBio is a finite element package. FEBio is not provided with GIBBON. Install a desired release (see: http://febio.org/) and change the path name to FEBio in the configuration file `FEBioPath.txt`found in the `config` folder. However the config path can also be ignored by always directly specifying the location in the code (see FEBio related demos).
* **TetGen** TetGen is a tetrahedral meshing package. TetGen is provided with GIBBON and can be found in `...\GIBBON\lib_ext\tetGen`. No special treatment should be required to run and use TetGen. If an alternative release is required visit the [TetGen website](http://wias-berlin.de/software/tetgen/), and replace the existing files as desired. 
* **export_fig** export_fig is not included with GIBBON. These figure exporting MATLAB tools (used by the GIBBON 'efw' function) can be found on [GitHub](https://github.com/altmany/export_fig) or on [the MathWorks file exchange](http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig). 

## Getting started
* Study the GIBBON folder structure. For example, nottice how the `lib` folder contains all GIBBON's functions and that `lib_ext` contains "external functions" i.e. functions developed by others included with GIBBON. The `html` folder contains the help and documentation .html files which are integrated in MATLAB. 
* **Help and documentation** see the installation instructions on accessing the help and documentation. Through the integrated help and documentation the user can explore variations function descriptions and also demo entries. The codes that generate all the help and documention can be found in GIBBON's main folder. **The source for the help information for any function `functionName` is named `HELP_functionName`**, and  **The source for demos have `DEMO_` as part of the name**. Therefore if one is interested in reproducing or starting off from codes in the help and documentation simply start typing code names starting in `HELP_` or `DEMO_` in the MATLAB command window, e.g. `HELP_ind2patch` can be used to generate the help information for the `ind2patch` function. Users can start editing the file by typing `open HELP_ind2patch` in the command window. By publishing (MATLAB publish functionality) the HELP_ or DEMO_ files .html files are created in the `html` folder. As such if users alter/contribute code in the `lib` folder and generate associated `HELP_` or `DEMO_` files new help and documentation is added. For new help and documentation to become known and visible to MATLAB run the `createHelpDemoDocumentation` function and restart MATLAB. 

## [License (BSD-3-Clause)](https://github.com/Kevin-Mattheus-Moerman/GIBBON/blob/master/LICENSE)

##Contributing
Coming soon

##Roadmap
Coming soon
