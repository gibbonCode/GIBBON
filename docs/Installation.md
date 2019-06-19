---
layout: page
title: "Installation"
description: "How to install GIBBON"
header-img: "img/home-bg.jpg"
ordernumber: 1
---
# Installation

#### Table of content
* [Download](#Download)  
* [Installation](#Installation)
* [Setting up Git](#Git)

***

## Download <a name="Download"></a>
There are two ways to download the software, i.e. by using Git and GitHub or by downloading a zip file:   

|[Setup GIT/GitHub (recommended)](https://github.com/gibbonCode/GIBBON)|[Download .zip file directly](https://github.com/gibbonCode/GIBBON/archive/master.zip)|
|:-|-:|
|[![Github repository](img/gibbon_github.png){:height="125px"}](https://github.com/gibbonCode/GIBBON)|[![Download](img/gibbonDownLoadBox.png){:height="125px"}](https://github.com/gibbonCode/GIBBON/archive/master.zip)|  

Downloading the zip file will create a local copy on your computer which is fully decoupled from the new developments taking place. It is therefore recommended to use Git (a software version control system) instead, which allows you to 1) pull in the latest changes as they become available, 2) potentially join in with development and push changes to the project. Read more on how to use git [below](#Git).

__Make sure the GIBBON folder is in a location where you have full read/write permissions!__

## Installation <a name="Installation"></a>  
The steps below guide you through a streamlined installation procedure using the `installGibbon.m` function<sup>\*</sup>.   
\*<sub>If you prefer manual installation do the following: 1) Add the GIBBON folder (with subfolders) to the path and save the path definitions, 2) Run `createHelpDemoDocumentation.m` to integrate the help and documentation, 3) For the 3rd party packages: 3a) Add the `export_fig` folder to the path and save the path definitions, 3b) Go to the config folder in _../GIBBON/config_ and edit the _FEBioPath.txt_ file to contain the full path to the FEBio executable </sub>

#### 1. Installing 3rd party packages
Skip this step if finite element analysis and figure exporting are not required.

| Package | Purpose | Included? | Download |
|:--|:--|:--:|--:|
|[__FEBio__](https://febio.org) <br/> [![FEBio](/img/logos/febioLogo.png){:height="100px"}](https://febio.org)|FEBio is a finite element solver and is used in GIBBON for all finite element analysis. Use of FEBio is featured in many of the `DEMO_febio...` files |__No__|[__FEBio website__](https://febio.org) |
|[__export_fig__](https://github.com/altmany/export_fig) <br/> [![export_fig](/img/logos/export_fig_logo.jpg){:height="100px"}](https://github.com/altmany/export_fig)| <br/> `export_fig` helps to export publication quality images (e.g. .png, .jpg, .pdf, .eps), in GIBBON it is integrated in the export figure widget `efw` to export such images from the `cFigure` window directly. `export_fig` is also used for exporting images for creation of .gif animations with the GIBBON `anim8` function |__No__|[__Get via GitHub__](https://github.com/altmany/export_fig) <br/> <br/> [__Download zip__](https://github.com/altmany/export_fig/archive/master.zip)|
|<br/> [__TetGen__]() <br/> [![tetGen](/img/logos/tetgenLogo.gif){:height="100px"}](https://wias-berlin.de/software/tetgen/)| <br/> Is used for tetrahedral meshing (and possibly constrained 3D Delaunay tessellation). See for instance `HELP_runTetGen.m`|__Yes__| For other versions: [__TetGen website__](https://wias-berlin.de/software/tetgen/)|

#### 2. Run `installGibbon.m`
By running `installGibbon.m` the GIBBON, FEBio, and export_fig path definitions will be added and saved to MATLAB. The help and documentation will also be integrated. Once finished you will be asked to __restart MATLAB__. `installGibbon.m` can be found in the main GIBBON folder.

---

## Setting up and using Git <a name="Git"></a>
Below is a very basic set of instructions for using Git (which skips over advanced Git features). If you are new to Git you may want to search for tutorials online, also check out these [learning resources](https://help.github.com/articles/git-and-github-learning-resources/) and this [10 min. GIT tutorial](https://try.github.io/levels/1/challenges/1).   

#### Setting up the local repository
Make sure you have git [installed](https://git-scm.com/downloads), browse to a suitable folder location (one where you have read/write permissions) and run:
```
git init   
git clone https://github.com/gibbonCode/GIBBON.git   
```  
This will create a local clone on your system. Git will keep track of new developments on the remote version (GitHub repository) version and also changes you make locally.

#### Check what's been changed on the remote or your local copy
To view changes, e.g. new files, modifications, rename events, and deletions, use:
```
git status
```

#### Removing the changes you've made
If you want to remove all changes you've made to GIBBON locally (e.g. to avoid merge conflicts) you can force Git to remove your changes and revert back to the version you cloned by using:
```
git reset --hard HEAD
```
Note the above operation only removes your changes. To get the latest version from the remote repository you still need to do `git pull`.

#### Pulling in new GIBBON changes from the remote to your local copy
To pull changes available on the remote version (the GitHub repository). Run the following:
```
git pull
```
Pulling in changes is quick and easy if you haven't made conflicting changes, i.e. if you made changes to files locally that are also changed on the remote version. If changes have occurred there may be so called "merge conflicts" (note that changes may naturally in your GIBBON folder, e.g. as part of installation, or by making minor or even accidental changes to files as one uses GIBBON). Users comfortable with Git may be able to easy resolve such conflicts (by merging the changes from the remote and the local). If the changes listed are not relevant and can be reverted you can do a reset as described in the previous step. If this is done you may have to run `installGibbon.m` again (as some of the changes may be due to documentation integration on your system). You can verify everything is up to date with the remote repository by running `git status`.
