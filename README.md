Dentist2
========
This repository contains tools to help analyze RNA FISH data from large image scans. Example input and output data are available on Dropbox [here](https://www.dropbox.com/sh/djy9mock6vp5j2d/AADqZDtkAEROWK7Jh85oCE68a?dl=0).  

Table of contents
=================
* [Dentist2](#dentist2)
* [Dependencies](#dependencies)
* [Expected input](#expected-input)
* [Quick start](#quick-start)
* [Stitching](#stitching)
* [D2ThresholdGUI](#d2ThresholdGUI)
  * [Main axes](#main-axes)
  * [Thumbnail axes](#thumbnail-axes)
  * [Threshold axes](#threshold-axes)
  * [Masks](#masks)
* [Importing CellPose masks](#importing-cellpose-masks)
* [Description of objects and default parameters](#description-of-objects-and-default-parameters)
* [Troubleshooting](#troubleshooting)

Dependencies
=============
Dentist2 uses the Bio-Formats MATLAB toolbox (bfmatlab) to read .nd2 files. This software can be downloaded from the [Open Microscopy Consortium](https://www.openmicroscopy.org/) [here](https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/). For auto-contrating DAPI images, Dentist2 uses the scale.m function available in the [rajlabimagetools repository](https://github.com/arjunrajlaboratory/rajlabimagetools). 

Before running Dentist2, please download bfmatlab and rajlabimagetools (or just scale.m) and add them to your path by typing the following in your MATLAB console. 

```matlab
>>addpath('~/path/to/dentist2/')
>>addpath('~/path/to/bfmatlab/')
>>addpath('~/path/to/rajlabimagetools/')
>>savepath
```
Of note, Dentist2 was written and tested using MATLAB version 2020b. We have not tested Dentist2 using versions prior to 2018b.  

Expected input
==============
Dentist2 expects images from a **tiled rectangular scan** in a **single z-plane**. There should be one DAPI channel and one or more FISH channels. The FISH channels may be named 'YFP', 'GFP', 'CY3', 'A594', 'CY5', 'A647' '700' or 'CY7'. If you'd like to process other fluorescence channels, you may need to modify the map in +d2utils/readND2Channels.m.  

Quick start
============
Change your working directory to the folder containing your scan file (`>>cd('~/path/to/scan/'`). If you have multiple scans you wish to analyze with Dentist2, we recommend moving each scan file into separate folders. If your scan needs to be stitched, launch the stitching GUI by typing the following into your console. 
```matlab
>>h = d2stitchingGUI(scanDimensions, 'scanFileName');
```
The scanDimensions should be formatted as \[number of rows, numbers of columns\]. Use the stitching GUI to select the scan layout and control points as described [below](#stitching). When you close the GUI window, two files will be written to your working directory: 'scanSummary.txt' and 'tilesTable.csv' (see description below).
**Annotated screenshot of stitching GUI**

<img src="https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/d2StitchingGUI.png" width="1100">

**Example scan patterns**

<img src="https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/stitchGUIeg.png" width="500">

Next, launch the thresholding GUI by typing the following:
```matlab
>>h = launchD2ThresholdGUI();
```
or if you're loading a pre-stitched scan:
```matlab
>>h = launchD2ThresholdGUI('preStitchedScan', 'path/to/scanFile.nd2');
```
When first launching the threshold GUI, it may take several minutes for the software to stitch the scan, segment and identify nuclei, and identify RNA FISH spots.  The processed data are automatically saved to your working directory and can be loaded more quickly if you need to close MATLAB and restart the GUI. In addition, if the GUI handle (h) is still in your workspace, you can relaunch the GUI by typing `>>h.relaunchGUI;` without having to reload the data.

As described [below](#d2ThresholdGUI), use the threshold GUI to adjust the spot intensity threshold and mask, add, or delete erroneous spots and cells. When you are finished, you can export the data into a summarized table of spots per cell ('spotsSummary.csv') by clicking on the export button. When you close the GUI window, data for all spot calls, nuclei and masks will be saved to 'spots.csv', 'nuclei.csv' and 'masks.csv', respectively.

**Annotated screenshot of threshold GUI**
<img src="https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/D2threshGUI.png" width="1200">

Stitching
==========

D2ThresholdGUI
==============

Main axes
---------
### Scatter
When starting the threshold GUI, the main axes will display a scatter plot of the spatial coordinates of all identified nuclei. Each nucleus (i.e. point) will be colored according to the number of spots assigned to it. You can change the color scheme by selecting a different colormap in the colormap drop-down menu &#9318;. The nuclei colors will adjust automatically as you change the FISH channel, modify the spot intensity threshold or create/delete masks. 

If your scans is larger than 20,000 pixels in any dimension, the main axes will open with a zoomed-in view of your scan. The position of this zoomed-in view will be indicated in the thumbnail axes. You may zoom-out by 2X by right-clicking (or control-clicking on a Mac) the main axes. Note that the more zoomed-out view may slow down the GUI.

You can toggle between the scatterplot and image overlay by selecting the scatter checkbox (5). However, if your current view is >64,000,000 pixels (e.g. >8,000 x 8,000), Dentist2 will not show the image overlay.  
### Image overlay
When toggling off the scatterplot, the main axes will show an overlay of DAPI and the current FISH channel cropped to your current view (note the exception above for very large views). You can use the sliders (4)  


Thumbnail axes
---------------

Threshold axes
--------------

Masks
-----

Importing CellPose masks
========================
If you use CellPose to segment nuclei, Dentist2 can use these outlines instead of it's default algorithm to identify cells. If you have a single CellPose outlines file (i.e. you ran cellpose on the pre-stitched scan), include the file name when running launchD2ThresholdGUI() as below:
```matlab
>>h = launchD2ThresholdGUI('cellPose', 'path/to/cp_outlines.txt');
```
Alternatively, if you ran CellPose on multiple image tiles that need to be stitched, specify the path to the folder containing the CellPose outlines as well as the name of a file that lists the position of each tile. For example: 
```matlab
>>h = launchD2ThresholdGUI('cellPose', 'path/to/cp_directory/', 'cellPoseTileTable', 'cellPoseTilePositions.csv');
```
An example of the cellPoseTilePositions.csv file can be found [here](). In addition, you may use +d2utils/splitStitchedScan.m to generate both the overlapping image tiles from a pre-stitched scan and the cellPoseTilePositions.csv file. For each outline in tile i, Dentist2 will determine if the outline overlaps any outlines from previous neighboring tiles, and if so, the overlapping outlines will be merged. This is to avoid duplicating nuclei that fall on tile boundaries. If your scan has a lot of nuclei, this step may take a very long time. You may be better off just running CellPose on a larger (stitched) image.  

Note that if you ran CellPose on a resized image (to save on memory), you will want to specify the resize factor when running launchD2ThresholdGUI() as below:
```matlab
>>h = launchD2ThresholdGUI('cellPose', 'path/to/cp_outlines.txt', 'maskResizeFactor', 4); %Indicates a resize factor of 4. 
```
or 
```matlab
>>h = launchD2ThresholdGUI('cellPose', 'path/to/cp_directory/', 'cellPoseTileTable', 'cellPoseTilePositions.csv', 'maskResizeFactor', 4);
```

Description of objects and default parameters
=============================================

Troubleshooting
===============




