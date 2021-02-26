Dentist2
========
This repository contains tools to help analyze RNA FISH data from large image scans. Example input and output data are available on Dropbox [here].  

Table of contents
=================
* [Dentist2](#dentist2)
* [Dependencies](#dependencies)
* [Expected input](#expected-input)
* [Quick start](#quick-start)
* [Stitching](#stitching)
* [D2ThresholdGUI](#d2ThresholdGUI)
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
Of note, Dentist2 was written and tested using MATLAB version 2020b. We have not tested Dentists using versions prior to 2018b.  

Expected input
==============
Dentist2 expects images from a **tiled rectangular scan** in a **single z-plane**. There should be one DAPI channel and one or more FISH channels. The FISH channels may be named 'YFP', 'GFP', 'CY3', 'A594', 'CY5', 'A647' '700' or 'CY7'. If you'd like to process other fluorescence channels, you may need to modify the map in +d2utils/readND2Channels.m.  

Quick start
============
Change your working directory to the folder containing your scan file (`>>cd('~/path/to/scan/'`). If you have multiple scans you wish to analyze with Dentist2, we recommend moving each scan file into separate folders. Launch the stitching GUI by typing the following into your console. 
```matlab
>>h = d2stitchingGUI(scanDimensions, 'scanFileName');
```
The scanDimensions should be formatted as \[number of rows, numbers of columns\]. Use the stitching GUI to select the scan layout and control points as described [below](#stitching). When you close the GUI window, two files will be written to your working directory: 'scanSummary.txt' and 'tilesTable.csv' (see description below).

Next, launch the thresholding GUI by typing the following:
```matlab
>>h = launchD2ThresholdGUI();
```
When first launching the threshold GUI, it may take several minutes for the software to stitch the scan, segment and identify nuclei, and identify RNA FISH spots.  The processed data are automatically saved to your working directory and can be reloaded more quickly if you need to close MATLAB and restart the GUI. In addition, if the GUI handle (h) is still in your workspace, you can relaunch the GUI by typing `>>h.relaunchGUI;` without having to reload the data.

As described [below](#d2ThresholdGUI), use the threshold GUI to adjust the spot intensity threshold and mask, add, or delete erroneous spots and cells. When you are finished, you can export the data into a summarized table of spots per cell ('spotsSummary.csv') by clicking on the export button. When you close the GUI window, data for all spot calls, nuclei and masks will be saved to 'spots.csv', 'nuclei.csv' and 'masks.csv', respectively.

Stitching
==========

D2ThresholdGUI
==============

Description of objects and default parameters
=============================================

Troubleshooting
===============




