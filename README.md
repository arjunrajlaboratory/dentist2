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
Dentist2 uses the Bio-Formats MATLAB toolbox (bfmatlab) to read .nd2 files. This software can be downloaded from the [Open Microscopy Consortium](https://www.openmicroscopy.org/) [here](https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/). Dentist2 also uses a few functions from the [rajlabimagetools repository](https://github.com/arjunrajlaboratory/rajlabimagetools). 

Before running Dentist2, please download bfmatlab and rajlabimagetools and add them to your path by typing the following in your MATLAB console. 

```matlab
>>addpath('~/path/to/dentist2/')
>>addpath('~/path/to/bfmatlab/')
>>addpath('~/path/to/rajlabimagetools/')
>>savepath
```
Of note, Dentist2 was written and tested using MATLAB version 2020b. We have not tested Dentist2 using versions prior to 2018b. In addition, Dentist2 uses functions from a variety of built-in MATLAB toolboxes. We recommend installing all available toolboxes when downloading/updating MATLAB. 
 
Expected input
==============
Dentist2 expects images from a **tiled rectangular scan** in a **single z-plane**. There should be one DAPI channel and one or more FISH channels. The FISH channels may be named 'YFP', 'GFP', 'CY3', 'A594', 'CY5', 'A647' '700' or 'CY7'. If you'd like to process other fluorescence channels, you may need to modify the map in +d2utils/readND2Channels.m. 

At the moment, Dentist2 identifies nuclei based on DAPI signal and assigned RNA FISH spots to the nearest nucleus. In the future, we may add the option to input whole-cell segmentations and assign spots to the surrounding cell or, for IF, average fluorescence intensity per cell.

Quick start
============
Change your working directory to the folder containing your scan file (`>>cd('~/path/to/scan/'`). If you have multiple scans you wish to analyze with Dentist2, we recommend moving each scan file into separate folders. If your scan needs to be stitched, launch the stitching GUI by typing the following into your console. 
```matlab
>>h = d2stitchingGUI(scanDimensions, 'scanFileName');
```
The scanDimensions should be formatted as \[number of rows, numbers of columns\]. Use the stitching GUI to select the scan layout and control points as described [below](#stitching). When you close the GUI window, two files will be written to your working directory: 'scanSummary.txt' and 'tilesTable.csv'.

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

![thresholdGUI](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/D2threshGUI.png)

Note these hotkeys:
- 'up arrow': move up centroid list
- 'down arrow': move down centroid list
- 'z': zoom-mode
- 'p': pan-view mode
- 's': overlay spots checkbox
- 'n': overlay nuclei checkbox
- 'c': scatter checkbox
- 'm': add spot mask
- 'shift-m': add cell mask
- 'd': delete mask
- 'shift-s': save tables
- 'shift-e': export spots summary

Stitching
==========
Use the scan pattern checkboxes &#9313; to select the pattern corresponding to how your scan images were acquired. If unsure, use the show new positions button &#9314; to display a random set of tiles consistent with the selected pattern. Note that some patterns may place the same tiles at the same positions in a scan, and it can be helpful to check several scan positions. 

Once you've chosen your scan pattern, use the show new positions button &#9314; or textbox underneath to select a position with plenty of nuclei on the border of the images. Click select row control points &#9315; to launch MATLAB's [cpselect](https://www.mathworks.com/help/images/ref/cpselect.html) tool. This should create a new figure window. Select matching pixels on the top and bottom borders of the images like so: 

![cpSelect Row](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/cpSelectRow.png)

Close the window when finished. Click select column control points &#9315; and select matching pixels on the right and left borders of the images like so: 

![cpSelect Col](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/cpSelectColumn.png)

Close the window when finished.

You may click on preview stitch &#9316; to view a stitched image based on the stored control points. If you are satisfied with the stitch, you may exit out of the GUI. Otherwise, you can change scan positions and select additional control points. 

When you exit out of the GUI, the files 'tilesTable.csv' and 'scanSummary.txt' will be saved to your working directory. TilesTable.csv contains the relative coordinates of all tiles in your scan.

D2ThresholdGUI
==============

Main axes
---------
### Scatter
When starting the threshold GUI, the main axes will display a scatter plot of the spatial coordinates of all identified nuclei. Each nucleus (i.e. point) will be colored according to the number of spots assigned to it. You can change the color scheme by selecting a different colormap in the colormap drop-down menu &#9318;. The nuclei colors will adjust automatically as you change the FISH channel, modify the spot intensity threshold or create/delete masks. 

If your scans is larger than 20,000 pixels in any dimension, the main axes will open with a zoomed-in view of your scan. The position of this zoomed-in view will be indicated in the thumbnail axes. You may zoom-out by 2X by right-clicking (or control-clicking on a Mac) the main axes. Note that the more zoomed-out view may slow down the GUI.

You can toggle between the scatterplot and image overlay by selecting the scatter checkbox &#9316;. However, if your current view is >64,000,000 pixels (e.g. >8,000 x 8,000), Dentist2 will not show the image overlay.  
### Image overlay
When toggling off the scatterplot, the main axes will show an overlay of DAPI and the current FISH channel cropped to your current view (note the exception above for very large views). You can use the sliders &#9315;.   


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
* Problem: "Error using X (line #). Undefined function 'A' for input arguments of type 'B'." 
  * Possible Reason: Some MATLAB toolkit is missing. Your MATLAB installation is trying to use a different function named 'A' than what was used for writing Dentist2. 
  * Possible Solution: If you can figure out which toolkit is missing, you can install that toolkit. Alternatively, you can try reinstalling MATLAB with all available toolkits.
* Problem: "Java exception occured: java.lang.OutOfMemoryError: Java heap space". 
  * Possible Reason: You are trying to read a very large image or scan into MATLAB using bfmatlab functions and your Java heap size is set too low. 
  * Possible Solution: Increase your Java heap size by selecting Preferences > MATLAB > General > Java Heap Memory. Instructions available ![here](https://www.mathworks.com/help/matlab/matlab_external/java-heap-memory-preferences.html). 
* Problem: Dentist2 is missing FISH spots that are clearly present in the image. Oddly, this is happening in blocks like so: 
![missingSpotsExample1](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/missingSpotsExample1.png).
  * Possible Reason: When finding spots, dentist2 automatically sets a minimum intensity threshold in order to limit the size of the spotsTable (see the function d2utils.findSpotsaTrous.m). The default is max(spotIntensity)/4. If your FISH signal intensity is low and there is bright schmutz in your images, then your FISH signal may fall below the minimum intensity threshold. Also note that denists2 finds spots in blocks in parallel to speed up processing. Thus different blocks may have different max(spotIntensity) values. 
  * Possible Solution: You can decrease the minimum intensity threshold for a specific channel using the method changeSpotFilterThreshold in the d2MainAxesController. For example, by typing `h.mainAxesCntrlr.changeSpotFilterThreshold('cy', 10);`  into your MATLAB console, you can set the minimum intensity threshold for 'cy' to max(spotIntensity)/10 (assuming h is your GUI object). Alternatively, you can change the minimum intensity threshold using the 'aTrousMinThreshFactor' parameter in the launchD2ThresholdGUI.m function, for example: `h = launchD2ThresholdGUI('aTrousMinThreshFactor', 10)`. 
    * Note that if you close the GUI and rerun `h = launchD2ThresholdGUI()` it will set the minimum intensity threshold to it's default max(spotIntensity)/4. If available, `h = launchD2ThresholdGUI()` will load the spots.csv table even if those spots were found using a different intensity threshold. However, if you try to refind spots (because spots.csv is missing or because you ran subtractBackground as described below), dentist2 will use the default threshold (or whatever value is stored in h.spotTable.threshFactor).

* Problem: Dentist2 is calling extra spots. Oddly, these spots are falling in straight lines in what appears to be the boundary between tiles like so:
![]()
  * Possible Reason:
  * Possible Solution:

