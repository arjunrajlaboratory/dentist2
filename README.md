Dentist2
========
This repository contains tools to help analyze RNA FISH data from large image scans. Example input and output data are available on Dropbox here.  

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
Dentist2 uses the Bio-Formats MATLAB toolbox (bfmatlab) to read .nd2 files. This software can be downloaded from the [Open Microscopy Consortium](https://www.openmicroscopy.org/)[here](https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/). For auto-contrating DAPI images, Dentist2 uses the scale.m function available in the [rajlabimagetools repository](https://github.com/arjunrajlaboratory/rajlabimagetools). 

Before trying to run Dentist2, please download bfmatlab and rajlabimagetools (or just scale.m) and add them to your path by typing the following in your MATLAB console. 

```matlab
test
```
Of note, Dentist2 was written and tested using MATLAB version 2020b. We have not tested Dentists using versions prior to 2018b.  

Expected input
==============


Quick start
============

Stitching
==========
