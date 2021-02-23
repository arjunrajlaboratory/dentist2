function guiHandle = launchD2ThresholdGUI(varargin)
    n = inputParser;
    n.addParameter('scanSummary', 'scanSummary.txt', @ischar);
    n.addParameter('masksFile', 'masks.csv', @ischar); 
    n.addParameter('nucleiFile', 'nuclei.csv', @ischar); 
    n.addParameter('spotsFile', 'spots.csv', @ischar); 
    
    n.parse(varargin{:});
%----------------------------------------------------------------
%
    if isfile(n.Results.scanSummary)
        scanObj = scanObject('scanSummary', n.Results.scanSummary);
    else
        fprintf('Unable to detect %s in your current directory.\n. Make sure to run the d2StitchingGUI before launching the d2ThresholdGUI.\n You may also want to check your path and %s and try again. ', n.Results.scanSummary, n.Results.scanSummary)
        return
    end
    
    %Check if stitches have been saved. If not, stitch and save to default
    %files. 
    if isempty(scanObj.tilesTable)
        disp('The scan object does not contain a tiles table. Creating a new tiles table.')
        scanObj.loadTiles();
        scanObj.savetilesTable();
    end
    
    if isfile('stitchedScans.mat')
        fprintf('Loading stitched scans.\nThis may take several minutes.\n')
        scanObj.loadStitches();
    else
        disp('Stitching DAPI channel. This may take a few minutes.')
        scanObj.stitchDAPI();
        disp('Stitching FISH channels. This may take a few minutes.')
        scanObj.stitchChannels();
        disp('Saving stitched scans. This may take several minutes.')
        scanObj.saveStitches();
    end
    
    disp('Auto-contrasting stitched scans. This may take several minutes.')
    scanObj.contrastDAPIstitch();
    scanObj.contrastStitchedScans([1 99], [0.9 3]);
    disp('Resizing stitched scans')
    scanObj.resizeStitchedScans();
%----------------------------------------------------------------
%   
    if isfile(n.Results.masksFile)
        maskObj = maskTable(scanObj, n.Results.masksFile);
    else
        fprintf('Unable to detect %s in your current directory. Creating a new masks object\n', n.Results.masksFile)
        maskObj = maskTable(scanObj);
    end
%----------------------------------------------------------------
%   
    if isfile(n.Results.nucleiFile)
        disp('Loading nuclei file.')
        nucleiObj = nucleiTable(scanObj, maskObj, n.Results.nucleiFile);
        nucleiObj.addColors();
    else
        fprintf('Unable to detect %s in your current directory. Creating a new nuclei object\n', n.Results.nucleiFile)
        nucleiObj = nucleiTable(scanObj, maskObj);
        disp('Finding nuclei. This may take a few minutes.')
        nucleiObj.stitchDAPImask();
        nucleiObj.findNuclei();
        nucleiObj.addColors();
        disp('Saving nuclei file.')
        nucleiObj.saveNucleiTable();
    end
    
%----------------------------------------------------------------
% 
    if isfile(n.Results.spotsFile)
        spotsObj = spotTable(scanObj, maskObj, nucleiObj, n.Results.scanSummary, n.Results.spotsFile);
    else
        fprintf('Unable to find %s in your current directory. Creating a new spots object\n', n.Results.spotsFile)
        spotsObj = spotTable(scanObj, maskObj, nucleiObj, n.Results.scanSummary);
        disp('Finding spots. This may take several minutes.')
        spotsObj.findSpots();
        disp('Finished finding spots')
        spotsObj.assignSpotsToNuclei();
    end
    
    if isempty(spotsObj.thresholds)
        spotsObj.defaultThresholds();
    end

    spotsObj.updateScanSummary();
    spotsObj.updateAllSpotStatus();
    spotsObj.makeCentroidList();
    
    guiHandle = d2ThresholdView2(scanObj, maskObj, nucleiObj, spotsObj);

end