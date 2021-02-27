function guiHandle = launchD2ThresholdGUI(varargin)
    n = inputParser;
    n.addParameter('scanSummary', 'scanSummary.txt', @ischar);
    n.addParameter('masksFile', 'masks.csv', @ischar); 
    n.addParameter('nucleiFile', 'nuclei.csv', @ischar); 
    n.addParameter('spotsFile', 'spots.csv', @ischar); 
    n.addParameter('preStitchedScan', '', @ischar);
    n.addParameter('cellPose', '', @ischar);
    n.addParameter('maskResizeFactor', 1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>', 0}));
    n.addParameter('cellPoseTileTable', 'cellPoseTilePositions.csv', @ischar);
    
    n.parse(varargin{:});
%----------------------------------------------------------------
%
    if isempty(n.Results.preStitchedScan)
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
    else
        fprintf('Loading pre-stitched scans.\nThis may take several minutes.\n')
        scanObj = scanObject('scanFile', n.Results.preStitchedScan);
        scanObj.saveScanSummary();
    end
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
    else
         fprintf('Unable to detect %s in your current directory. Creating a new nuclei object\n', n.Results.nucleiFile)
         nucleiObj = nucleiTable(scanObj, maskObj);
         if isempty(n.Results.cellPose) %No cellpose
            disp('Finding nuclei. This may take a few minutes.')
            if isempty(n.Results.preStitchedScan)
                nucleiObj.stitchDAPImask();
            else
                nucleiObj.stitchDAPImask2();
            end
            nucleiObj.findNuclei();
%             disp('Saving nuclei file.') Table will be automatically saved when closing GUI.
%             nucleiObj.saveNucleiTable(); 
        elseif isfile(n.Results.cellPose) %Load pre-stitched cellpose mask
            fprintf('Loading pre-stitched cellpose mask:%s.\nResize factor is %d.\n', n.Results.cellPose, n.Results.maskResizeFactor)
            nucleiObj.loadCellPoseMasks(n.Results.cellPose, n.Results.maskResizeFactor);
%           disp('Saving nuclei file.') 
%           nucleiObj.saveNucleiTable();
        elseif isfolder(n.Results.cellPose) %Load & stitch cellpose masks
            if isfile(n.Results.cellPoseTileTable)
                fprintf('Stitching cellpose masks in directory:%s.\nResize factor is %d.\n', n.Results.cellPose, n.Results.maskResizeFactor)
                nucleiObj.stitchCellPoseMasks(n.Results.cellPoseTileTable, n.Results.cellPose, n.Results.maskResizeFactor);
%               nucleiObj.saveNucleiTable();
            else
                fprintf('Unable to find cellpose file table: %s.\nThe file is need for stitching maks in %s\n', n.Results.cellPoseTileTable, n.Results.cellPose)
                disp('If you intended to input a pre-stitched mask, please specify the full filename rather than the name of a directory.')
            end
         else
            fprintf('Unable to find file or folder %s', n.Results.cellPose)
            return
         end
    end
    nucleiObj.addColors();
    nucleiObj.updateAllMasks();
%----------------------------------------------------------------
%   
    if isfile(n.Results.spotsFile)
        spotsObj = spotTable(scanObj, maskObj, nucleiObj, n.Results.scanSummary, n.Results.spotsFile);
    else
        fprintf('Unable to find %s in your current directory. Creating a new spots object\n', n.Results.spotsFile)
        spotsObj = spotTable(scanObj, maskObj, nucleiObj, n.Results.scanSummary);
        disp('Finding spots. This may take several minutes.')
        spotsObj.findSpots4(); %Only run findSpots4 on non-contrasted stitched scans!
        spotsObj.maskBorderSpots();
        disp('Finished finding spots')
        spotsObj.assignSpotsToNuclei();
    end
    
    if isempty(spotsObj.thresholds)
        spotsObj.defaultThresholds();
    end

    spotsObj.updateScanSummary();
    spotsObj.updateAllMasks();
    spotsObj.updateAllSpotStatus();
    spotsObj.makeCentroidList();

%----------------------------------------------------------------
% 
    disp('Auto-contrasting stitched scans. This may take several minutes.')
    scanObj.contrastDAPIstitch();
    scanObj.contrastStitchedScans([1 99], [0.9 3]);
    disp('Resizing stitched scans')
    scanObj.resizeStitchedScans();
    
    guiHandle = d2ThresholdView2(scanObj, maskObj, nucleiObj, spotsObj);
end