function guiHandle = launchD2IFGUI(varargin)
    n = inputParser;
    n.addParameter('scanSummary', 'scanSummary.txt', @ischar);
    n.addParameter('imgMasksFile', 'imgMasks.csv', @ischar); 
    n.addParameter('cellMasksFile', 'cellMasks.csv', @ischar); 
    n.addParameter('nucBoundariesFile', 'nucBoundariesIF.csv', @ischar); 
    n.addParameter('cellBoundariesFile', 'cellBoundariesIF.csv', @ischar); 
    n.addParameter('IFquantFile', 'IFquantTable.csv', @ischar); 
    n.addParameter('preStitchedScan', '', @ischar);
    n.addParameter('cellPoseNuclei', '', @ischar);
    n.addParameter('cellPoseCyto', '', @ischar);
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
        scanObj.scanSummaryFile = n.Results.scanSummary;
        scanObj.saveScanSummary();
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.imgMasksFile)
        maskObjImg = maskTable(scanObj, n.Results.imgMasksFile);
    else
        fprintf('Unable to detect %s in your current directory. Creating a new masks object\n', n.Results.imgMasksFile)
        maskObjImg = maskTable(scanObj);
        maskObjImg.masksFile = n.Results.imgMasksFile;
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.cellMasksFile)
        maskObjBoundaries = maskTable(scanObj, n.Results.cellMasksFile);
    else
        fprintf('Unable to detect %s in your current directory. Creating a new masks object\n', n.Results.cellMasksFile)
        maskObjBoundaries = maskTable(scanObj);
        maskObjBoundaries.masksFile = n.Results.cellMasksFile;
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.nucBoundariesFile) && isfile(n.Results.cellBoundariesFile) && isfile(n.Results.IFquantFile) %%
        disp('Loading IFboundaries and IFquant objects.')
        IFboundariesObj = d2IF.IFboundaries(scanObj, maskObjBoundaries, n.Results.nucBoundariesFile, n.Results.cellBoundariesFile);
        IFquantObj = d2IF.IFtable(scanObj, maskObjImg, IFboundariesObj, n.Results.IFquantFile);
        IFboundariesObj.addColors();
        IFboundariesObj.addEmptyRows(1000);
        IFquantObj.addEmptyRows(1000);
    else
        disp('Making new IFboundaries and IFquant objects.')
        IFboundariesObj = d2IF.IFboundaries(scanObj, maskObjBoundaries);
        IFquantObj = d2IF.IFtable(scanObj, maskObjImg, IFboundariesObj);
        if isfile(n.Results.cellPoseNuclei) %Load pre-stitched cellpose mask
            disp('Loading cellpose nuclei boundaries.')
            IFboundariesObj.loadCellPoseDapi(sprintf('%s_masks.tif', n.Results.cellPose), sprintf('%s_outlines.txt',  n.Results.cellPose));
            IFboundariesObj.labelMat2nucTable();
        else
            disp('Calculating nuclei boundaries.')
            IFboundariesObj.makeNucleiLabelMat();
            IFboundariesObj.labelMat2nucTable();
        end
        
        if isfile(n.Results.cellPoseCyto)
        else
            disp('Quantifying nuclear and cytoplasmic signal.')
            IFquantObj.quantAllLabelMat2();
        end
    end
    IFquantObj.makeCentroidList('meanNuc');
%----------------------------------------------------------------
% 
    disp('Auto-contrasting stitched scans. This may take several minutes.')
    scanObj.contrastDAPIstitch();
    scanObj.contrastStitchedScans([1 99], [0.9 3]);
    disp('Resizing stitched scans')
    scanObj.resizeStitchedScans();
    
    guiHandle = d2IF.d2IFView(scanObj, IFboundariesObj, IFquantObj);
end