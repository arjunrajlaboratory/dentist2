classdef IFtable < handle
    
    properties (Access = public)
        
        scanObj
        maskObj
        IFquant
        IFquant2
        nucBoundaries %based on cellpose outlines
        nucBoundaries2 %based on labelmat
        cellBoundaries %based on cellpose outlines
        cellBoundaries2 %based on labelmat
        channels
        dapiMask
        dapiLabelMat
        dapiCC %Might not need to be stored as property
        dapiRP
        minNucleusSize = 1000;
        radius = 20;
        IFfile
    end
    
    methods
        function p = IFtable(scanObject, maskObj, varargin)
            p.scanObj = scanObject;
            p.maskObj = maskObj;
            p.channels = p.scanObj.stitchedScans.labels;
            if nargin == 2
                fprintf('New Table\n');
                p.IFquant = table('size', [0,numel(p.channels)+ 9],... %Possibly unnecessary since new table created with p.findNuclei. 
                    'VariableNames', [{'cellID', 'x', 'y', 'status'}, p.channels, {'meanNuc', 'meanCyto',  'sumNuc', 'sumCyto','maskID'}],...
                    'VariableTypes', [repmat({'single'}, 1, 3), repmat({'logical'}, 1, 1+numel(p.channels)), repmat({'single'}, 1, 5)]);
            elseif nargin == 3
                fprintf('Loading Table\n');
                p.IFfile = varargin{1};
                opts = detectImportOptions(varargin{1});
                opts = setvartype(opts, 'single');
                p.IFquant = readtable(varargin{1},opts);
                p.IFquant(:,[{'status'},  p.channels]) = logical(p.cells(:,[{'status'},  p.channels])); %For some reason, when I set 'status' to 'logical' they all go to false. So doing this instead
            end
        end
        
        function p = stitchDAPImask(p, varargin)  
            tileTable = p.scanObj.tilesTable;
            tilesTmp = transpose(p.scanObj.scanMatrix);
            tiles = tilesTmp(:);
            height = p.scanObj.tileSize(1);
            width = p.scanObj.tileSize(2);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'logical');
            channel = find(ismember(p.scanObj.channels,'dapi'));
            reader = bfGetReader(p.scanObj.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            
            if nargin == 1
                s = 0.1;
            elseif nargin == 2
                s = varargin{1};
            end
                
            for i = 1:numel(tiles)
                
                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);
                
                tileMask = d2utils.makeDAPImask(scale(tmpPlane), 'sensitivity', s);
                
                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tileMask;
            end
            reader.close()
            tmpCC = bwconncomp(tmpStitch);
            tmpArea = regionprops(tmpCC, 'area');
            schmutzIdx = [tmpArea.Area] < p.minNucleusSize;
            schmutz = tmpCC.PixelIdxList(schmutzIdx);
            for i = 1:numel(schmutz)
                tmpStitch(schmutz{i}) = false;
            end
            p.dapiMask = tmpStitch;
        end
        
        function p = stitchDAPImask2(p, varargin) %For pre-stitched scans
            
            n = inputParser;
            n.addParameter('sensitivity', 0.1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',0,'<=',1.0}));    
            n.addParameter('blockSize', [1500 1500], @(x)validateattributes(x,{'numeric'}, {'size', [1 2]}));    
            n.parse(varargin{:});
            %Should maybe check that the block size is >[1 1] and < scanDim.
            s = n.Results.sensitivity;
            block = n.Results.blockSize;
            
            function_mask = @(block_struct) imbinarize(scale(block_struct.data),...
                adaptthresh(scale(block_struct.data), s, 'ForegroundPolarity','bright'));
            
            tmpStitch = blockproc(p.scanObj.dapiStitch, block, function_mask, 'BorderSize', [0 0], 'UseParallel', true);
            tmpCC = bwconncomp(tmpStitch);
            tmpArea = regionprops(tmpCC, 'area');
            schmutzIdx = [tmpArea.Area] < p.minNucleusSize;
            schmutz = tmpCC.PixelIdxList(schmutzIdx);
            for i = 1:numel(schmutz)
                tmpStitch(schmutz{i}) = false;
            end
            p.dapiMask = tmpStitch;
        end
        
        function p = makeNucleiLabelMat(p)
            %First rounds out dapi masks then finds CC and makes
            %label matrix
            tmpBoundaries = bwboundaries(p.dapiMask);
            tmpBoundariesArea = cellfun(@(x) polyarea(x(:,1), x(:,2)), tmpBoundaries);
            tmpBoundaries = tmpBoundaries(tmpBoundariesArea>= 1000); 
            
            boundaryStitch = zeros(size(p.dapiMask)); %Could make this max of boundaryArray. Or stitchDim
           
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            for i = 1:numel(tmpBoundaries)
                [tmpXlim,tmpYlim]  = boundingbox(polyshape(tmpBoundaries{i}));
                tmpMask = poly2mask(tmpBoundaries{i}(:,2)-tmpYlim(1), tmpBoundaries{i}(:,1)-tmpXlim(1), diff(tmpXlim)+1, diff(tmpYlim)+1);
                boundaryStitch(tmpXlim(1):tmpXlim(2), tmpYlim(1):tmpYlim(2)) = imclose(tmpMask, strel('disk', 30)); %Consider adjusting operation or strel. 
            end
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
            p.dapiCC = bwconncomp(boundaryStitch);
            p.dapiRP = regionprops(p.dapiCC);
            p.dapiLabelMat = labelmatrix(p.dapiCC);
        end
        
        function p = quantAllChannels(p, varargin)
        %Uses use polyshapes to mask nuclei and cells.
            if nargin == 1
                channelsToQaunt = p.channels;
            elseif nargin ==2
                channelsToQaunt = varargin{1};
            end
            
%             cellIDs = unique(p.dapiLabelMat);
%             cellIDs(cellIDs == 0) = [];
            nRows = prod(numel(p.dapiRP), numel(p.channels));
            meanNuc = zeros(nRows,1);
            meanCyto = zeros(nRows,1);
            sumNuc = zeros(nRows,1);
            sumCyto = zeros(nRows,1);
            channelStatus = false(nRows,numel(p.channels));
            nucBoundariesTmp = cell(0, numel(p.dapiRP));
            cellBoundariesTmp = cell(0, numel(p.dapiRP));
            cellCoords = zeros(nRows,2);
            for i = 1:numel(p.dapiRP) %can make this parfor?
%                 tmpRP = regionprops(p.dapiLabelMat == cellIDs(i));
%                 if size(tmpRP, 1) == 1 %Allows for disconnected masks which can occur with cellpose. 
%                     tmpBB = tmpRP.BoundingBox;
%                     tmpCentroid = tmpRP.Centroid;
%                 else
%                     tmpBB2 = vertcat(tmpRP.BoundingBox);
%                     tmpBB = [min(tmpBB2(:,1:2)), max(tmpBB2(:,3:4))];
%                     tmpCentroid = mean(vertcat(tmpRP.Centroid));
%                 end
                tmpBB = d2utils.rotateRectROI(p.dapiRP(i).BoundingBox);
                tmpSize = tmpBB(3:4) + (2*p.radius);
                tmpStart = max([1,1], tmpBB(1:2)-p.radius-5); %Add a buffer beyond the p.radius
                tmpEnd = min(tmpStart+tmpSize+10, size(p.dapiLabelMat)); %Could make this stitchDim
                
                tmpRegionMask = p.dapiLabelMat(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
                tmpDapiMask = tmpRegionMask == i;
                tmpDapiMaskDilated = imdilate(tmpDapiMask, strel('disk', p.radius)); %Consider other strel. Including tmpDapiMask. 
                tmpCytoMask = tmpDapiMaskDilated & ~logical(tmpRegionMask);
                tmpNucBoundary = bwboundaries(tmpDapiMask, 'noholes');
                tmpCellBoundary = bwboundaries(tmpDapiMaskDilated, 'noholes');
                nucBoundariesTmp{i} = cellfun(@(x) x+tmpStart, tmpNucBoundary, 'UniformOutput', false);
                cellBoundariesTmp{i} = cellfun(@(x) x+tmpStart, tmpCellBoundary, 'UniformOutput', false);
                tmpRP = regionprops(tmpDapiMask);
                if size(tmpRP, 1) == 1 %Allows for disconnected masks which can occur with cellpose. 
                    tmpCentroid = fliplr(round(tmpRP.Centroid))+tmpStart;
                    cellCoords(i,:) = single(round(tmpCentroid));
                elseif size(tmpRP, 1) > 1
                    tmpCentroid = fliplr(mean(vertcat(tmpRP.Centroid)))+tmpStart; %For weighted centroid, can convert to polyshape then calc centroid 
                    cellCoords(i,:) = single(round(tmpCentroid));
                end
                for ii = 1:numel(channelsToQaunt)
                    stitchIdx = ismember(p.channels, channelsToQaunt(ii)); %Redundant in some cases but flexible
                    tmpImage = p.scanObj.stitchedScans.stitches{stitchIdx}(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
%                     tmpImage = im2single(tmpImage);
                    tmpIdx = (i-1)*numel(channelsToQaunt)+ii;
                    meanNuc(tmpIdx) = single(mean(tmpImage(tmpDapiMask)));
                    sumNuc(tmpIdx) = sum(single(tmpImage(tmpDapiMask)));
                    meanCyto(tmpIdx) = single(mean(tmpImage(tmpCytoMask)));
                    sumCyto(tmpIdx) = sum(single(tmpImage(tmpCytoMask)));
                    channelStatus(tmpIdx,:) = stitchIdx;
                end
            end
            status = true(nRows,1);
            maskID = single(zeros(nRows,1));

            p.IFquant = array2table([(1:nRows)',cellCoords, status, channelStatus, meanNuc, meanCyto, sumNuc, sumCyto, maskID],...
                'VariableNames', [{'cellID', 'x', 'y', 'status'}, p.channels, {'meanNuc', 'meanCyto',  'sumNuc', 'sumCyto','maskID'}]);
            p.IFquant = convertvars(p.IFquant, [{'status'}, p.channels], 'logical');
            
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            nucBoundariesArray = cellfun(@(x) polyshape(vertcat(x{:})), nucBoundariesTmp, 'UniformOutput', false); %Vertcat in case there are disconnected regions
            nucBoundariesArray = cellfun(@(x) single(x.Vertices), nucBoundariesArray, 'UniformOutput', false);
            nucBoundariesHeight = cellfun(@(x) height(x), nucBoundariesArray, 'UniformOutput', true);
            nucIDArray = repelem(1:numel(nucBoundariesHeight), nucBoundariesHeight);
            p.nucBoundaries2 = array2table([single(nucIDArray'), vertcat(nucBoundariesArray{:})], 'VariableNames', {'cellID', 'x', 'y'});
            
            cellBoundariesArray = cellfun(@(x) polyshape(vertcat(x{:})), cellBoundariesTmp, 'UniformOutput', false);
            cellBoundariesArray = cellfun(@(x) x.Vertices, cellBoundariesArray, 'UniformOutput', false);
            cellBoundariesHeight = cellfun(@(x) height(x), cellBoundariesArray, 'UniformOutput', true);
            cellIDArray = repelem(1:numel(cellBoundariesHeight), cellBoundariesHeight);
            p.cellBoundaries2 = array2table([single(cellIDArray'), vertcat(cellBoundariesArray{:})], 'VariableNames', {'cellID', 'x', 'y'});
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function p = loadCellPoseDapi(p, labelMatFile, outlineFile)
            tmpLabelMat = imread(labelMatFile);
            if all(size(tmpLabelMat) ==  p.scanObj.stitchDim)
                p.dapiLabelMat =  tmpLabelMat;
            else
                p.dapiLabelMat = imresize(tmpLabelMat, p.scanObj.stitchDim, 'nearest');
                scaleFactor = p.scanObj.stitchDim./size(tmpLabelMat);
            end
            
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            warning('off', 'MATLAB:polyshape:boundary3Points')
            polyVect = d2utils.parseCellposeOutlines(outlineFile, 'scaleFactor', fliplr(scaleFactor));
            polyVect = scale(polyVect, fliplr(scaleFactor));
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
            warning('on', 'MATLAB:polyshape:boundary3Points')
            tmpArea = num2cell([polyVect.area]);
            [tmpX, tmpY] = centroid(polyVect);
            tmpCentroid = num2cell([tmpX', tmpY'],2);
            nucBoundariesArray = num2cell(polyVect);
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
%             tmpArea = cellfun(@(x) polyarea(x(:,1), x(:,2))*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpCentroid = cellfun(@(x) d2utils.poly2centroid(x)*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpBB = cellfun(@(x) d2utils.polygonBoundingBox2(x)*scaleFactor,polyArray, 'UniformOutput', false);
            p.dapiRP = cell2struct([tmpArea', tmpCentroid, tmpBB'], {'Area', 'Centroid', 'BoundingBox'}, 2);
            
            nucBoundariesArray = cellfun(@(x) single(x.Vertices), nucBoundariesArray, 'UniformOutput', false);
            nucBoundariesHeight = cellfun(@(x) height(x), nucBoundariesArray, 'UniformOutput', true);
            nucIDArray = repelem(1:numel(nucBoundariesHeight), nucBoundariesHeight);
            p.nucBoundaries = array2table([single(nucIDArray'), fliplr(vertcat(nucBoundariesArray{:}))], 'VariableNames', {'cellID', 'x', 'y'});
            
            polyVectDilated = polybuffer(polyVect, p.radius);
            cellBoundariesArray = num2cell(polyVectDilated);
            cellBoundariesArray = cellfun(@(x) x.Vertices, cellBoundariesArray, 'UniformOutput', false);
            cellBoundariesHeight = cellfun(@(x) height(x), cellBoundariesArray, 'UniformOutput', true);
            cellIDArray = repelem(1:numel(cellBoundariesHeight), cellBoundariesHeight);
            p.cellBoundaries = array2table([single(cellIDArray'), fliplr(vertcat(cellBoundariesArray{:}))], 'VariableNames', {'cellID', 'x', 'y'});
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function outCells = getCellsInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.IFquant.x >= rect(1) & p.IFquant.x < rect(1) + rect(3) ...
                & p.IFquant.y >= rect(2) & p.IFquant.y < rect(2) + rect(4) ...
                & p.IFquant.status;
            
            outCells = p.IFquant(idx,:);
        end
        
    end
end
