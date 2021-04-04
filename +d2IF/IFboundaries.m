classdef IFboundaries < handle
    
    properties (Access = public)
       scanObj 
       maskObj %Masks for cells and nuclei. Not for IF signal. 
       nucBoundaries %based on cellpose outlines
       nucBoundaries2 %based on labelmat
       cellBoundaries %based on cellpose outlines
       cellBoundaries2 %based on labelmat
       channels
       dapiMask
       dapiLabelMat
       dapiRP
       randomColors
       nucBoundariesFile = 'nucBoundariesIF.csv'
       cellBoundariesFile = 'cellBoundariesIF.csv'
       minNucleusSize = 1000;
       radius = 20;
    end
    
    methods
        function p = IFboundaries(scanObject, maskObj, varargin)
            n = inputParser;
            n.addRequired('scanObject', @(x) isa(x, 'scanObject')) 
            n.addRequired('maskObj', @(x) isa(x, 'maskTable'))
            n.addOptional('nucleiFile', '', @ischar)
            n.addOptional('cellsFile', '', @ischar)
            n.parse(scanObject, maskObj, varargin{:})
            
            p.scanObj = n.Results.scanObject;
            p.maskObj = n.Results.maskObj;
            p.channels = p.scanObj.stitchedScans.labels;
            p.randomColors = single(d2utils.distinguishable_colors(50));
            if isfile(n.Results.nucleiFile)
                p.nucBoundariesFile = n.Results.nucleiFile;
                tmpNucBoundaries = readtable(p.nucBoundariesFile);
                tmpNucBoundaries = rowfun(@(x) polyshape(x), tmpNucBoundaries, 'InputVariables', {'x', 'y'}, 'SeparateInputs', false, 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'nucBoundary');
                tmpNucBoundaries(:,'GroupCount') = [];
                tmpNucBoundaries.nucBB = cell2mat(arrayfun(@(x) d2utils.polyshapeBoundingBox(x), tmpNucBoundaries.nucBoundary, 'UniformOutput', false));
                tmpNucBoundaries = convertvars(tmpNucBoundaries, p.channels,'logical');
                p.nucBoundaries2 = movevars(tmpNucBoundaries, {'nucBoundary', 'nucBB'}, 'After', 'cellID');
            end
            
            if isfile(n.Results.cellsFile)
                p.cellBoundariesFile = n.Results.cellsFile;
                tmpCellBoundaries = readtable(p.cellBoundariesFile);
                tmpCellBoundaries = rowfun(@(x) polyshape(x), tmpCellBoundaries, 'InputVariables', {'x', 'y'}, 'SeparateInputs', false, 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'cellBoundary');
                tmpCellBoundaries(:,'GroupCount') = [];
                tmpCellBoundaries.cellBB = cell2mat(arrayfun(@(x) d2utils.polyshapeBoundingBox(x), tmpCellBoundaries.cellBoundary, 'UniformOutput', false));
                tmpCellBoundaries = convertvars(tmpCellBoundaries, p.channels,'logical');
                p.cellBoundaries2 = movevars(tmpCellBoundaries, {'cellBoundary', 'cellBB'}, 'After', 'cellID');
            end
        end
        
        function p = addEmptyRows(p, n)
            newNuclei = table('Size', [n, width(p.nucBoundaries2)],...
                'VariableNames', p.nucBoundaries2.Properties.VariableNames,...
                'VariableTypes', varfun(@class, p.nucBoundaries2, 'output','cell'));
            newNuclei{:,p.channels} = false(n, numel(p.channels));
            newNuclei.nucBB = zeros(n, 4);
            if ismember('colors', p.nucBoundaries2.Properties.VariableNames)
                newNuclei.colors = zeros(n, 3);
            end
            p.nucBoundaries2 = [p.nucBoundaries2; newNuclei];
            
            newCells = table('Size', [n, width(p.cellBoundaries2)],...
                'VariableNames', p.cellBoundaries2.Properties.VariableNames,...
                'VariableTypes', varfun(@class, p.cellBoundaries2, 'output','cell'));
            newCells{:,p.channels} = false(n, numel(p.channels));
            newCells.cellBB = zeros(n, 4);
            if ismember('colors', p.cellBoundaries2.Properties.VariableNames)
                newCells.colors = zeros(n, 3);
            end
            p.cellBoundaries2 = [p.cellBoundaries2; newCells];
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
            dapiCC = bwconncomp(boundaryStitch);
            p.dapiRP = regionprops(dapiCC);
            p.dapiLabelMat = labelmatrix(dapiCC);
            
            tmpBoundaries = bwboundaries(boundaryStitch);
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            nucBoundariesArray = cellfun(@(x) polyshape(cell2mat(x)), tmpBoundaries, 'UniformOutput', false); 
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
            status = true(numel(nucBoundariesArray),numel(p.channels));
            p.nucBoundaries = cell2table([num2cell((1:numel(nucBoundariesArray))'), nucBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
        end
        
        function p = loadCellPoseDapi(p, labelMatFile, outlineFile)
            tmpLabelMat = imread(labelMatFile);
            if all(size(tmpLabelMat) ==  p.scanObj.stitchDim)
                p.dapiLabelMat =  tmpLabelMat;
                scaleFactor = [1, 1];
            else
                p.dapiLabelMat = imresize(tmpLabelMat, p.scanObj.stitchDim, 'nearest');
                scaleFactor = p.scanObj.stitchDim./size(tmpLabelMat);
            end
            
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            warning('off', 'MATLAB:polyshape:boundary3Points')
            polyVect = d2utils.parseCellposeOutlines(outlineFile, 'scaleFactor', scaleFactor, 'flip', true);
            polyVect = scale(polyVect, scaleFactor);

            tmpArea = num2cell([polyVect.area]);
            [tmpX, tmpY] = centroid(polyVect);
            tmpCentroid = num2cell([tmpX', tmpY'],2);
            nucBoundariesArray = num2cell(polyVect);
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
%             tmpArea = cellfun(@(x) polyarea(x(:,1), x(:,2))*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpCentroid = cellfun(@(x) d2utils.poly2centroid(x)*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpBB = cellfun(@(x) d2utils.polygonBoundingBox2(x)*scaleFactor,polyArray, 'UniformOutput', false);
            p.dapiRP = cell2struct([tmpArea', tmpCentroid, tmpBB'], {'Area', 'Centroid', 'BoundingBox'}, 2);
            
%             nucBoundariesArray = cellfun(@(x) single(x.Vertices), nucBoundariesArray, 'UniformOutput', false);
%             nucBoundariesHeight = cellfun(@(x) height(x), nucBoundariesArray, 'UniformOutput', true);
%             nucIDArray = repelem(1:numel(nucBoundariesHeight), nucBoundariesHeight);
%             p.nucBoundaries = array2table([single(nucIDArray'), vertcat(nucBoundariesArray{:})], 'VariableNames', {'cellID', 'x', 'y'});
            %Try storying boundaries as polyshape.
            status = true(numel(nucBoundariesArray),numel(p.channels));
            p.nucBoundaries = cell2table([num2cell((1:numel(nucBoundariesArray))'), nucBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
%             p.nucBoundaries.status = true(numel(nucBoundariesArray), 1);
%             
            polyVectDilated = polybuffer(polyVect, p.radius);
            cellBoundariesArray = num2cell(polyVectDilated);
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), cellBoundariesArray, 'UniformOutput', false);
%             cellBoundariesArray = cellfun(@(x) x.Vertices, cellBoundariesArray, 'UniformOutput', false);
%             cellBoundariesHeight = cellfun(@(x) height(x), cellBoundariesArray, 'UniformOutput', true);
%             cellIDArray = repelem(1:numel(cellBoundariesHeight), cellBoundariesHeight);
%             p.cellBoundaries = array2table([single(cellIDArray'), vertcat(cellBoundariesArray{:})], 'VariableNames', {'cellID', 'x', 'y'});
            status = true(numel(cellBoundariesArray),numel(p.channels));
            p.cellBoundaries = cell2table([num2cell((1:numel(cellBoundariesArray))'), cellBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
            p.cellBoundaries.status = true(numel(cellBoundariesArray), 1);
%             warning('on', 'MATLAB:polyshape:repairedBySimplify')
            warning('on', 'MATLAB:polyshape:boundary3Points')
        end
        
        function p = labelMat2nucTable(p, varargin)
            if nargin == 2
                tmpLabelMat = imread(labelMatFile);
                if all(size(tmpLabelMat) ==  p.scanObj.stitchDim)
                    p.dapiLabelMat =  tmpLabelMat;
                else
                    p.dapiLabelMat = imresize(tmpLabelMat, p.scanObj.stitchDim, 'nearest');
                end
            end
            nucBoundariesTmp = cell(0, numel(p.dapiRP));
            for i = 1:numel(p.dapiRP) %can make this parfor?
                tmpBB = p.dapiRP(i).BoundingBox;
                %Create buffer
                tmpShift = ceil(tmpBB(3:4)*0.1);
                tmpStart = max([1,1], tmpBB(1:2)-tmpShift);
                tmpEnd = min(tmpStart+tmpBB(3:4)+(tmpShift*2), size(p.dapiLabelMat)); %Could make this stitchDim
                
                tmpRegionMask = p.dapiLabelMat(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
                tmpDapiMask = tmpRegionMask == i;
                tmpNucBoundary = bwboundaries(tmpDapiMask, 'noholes');
                nucBoundariesTmp{i} = cellfun(@(x) x+tmpStart-1, tmpNucBoundary, 'UniformOutput', false);
            end
            %Update nucBoundaries
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            nucBoundariesArray = cellfun(@(x) polyshape(cell2mat(x)), nucBoundariesTmp, 'UniformOutput', false); 
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
            status = true(numel(nucBoundariesArray), numel(p.channels));
            p.nucBoundaries2 = cell2table([num2cell((1:numel(nucBoundariesArray))'), nucBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function p = removeNucHoles(p)
            numHoles = [p.nucBoundaries2.nucBoundary.NumHoles];
            p.addEmptyRows(sum(numHoles)); %Technically don't need to add rows to cell boundaries.
            tmpNuc = p.nucBoundaries2(numHoles>0, :);
            for i = height(tmpNuc)
                tmpHoles = holes(tmpNuc.nucBoundary(i));
                for ii = 1:numel(p.channels)
                    for iii = numel(tmpHoles)
                        p.maskObj.addMaskLocalCoords(fliplr(tmpHoles(iii).Vertices), p.channels{ii})
                    end
                end
                tmpNuc.nucBoundary(i) = rmholes(tmpNuc.nucBoundary(i));    
            end
        end
        
        function outCellBoundaries = getCellBoundariesInRect(p, channel, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.cellBoundaries2.cellBB,rect)>0 & p.cellBoundaries2{:,channel};
            
            outCellBoundaries = p.cellBoundaries2(idx,:);
        end
        
        function outNucBoundaries = getNucBoundariesInRect(p,channel, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.nucBoundaries2.nucBB,rect)>0 & p.nucBoundaries2{:,channel};
            
            outNucBoundaries = p.nucBoundaries2(idx,:);
        end
        
        function outNucBoundaries = getAllNucBoundariesInRect(p, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.nucBoundaries2.nucBB,rect)>0;
            
            outNucBoundaries = p.nucBoundaries2(idx,:);
        end
        
        function outCellBoundaries = getAllCellBoundariesInRect(p, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.cellBoundaries2.cellBB,rect)>0;
            
            outCellBoundaries = p.cellBoundaries2(idx,:);
        end
        
        function p = addColors(p)
            %Note: 1) assumes equal number and order of cellIDs in nucBoundaries2 and cellBoundaries2 and 
            % 2) assigns colors to all cells and nuclei, even if status = false. 
            colorArray = [repmat(p.randomColors, floor(height(p.nucBoundaries2)/height(p.randomColors)), 1);...
                p.randomColors(1:mod(height(p.nucBoundaries2), height(p.randomColors)), :)];
           
            p.nucBoundaries2.colors = colorArray;
            [~, cellIdx] = ismember(p.nucBoundaries2.cellID, p.cellBoundaries2.cellID);
            if ~isempty(p.cellBoundaries2)
%                 p.cellBoundaries2.colors = colorArray;
                p.cellBoundaries2.colors = zeros(numel(p.cellBoundaries2.cellID), 3);
                p.cellBoundaries2.colors(cellIdx,:) = colorArray;
            end
        end
        
        function p = addNewCell(p, newCellID, nucPoly, cellPoly)
            nucBB = d2utils.polyshapeBoundingBox(nucPoly);
            cellBB = d2utils.polyshapeBoundingBox(cellPoly);
            idx = find(p.nucBoundaries2.cellID == 0, 1, 'first');
            randomColor = single(p.randomColors(randi(50),:));
            p.nucBoundaries2(idx,:) = [{newCellID, nucPoly, nucBB}, num2cell(true(1, numel(p.channels))), {randomColor}];
            p.cellBoundaries2(idx,:) = [{newCellID, cellPoly, cellBB}, num2cell(true(1, numel(p.channels))), {randomColor}];
        end
        
        function p = updateNucPoly(p, cellID, nucPoly)
            nucBB = d2utils.polyshapeBoundingBox(nucPoly);
            nucBoundaryIdx = p.nucBoundaries2.cellID == cellID; 
            p.nucBoundaries2(nucBoundaryIdx, {'nucBoundary', 'nucBB'}) = {nucPoly, nucBB};
            p.nucBoundaries2{nucBoundaryIdx, p.channels} = true(1, numel(p.channels));
        end
        
        function p = updateCellPoly(p, cellID, cellPoly)
            cellBB = d2utils.polyshapeBoundingBox(cellPoly);
            cellBoundaryIdx = p.cellBoundaries2.cellID == cellID; 
            p.cellBoundaries2(cellBoundaryIdx, {'cellBoundary', 'cellBB'}) = {cellPoly, cellBB};
            p.cellBoundaries2{cellBoundaryIdx, p.channels} = true(1, numel(p.channels));
        end
        
        function cellIDToMask = addCellMask(p, channel, maskPosition)
            tmpBB = d2utils.polygonBoundingBox(fliplr(maskPosition));
            tmpMaskPoly = polyshape(fliplr(maskPosition));
            
            nucTableTmp = p.getNucBoundariesInRect(channel, tmpBB);
            cellTableTmp = p.getCellBoundariesInRect(channel, tmpBB);
            
            nucIDs = nucTableTmp.cellID(overlaps(tmpMaskPoly, nucTableTmp.nucBoundary));
            cellIDs = cellTableTmp.cellID(overlaps(tmpMaskPoly, cellTableTmp.cellBoundary));
            cellIDToMask = union(nucIDs, cellIDs);
            
            nucIdx = ismember(p.nucBoundaries2.cellID, cellIDToMask);
            p.nucBoundaries2{nucIdx, channel} = false;
            
            cellIdx = ismember(p.cellBoundaries2.cellID, cellIDToMask);
            p.cellBoundaries2{cellIdx, channel} = false;
        end
        
        function p = saveBoundaries(p)
            if isempty(p.nucBoundaries2) || ~any(p.nucBoundaries2{:,p.channels}, 'all')
                fprintf("There are no valid nuclei boundaries to save")
            else
                outTable = p.nucBoundaries2(~p.nucBoundaries2.cellID == 0,:); %Remove empty rows but keep masked cells (status = false). 
                outTable = rowfun(@(x) x.Vertices, outTable, 'InputVariables', 'nucBoundary', 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'Vertices');
                outTable(:,'GroupCount') = [];
                outTable = splitvars(outTable, 'Vertices', 'NewVariableNames', {'x', 'y'});
                writetable(outTable, p.nucBoundariesFile);
            end
            
            if isempty(p.cellBoundaries2) || ~any(p.cellBoundaries2{:,p.channels}, 'all')
                fprintf("There are no valid cell boundaries to save")
            else
                outTable = p.cellBoundaries2(~p.cellBoundaries2.cellID == 0,:); %Remove empty rows but keep masked cells (status = false). 
                outTable = rowfun(@(x) x.Vertices, outTable, 'InputVariables', 'cellBoundary', 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'Vertices');
                outTable(:,'GroupCount') = [];
                outTable = splitvars(outTable, 'Vertices', 'NewVariableNames', {'x', 'y'});
                writetable(outTable, p.cellBoundariesFile);
            end
        end

%         function p = maskCellsInChannel(p, channel)
%             maxCellMask = max(p.maskObj.masksBB{p.maskObj.masksBB.dapi,'maskID'});
%             maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxCellMask,'BB'}; %Only query nuclei within mask bouding box
%             tmpMaskPoly = polyshape(p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask,{'x', 'y'}});
%             nucTableTmp = p.getNucBoundariesInRect(channel, maskBB);
%             cellTableTmp = p.getCellBoundariesInRect(channel, maskBB);
%             
%             nucIdx = p.getNucleiInRectIdx(maskBB);
%             polyIdx = inpolygon(p.nuclei.x(nucIdx), p.nuclei.y(nucIdx),...
%                 p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'y'});
%             
%             nucIdx(nucIdx) = polyIdx; %Only nuclei in polygon remain true
%             p.nuclei.maskID(nucIdx) = maxCellMask;
%             p.nuclei.status(nucIdx) = false;
%             
%         end
%         
%         function p = addMask(p, channel, rect, maskCoords)
%             tmpMaskPoly = polyshape(maskCoords);
%             nucTableTmp = p.getNucBoundariesInRect(channel, rect);
%             cellTableTmp = p.getCellBoundariesInRect(channel, rect);
%             
%             nucIDs = nucTableTmp.cellID(overlaps(tmpMaskPoly, nucTableTmp.nucBoundary));
%             cellIDs = cellTableTmp.cellID(overlaps(tmpMaskPoly, cellTableTmp.cellBoundary));
%             cellIDToMask = union(nucIDs, cellIDs);
%             mask
%                 tmpCellIDs = union(nucTableTmp.cellID, cellTableTmp.cellID);
%         end
    end
end