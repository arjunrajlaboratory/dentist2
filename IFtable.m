classdef IFtable < handle
    
    properties (Access = public)
        
        scanObj
        maskObj %Masks for IF signal. Not for cell/nuclei boundaries. 
        IFboundaries
        
        IFquant
        IFquant2
        channels
        centroidLists
        expressionColorPal = {'BuYlRd', 'YlOrRd', 'GrBu', 'BuGnYlRd'}
        paletteIdx = 1;
        minNucleusSize = 1000;
        radius = 20;
        IFfile
    end
    
    methods
        function p = IFtable(scanObject, maskObj, IFboundaries, varargin)
            p.scanObj = scanObject;
            p.maskObj = maskObj;
            p.IFboundaries = IFboundaries;
            p.channels = p.scanObj.stitchedScans.labels;
            if nargin == 3
                fprintf('New Table\n');
                p.IFquant = table('size', [0,numel(p.channels)+ 9],... %Possibly unnecessary since new table created with p.findNuclei. 
                    'VariableNames', [{'cellID', 'x', 'y', 'status'}, p.channels, {'meanNuc', 'meanCyto',  'sumNuc', 'sumCyto','maskID'}],...
                    'VariableTypes', [repmat({'single'}, 1, 3), repmat({'logical'}, 1, 1+numel(p.channels)), repmat({'single'}, 1, 5)]);
            elseif nargin == 4
                fprintf('Loading Table\n');
                p.IFfile = varargin{1};
                opts = detectImportOptions(varargin{1});
                opts = setvartype(opts, 'single');
                p.IFquant = readtable(varargin{1},opts);
                p.IFquant(:,[{'status'},  p.channels]) = logical(p.cells(:,[{'status'},  p.channels])); %For some reason, when I set 'status' to 'logical' they all go to false. So doing this instead
            end
        end
        
        function p = addEmptyRows(p, n)
            newRows = table('Size', [n, width(p.IFquant)],...
                'VariableNames', p.IFquant.Properties.VariableNames,...
                'VariableTypes', varfun(@class, p.IFquant, 'output','cell'));
            newRows.status = false(n, 1);
            p.IFquant = [p.IFquant; newRows];
        end
        
        function p = quantAllLabelMat(p, varargin)
        %Uses polyshapes to mask nuclei and cells.
            if nargin == 1
                channelsToQaunt = p.channels;
            elseif nargin ==2
                channelsToQaunt = varargin{1};
            end
            
            nRows = prod(numel(p.IFboundaries.dapiRP), numel(channelsToQaunt));
            meanNuc = zeros(nRows,1);
            meanCyto = zeros(nRows,1);
            sumNuc = zeros(nRows,1);
            sumCyto = zeros(nRows,1);
            channelStatus = false(nRows,numel(p.channels));
            nucBoundariesTmp = cell(0, numel(p.IFboundaries.dapiRP));
            cellBoundariesTmp = cell(0, numel(p.IFboundaries.dapiRP));
            for i = 1:numel(p.IFboundaries.dapiRP) %can make this parfor?
                tmpBB = p.IFboundaries.dapiRP(i).BoundingBox;
                tmpSize = tmpBB(3:4) + (2*p.radius);
                tmpStart = max([1,1], tmpBB(1:2)-p.radius-5); %Add a buffer beyond the p.radius
                tmpEnd = min(tmpStart+tmpSize+10, size(p.IFboundaries.dapiLabelMat)); %Could make this stitchDim
                
                tmpRegionMask = p.IFboundaries.dapiLabelMat(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
                tmpDapiMask = tmpRegionMask == i;
                %Consider switching from masks to polyshapes below (and
                %expanding polyshape with polybuffer instead of imdilate). 
                tmpDapiMaskDilated = imdilate(tmpDapiMask, strel('disk', p.radius)); %Consider other strel. Including tmpDapiMask. 
                tmpCytoMask = tmpDapiMaskDilated & ~logical(tmpRegionMask);
                tmpCellMask = or(tmpCytoMask, tmpDapiMask);
                tmpNucBoundary = bwboundaries(tmpDapiMask, 'noholes');
                tmpCellBoundary = bwboundaries(tmpCellMask, 'noholes');
                nucBoundariesTmp{i} = cellfun(@(x) x+tmpStart-1, tmpNucBoundary, 'UniformOutput', false);
                cellBoundariesTmp{i} = cellfun(@(x) x+tmpStart-1, tmpCellBoundary, 'UniformOutput', false);
                for ii = 1:numel(channelsToQaunt)
                    stitchIdx = ismember(p.channels, channelsToQaunt(ii)); %Redundant in some cases but flexible
                    tmpImage = p.scanObj.stitchedScans.stitches{stitchIdx}(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
                    tmpMasks = p.maskObj.getChannelMasksInRect([tmpStart, tmpEnd - tmpStart], channelsToQaunt(ii));
                    if ~isempty(tmpMasks)
                        maskIDs = unique(tmpMasks.maskID);
                        for iii = 1:numel(maskIDs)
                            bwmask = poly2mask(tmpMasks{tmpMasks.maskID == maskIDs(iii), 'x'}, tmpMasks{tmpMasks.maskID == maskIDs(iii), 'y'}, size(tmpDapiMask, 1), size(tmpDapiMask, 2));
                            tmpDapiMask(bwmask) = false;
                            tmpCytoMask(bwmask) = false;
                        end
                    end
                    tmpIdx = (i-1)*numel(channelsToQaunt)+ii;
                    meanNuc(tmpIdx) = mean(tmpImage(tmpDapiMask));
                    sumNuc(tmpIdx) = sum(tmpImage(tmpDapiMask));
                    meanCyto(tmpIdx) = mean(tmpImage(tmpCytoMask));
                    sumCyto(tmpIdx) = sum(tmpImage(tmpCytoMask));
                    channelStatus(tmpIdx,:) = stitchIdx;
                end
            end
            %Update nucBoundaries
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            nucBoundariesArray = cellfun(@(x) polyshape(cell2mat(x)), nucBoundariesTmp, 'UniformOutput', false); 
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
            status = true(numel(nucBoundariesArray), numel(p.channels));
            p.IFboundaries.nucBoundaries2 = cell2table([num2cell((1:numel(nucBoundariesArray))'), nucBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
            %Update cellBoundaries
            cellBoundariesArray = cellfun(@(x) polyshape(cell2mat(x)), cellBoundariesTmp, 'UniformOutput', false);
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), cellBoundariesArray, 'UniformOutput', false);
            status = true(numel(cellBoundariesArray),numel(p.channels));
            p.IFboundaries.cellBoundaries2 = cell2table([num2cell((1:numel(cellBoundariesArray))'), cellBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
            p.IFboundaries.addColors();
            p.IFboundaries.addEmptyRows(1000);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
            %Update IFquant
            [cellCoordsX, cellCoordsY]  = cellfun(@(x) centroid(x, 1:x.NumRegions), nucBoundariesArray, 'UniformOutput', false); 
            cellCoords = cellfun(@(x) single(round(mean(x, 1))), [cellCoordsX; cellCoordsY], 'UniformOutput', false); 
            cellCoords = repelem(cellCoords', numel(p.channels), 1);
            status = true(nRows,1);
            maskID = single(zeros(nRows,1));
            cellIDs = repelem(1:numel(p.IFboundaries.dapiRP), numel(p.channels));
            p.IFquant = array2table([cellIDs',cell2mat(cellCoords), status, channelStatus, meanNuc, meanCyto, sumNuc, sumCyto, maskID],...
                'VariableNames', [{'cellID', 'x', 'y', 'status'}, p.channels, {'meanNuc', 'meanCyto',  'sumNuc', 'sumCyto','maskID'}]);
            p.IFquant = convertvars(p.IFquant, [{'status'}, p.channels], 'logical');
            p.addEmptyRows(1000);
        end
        
        function p = quantAllLabelMat2(p, varargin)
            %Uses polybuffer instead of imdilate
            if nargin == 1
                channelsToQaunt = p.channels;
            elseif nargin ==2
                channelsToQaunt = varargin{1};
            end
            
            nRows = prod(numel(p.IFboundaries.dapiRP), numel(channelsToQaunt));
            meanNucArray = zeros(nRows,1);
            meanCytoArray = zeros(nRows,1);
            sumNucArray = zeros(nRows,1);
            sumCytoArray = zeros(nRows,1);
            channelStatus = false(nRows,numel(p.channels));
            cellBoundariesTmp = cell(0, numel(p.IFboundaries.dapiRP));
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            for i = 1:numel(p.IFboundaries.dapiRP) %can make this parfor?
                tmpNucPoly = p.IFboundaries.nucBoundaries2.nucBoundary(p.IFboundaries.nucBoundaries2.cellID == i);
                tmpCellPoly = tmpNucPoly.polybuffer(p.radius); %What happens if disjoint nuclei? 
                tmpCellPoly = polyshape(round(tmpCellPoly.Vertices));
                tmpCellPoly = rmholes(tmpCellPoly);
                tmpBB = d2utils.polyshapeBoundingBox(tmpCellPoly, p.scanObj.stitchDim);
                tmpNucInRect = p.IFboundaries.getAllNucBoundariesInRect(tmpBB);
                %Subtract neighboring nuclei (but not own nucleus).
                %Should merge neighboring nuclei if they fall entirely within bounds
                %of the cellPoly. 
                tmpNucInRect = subtract(union(tmpNucInRect.nucBoundary), tmpNucPoly);
                tmpCellPoly = subtract(tmpCellPoly,tmpNucInRect);
                tmpDapiMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpNucPoly), tmpBB(1:2), 'flip', true);
                tmpCellMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpCellPoly), tmpBB(1:2), 'flip', true);
                tmpCytoMask = tmpCellMask & ~ tmpDapiMask;
                cellBoundariesTmp{i} = tmpCellPoly;
                for ii = 1:numel(channelsToQaunt)
                    [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, channelsToQaunt(ii));
                    tmpIdx = (i-1)*numel(channelsToQaunt)+ii;
                    meanNucArray(tmpIdx) = meanNuc;
                    sumNucArray(tmpIdx) = sumNuc;
                    meanCytoArray(tmpIdx) = meanCyto;
                    sumCytoArray(tmpIdx) = sumCyto;
                    stitchIdx = ismember(p.channels, channelsToQaunt(ii)); %Redundant in some cases but flexible
                    channelStatus(tmpIdx,:) = stitchIdx;
                end                
            end
            %Update cellBoundaries
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), cellBoundariesTmp, 'UniformOutput', false);
            status = true(numel(cellBoundariesTmp),numel(p.channels));
            p.IFboundaries.cellBoundaries2 = cell2table([num2cell((1:numel(cellBoundariesTmp))'), cellBoundariesTmp', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
            p.IFboundaries.addColors();
            %Update IFquant
            [cellCoordsX, cellCoordsY]  = arrayfun(@(x) centroid(x, 1:x.NumRegions),p.IFboundaries.nucBoundaries2.nucBoundary, 'UniformOutput', false); 
            cellCoords = cellfun(@(x) single(round(mean(x, 1))), [cellCoordsX, cellCoordsY], 'UniformOutput', false); 
            cellCoords = repelem(cellCoords, numel(p.channels), 1);
            status = true(nRows,1);
            maskID = single(zeros(nRows,1));
            cellIDs = repelem(1:numel(p.IFboundaries.dapiRP), numel(p.channels));
            p.IFquant = array2table([cellIDs',cell2mat(cellCoords), status, channelStatus, meanNucArray, meanCytoArray, sumNucArray, sumCytoArray, maskID],...
                'VariableNames', [{'cellID', 'x', 'y', 'status'}, p.channels, {'meanNuc', 'meanCyto',  'sumNuc', 'sumCyto','maskID'}]);
            p.IFquant = convertvars(p.IFquant, [{'status'}, p.channels], 'logical');
            p.IFboundaries.addEmptyRows(1000);
            p.addEmptyRows(1000);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function outCells = getCellsInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.IFquant.x >= rect(1) & p.IFquant.x < rect(1) + rect(3) ...
                & p.IFquant.y >= rect(2) & p.IFquant.y < rect(2) + rect(4) ...
                & p.IFquant.status;
            
            outCells = p.IFquant(idx,:);
        end
       
        function p = addNucleus(p, polyXY, rect, andCyto)
            
            if ~any(p.IFquant.cellID == 0)
                p.addEmptyRows(1000);
                p.IFboundaries.addEmptyRows(1000)
            end
            %See if nucleus falls within cell/nuclei mask. If so, return. 
            
            %See if new nucleus falls within known cell boundary.
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            tmpPoly = polyshape(polyXY+rect(1:2));
            polyBB = d2utils.polyshapeBoundingBox(tmpPoly);
            cellBoundariesInView = p.IFboundaries.getAllCellBoundariesInRect(polyBB);
            nucBoundariesInView = p.IFboundaries.getAllNucBoundariesInRect(polyBB);

            if isempty(cellBoundariesInView)
                cellOverlapIdx = false;
            else
                cellOverlapIdx = overlaps(tmpPoly, cellBoundariesInView.cellBoundary);
            end
            
            if isempty(nucBoundariesInView)
                nucOverlapIdx = false;
            else
                nucOverlapIdx = overlaps(tmpPoly, nucBoundariesInView.nucBoundary);
            end
            if any(cellOverlapIdx)
                tempCellIDs = cellBoundariesInView.cellID(cellOverlapIdx);
                %For each cell, mask nuclei and save nuclei boundaries
                for i = 1:numel(tempCellIDs) 
                    p.addNucToCell(tmpPoly, tempCellIDs(i), cellBoundariesInView);
                end
            elseif any(nucOverlapIdx)
                p.addNucToNuc(tmpPoly, nucBoundariesInView(nucOverlapIdx,:));
            else
                if andCyto
                    cellPoly = tmpPoly.polybuffer(p.radius);
                    cellPoly = polyshape(round(cellPoly.Vertices));
                    %Need to subtract neighboring nuclei
                    tmpNucleiInRect = p.IFboundaries.getAllNucBoundariesInRect(tmpBB);
                    cellPoly = subtract(cellPoly,tmpNucleiInRect.nucBoundary);
                else
%                     tmpCytoMask = false(size(tmpDapiMask));
                    cellPoly = polyshape();
                end
                p.addNewCell(tmpPoly, cellPoly);
            end
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function p = addNewCell(p, nucPoly, cellPoly, varargin)
            if nargin > 3 %Option to input dapi mask and cyto mask and rect
                tmpDapiMask = varargin{1};
                tmpCytoMask = varargin{2};
                tmpBB = varargin{3};
            else
                tmpBB = d2utils.polyshapeBoundingBox(union(nucPoly, cellPoly)); %Union in case nuclei boundary extends beyond cell boundary
                tmpDapiMask = d2utils.polyshape2mask(nucPoly, tmpBB(1:2), tmpBB(3:4), 'flip', true);
                tmpCellMask = d2utils.polyshape2mask(cellPoly, tmpBB(1:2), tmpBB(3:4), 'flip', true);
                tmpCytoMask = tmpCellMask & ~tmpDapiMask;
            end
            nRows = numel(p.channels);
            meanNucArray = zeros(nRows, 1);
            meanCytoArray = zeros(nRows, 1);
            sumNucArray = zeros(nRows, 1);
            sumCytoArray = zeros(nRows, 1);
            channelStatus = false(nRows);
            [tmpX, tmpY] = centroid(nucPoly);
            cellCoords = repmat(round([tmpX, tmpY]), nRows,1); %Does this need to be flipped?
            for i = 1:nRows
                [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, p.channels{i});
                meanNucArray(i) = meanNuc;
                meanCytoArray(i) = meanCyto;
                sumNucArray(i) = sumNuc;
                sumCytoArray(i) = sumCyto;
                channelStatus(i,:) = ismember(p.channels(i), p.channels);
            end
            status = true(nRows,1);
            maskID = single(zeros(nRows,1));
            newCellID = max(p.IFquant.cellID)+1;
            cellIDs = repelem(newCellID, nRows, 1);
            newCell = array2table([cellIDs, cellCoords, status, channelStatus, meanNuc, meanCyto, sumNuc, sumCyto, maskID],...
                'VariableNames', [{'cellID', 'x', 'y', 'status'}, p.channels, {'meanNuc', 'meanCyto',  'sumNuc', 'sumCyto','maskID'}]);
            newCell =  convertvars(newCell, [{'status'}, p.channels], 'logical');
            startIdx = find(p.IFquant.cellID == 0, 1, 'first');
            p.IFquant(startIdx:startIdx+nRows-1,:) = newCell;
            %Update boundaries
            p.IFboundaries.addNewCell(newCellID, nucPoly, cellPoly); 
        end
        
        function p = addNucToCell(p, nucPoly, cellID, cellBoundariesInView)
            tmpCellBoundary = cellBoundariesInView{cellBoundariesInView.cellID == cellID, 'cellBoundary'};
            priorNuclei = p.IFboundaries.nucBoundaries2(p.IFboundaries.nucBoundaries2.cellID == cellID, :);
            tmpNucPoly = intersect(tmpCellBoundary, nucPoly);
            tmpNucPoly = union([priorNuclei.nucBoundary; tmpNucPoly]);
            tmpBB = d2utils.polyshapeBoundingBox(union(tmpCellBoundary, tmpNucPoly)); %Union in case nuclei boundary extends beyond cell boundary
            tmpCellMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpCellBoundary), tmpBB(1:2), 'flip', true);
            
%             priorDAPIMask = d2utils.polyshape2mask(priorNuclei.nucBoundary, tmpStart-1, tmpCellBoundary.cellBB(3:4), 'flip', true);
            tmpDapiMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpNucPoly), tmpBB(1:2), 'flip', true);
%             tmpDapiMask = or(tmpDapiMask, priorDAPIMask);
%             tmpCytoPoly = subtract(tmpCellBoundary, tmpNucPoly);
            tmpCytoMask = tmpCellMask & ~tmpDapiMask;
%             tmpCytoMask = d2utils.polyshape2mask(tmpCytoPoly, tmpStart-1, tmpCellBoundary.cellBB(3:4), 'flip', true);
            %Update nucBoundaries
            p.IFboundaries.updateNucPoly(cellID, tmpNucPoly);
            
            %Find new centroids
            [tmpX, tmpY] = centroid(tmpNucPoly, 1:tmpNucPoly.NumRegions);
            tmpCentroid = round(mean([tmpX, tmpY], 1)); %Does this need to be flipped? 

            %Update IF quant table
            cellIdx = p.IFquant.cellID == cellID;
            p.IFquant{cellIdx, {'x', 'y'}} = repmat(tmpCentroid, numel(p.channels), 1);
            for i = 1:numel(p.channels)
                [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, p.channels{i});
                cellChannelIdx = cellIdx & p.IFquant{:,p.channels{i}};
                p.IFquant{cellChannelIdx,{'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}} = [meanNuc, meanCyto, sumNuc, sumCyto];
            end
        end
        
        function p = addNucToNuc(p, nucPoly,nucBoundaries)
            %Adds new nucleus poly but doesn't create new cell
            %If we want an option to add nuc+cell here, can use polybuffer
            %to create cellpoly from nucpoly, then take union of cellpoly
            %with prior cellboundaries. 
            tmpNucPoly = union(nucPoly, nucBoundaries.nucBoundary);
            tmpCellIdx = ismember(p.IFboundaries.cellBoundaries2.cellID, nucBoundaries.cellID);
            tmpCellBoundary = union(p.IFboundaries.cellBoundaries2.cellBoundary(tmpCellIdx));
            cellID = nucBoundaries.cellID(1);

            tmpBB = d2utils.polyshapeBoundingBox(union(tmpCellBoundary, tmpNucPoly)); %Union in case nuclei boundary extends beyond cell boundary
        
            tmpCellMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpCellBoundary), tmpBB(1:2), 'flip', true);

            tmpDapiMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpNucPoly), tmpBB(1:2), 'flip', true);
            tmpCytoMask = tmpCellMask & ~tmpDapiMask;
            
            
            %Update IF quant table
            cellIdx = p.IFquant.cellID == cellID;
            for i = 1:numel(p.channels)
                [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, p.channels{i});
                cellChannelIdx = cellIdx & p.IFquant{:,p.channels{i}};
                p.IFquant{cellChannelIdx,{'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}} = [meanNuc, meanCyto, sumNuc, sumCyto];
            end
            %Find new centroids
            [tmpX, tmpY] = centroid(tmpNucPoly, 1:tmpNucPoly.NumRegions);
            tmpCentroid = round(mean([tmpX, tmpY], 1));
            p.IFquant{cellIdx, {'x', 'y'}} = repmat(tmpCentroid, numel(p.channels), 1);
                
            %Update nucBoundaries
            p.IFboundaries.updateNucPoly(cellID, tmpNucPoly);
          
            %If cell merges multiple cellIDs 
            if numel(nucBoundaries.cellID) > 1
                %Update cellBoundaries
                p.IFboundaries.updateCellPoly(cellID, tmpCellBoundary);
                
                %Delete redundant cellIDs
                cellIDsToDelete = nucBoundaries.cellID(2:end);
                p.IFboundaries.nucBoundaries2(ismember(p.IFboundaries.nucBoundaries2.cellID, cellIDsToDelete),:) = [];
                p.IFboundaries.cellBoundaries2(ismember(p.IFboundaries.cellBoundaries2.cellID, cellIDsToDelete),:) = [];
                p.IFquant(ismember(p.IFquant.cellID, cellIDsToDelete),:) = [];
            end
        end
        
        function p = addCell(p, polyXY, rect)
            
            if ~any(p.IFquant.cellID == 0)
                p.addEmptyRows(1000);
                p.IFboundaries.addEmptyRows(1000)
            end
            
            %See if new cell falls within known cell or nucleus boundary.
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            tmpPoly = polyshape(polyXY+rect(1:2));
            polyBB = d2utils.polyshapeBoundingBox(tmpPoly);
            cellBoundariesInView = p.IFboundaries.getAllCellBoundariesInRect(polyBB);
            nucBoundariesInView = p.IFboundaries.getAllNucBoundariesInRect(polyBB);
            
            if isempty(cellBoundariesInView)
                cellOverlapIdx = false;
            else
                cellOverlapIdx = overlaps(tmpPoly, cellBoundariesInView.cellBoundary);
            end
            
            if isempty(nucBoundariesInView)
                nucOverlapIdx = false;
            else
                nucOverlapIdx = overlaps(tmpPoly, nucBoundariesInView.nucBoundary);
            end
            
            if any(cellOverlapIdx)
                tempCellIDs = cellBoundariesInView.cellID(cellOverlapIdx);
                p.addCellToCell(tmpPoly, tempCellIDs, cellBoundariesInView);
            elseif any(nucOverlapIdx)
                 tempCellIDs = nucBoundariesInView.cellID(nucOverlapIdx);
                 p.addCellToNuc(tmpPoly, tempCellIDs, nucBoundariesInView);
            else
                emptyNuc = polyshape();
                p.addNewCell(emptyNuc, tmpPoly);
            end
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function p = addCellToCell(p, cellPoly, tempCellIDs, cellBoundariesInView)
            tmpCellBoundary = union([cellBoundariesInView{ismember(cellBoundariesInView.cellID, tempCellIDs), 'cellBoundary'};cellPoly]); %Consider searching all cell boundaries if new poly somehow extend beyond the view rect. 
            tmpNucleiIdx = ismember(p.IFboundaries.nucBoundaries2, tempCellIDs);
            tmpNucBoundary = union(p.IFboundaries.nucBoundaries2.nucBoundary(tmpNucleiIdx));
            cellID = tempCellIDs(1);
            
            tmpBB = d2utils.polyshapeBoundingBox(union(tmpCellBoundary, tmpNucBoundary)); %Union in case nuclei boundary extends beyond cell boundary
            tmpNucPoly = intersect(tmpNucBoundary, tmpCellBoundary);
            tmpCellMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpCellBoundary), tmpBB(1:2), 'flip', true);

            tmpDapiMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpNucPoly), tmpBB(1:2), 'flip', true);
            tmpCytoMask = tmpCellMask & ~tmpDapiMask;

            %Update IF quant table
            cellIdx = p.IFquant.cellID == cellID;
            for i = 1:numel(p.channels)
                [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, p.channels{i});
                cellChannelIdx = cellIdx & p.IFquant{:,p.channels{i}};
                p.IFquant{cellChannelIdx,{'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}} = [meanNuc, meanCyto, sumNuc, sumCyto];
            end
            
            %Update cellBoundaries
            p.IFboundaries.updateCellPoly(cellID, tmpCellBoundary);
          
            %If cell merges multiple cellIDs 
            if numel(tempCellIDs) > 1
                %Update nucBoundaries
                p.IFboundaries.updateNucPoly(cellID, tmpNucPoly);
                
                %Find new centroids
                [tmpX, tmpY] = centroid(tmpNucPoly, 1:tmpNucPoly.NumRegions);
                tmpCentroid = round(mean([tmpX, tmpY], 1)); 
                p.IFquant{cellIdx, {'x', 'y'}} = repmat(tmpCentroid, numel(p.channels), 1);

                %Delete redundant cellIDs
                cellIDsToDelete = tempCellIDs(2:end);
                p.IFboundaries.nucBoundaries2(ismember(p.IFboundaries.nucBoundaries2.cellID, cellIDsToDelete),:) = [];
                p.IFboundaries.cellBoundaries2(ismember(p.IFboundaries.cellBoundaries2.cellID, cellIDsToDelete),:) = [];
                p.IFquant(ismember(p.IFquant.cellID, cellIDsToDelete),:) = [];
            end
        end
        
        function p = addCellToNuc(p, tmpPoly, tempCellIDs, nucBoundariesInView)
            tmpNucBoundary = union(nucBoundariesInView{ismember(nucBoundariesInView.cellID, tempCellIDs), 'nucBoundary'}); %Consider searching all cell boundaries if new poly somehow extend beyond the view rect. 
            cellID = tempCellIDs(1);

            tmpBB = d2utils.polyshapeBoundingBox(union(tmpPoly, tmpNucBoundary)); %Union in case nuclei boundary extends beyond cell boundary
            tmpNucPoly = intersect(tmpNucBoundary, tmpPoly);
            tmpCellMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpPoly), tmpBB(1:2), 'flip', true);
            tmpDapiMask = d2utils.polyvect2mask(tmpBB(3:4), regions(tmpNucPoly), tmpBB(1:2), 'flip', true);
            tmpCytoMask = tmpCellMask & ~tmpDapiMask;

            %Update IF quant table
            cellIdx = p.IFquant.cellID == cellID;
            for i = 1:numel(p.channels)
                [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, p.channels{i});
                cellChannelIdx = cellIdx & p.IFquant{:,p.channels{i}};
                p.IFquant{cellChannelIdx,{'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}} = [meanNuc, meanCyto, sumNuc, sumCyto];
            end

            %Update cellBoundaries
            p.IFboundaries.updateCellPoly(cellID, tmpPoly);
            
            %If cell merges multiple cellIDs 
            if numel(tempCellIDs) > 1
                %Update nucBoundaries
                p.IFboundaries.updateNucPoly(cellID, tmpNucPoly);
                
                %Find new centroids
                [tmpX, tmpY] = centroid(tmpNucPoly, 1:tmpNucPoly.NumRegions);
                tmpCentroid = round(mean([tmpX, tmpY], 1)); %Does this need to be flipped?
                p.IFquant{cellIdx, {'x', 'y'}} = repmat(tmpCentroid, numel(p.channels), 1);

                %Delete redundant cellIDs
                cellIDsToDelete = tempCellIDs(2:end);
                p.IFboundaries.nucBoundaries2(ismember(p.IFboundaries.nucBoundaries2.cellID, cellIDsToDelete),:) = [];
                p.IFboundaries.cellBoundaries2(ismember(p.IFboundaries.cellBoundaries2.cellID, cellIDsToDelete),:) = [];
                p.IFquant(ismember(p.IFquant.cellID, cellIDsToDelete),:) = [];
            end
        end
        
        function p = deleteNuc(p, points, rect, inROI)
            if inROI
                tmpPoly = polyshape(points);
                tmpBB = d2utils.polygonBoundingBox(tmpPoly);
                nucsInView = p.IFboundaries.getAllNucBoundariesInRect(tmpBB);
                nucIdx = overlaps(tmpPoly, nucsInView.nucBoundary);
                cellIDsToDelete = nucsInView.cellID(nucIdx);
            else
                nucsInView = p.IFboundaries.getAllNucBoundariesInRect(rect);
                nucIdx = arrayfun(@(x) any(isinterior(x, points(:,2), points(:,1))), nucsInView.nucBoundary);
                cellIDsToDelete = nucsInView.cellID(nucIdx);
            end
            %Delete nuc poly, set all channels to false and requant cyto IF. 
            for i = 1:numel(cellIDsToDelete)
                if p.IFboundaries.cellBoundaries2{p.IFboundaries.cellBoundaries2.cellID == cellIDsToDelete(i), 'cellBoundary'}.NumRegions %Is the polyshape nonempty
                    nucIdx = p.IFboundaries.nucBoundaries2.cellID == cellIDsToDelete(i);
                    p.IFboundaries.nucBoundaries2{nucIdx,'nucBoundary'} = polyshape();
                    p.IFboundaries.nucBoundaries2{nucIdx,p.channels} = false(1,numel(p.channels));
                    p.requantCell(cellIDsToDelete(i));
                    %Update centroids. Consider functionalizing
                    tmpCellPoly =  p.IFboundaries.cellBoundaries2.cellBoundary(p.IFboundaries.cellBoundaries2.cellID == cellIDsToDelete(i));
                    [tmpX, tmpY] = centroid(tmpCellPoly, 1:tmpCellPoly.NumRegions);
                    tmpCentroid = round(mean([tmpX, tmpY], 1)); %Does this need to be flipped?
                    cellIdx = p.IFquant.cellID == cellIDsToDelete(i);
                    p.IFquant{cellIdx, {'x', 'y'}} = repmat(tmpCentroid, numel(p.channels), 1);
                else
                    %Delete cells entirely.
                    p.IFboundaries.nucBoundaries2(p.IFboundaries.nucBoundaries2.cellID == cellIDsToDelete(i), :) = [];
                    p.IFboundaries.cellBoundaries2(p.IFboundaries.cellBoundaries2.cellID == cellIDsToDelete(i), :) = [];
                    p.IFquant(p.IFquant.cellID == cellIDsToDelete(i),:) = [];
                end
            end
        end

        function p = deleteCell(p, points, rect, inROI)
            if inROI
                tmpPoly = polyshape(points);
                tmpBB = d2utils.polygonBoundingBox(tmpPoly);
                cellsInView = p.IFboundaries.getAllCellBoundariesInRect(tmpBB);
                cellIdx = overlaps(tmpPoly, cellsInView.cellBoundary);
                cellIDsToDelete = cellsInView.cellID(cellIdx);
            else
                cellsInView = p.IFboundaries.getAllCellBoundariesInRect(rect);
                cellIdx = arrayfun(@(x) any(isinterior(x, points(:,2), points(:,1))), cellsInView.cellBoundary);
                cellIDsToDelete = cellsInView.cellID(cellIdx);
            end
            for i = 1:numel(cellIDsToDelete)
                if p.IFboundaries.nucBoundaries2{p.IFboundaries.nucBoundaries2.cellID == cellIDsToDelete(i), 'nucBoundary'}.NumRegions %Is the polyshape nonempty
                    cellIdx = p.IFboundaries.cellBoundaries2.cellID == cellIDsToDelete(i);
                    p.IFboundaries.cellBoundaries2{cellIdx,'cellBoundary'} = polyshape();
                    p.IFboundaries.cellBoundaries2{cellIdx,p.channels} = false(1,numel(p.channels));
                    p.IFquant{p.IFquant.cellID == cellIDsToDelete(i),{'meanNuc','sumNuc'}} = zeros(numel(p.channels), 2);
                else
                    %Delete cells entirely.
                    p.IFboundaries.nucBoundaries2(p.IFboundaries.nucBoundaries2.cellID == cellIDsToDelete(i), :) = [];
                    p.IFboundaries.cellBoundaries2(p.IFboundaries.cellBoundaries2.cellID == cellIDsToDelete(i), :) = [];
                    p.IFquant(p.IFquant.cellID == cellIDsToDelete(i),:) = [];
                end
            end
        end
        
        function p = deleteNucAndCell(p, points, rect, inROI)
            if inROI
                tmpPoly = polyshape(points);
                tmpBB = d2utils.polygonBoundingBox(tmpPoly);
                cellsInView = p.IFboundaries.getAllCellBoundariesInRect(tmpBB);
                cellIdx = overlaps(tmpPoly, cellsInView.cellBoundary);
                nucsInView = p.IFboundaries.getAllNucBoundariesInRect(tmpBB);
                nucIdx = overlaps(tmpPoly, nucsInView.nucBoundary);
                cellIDsToDelete = union(nucsInView.cellID(nucIdx), cellsInView.cellID(cellIdx));
            else
                nucsInView = p.IFboundaries.getAllNucBoundariesInRect(rect);
                nucIdx = arrayfun(@(x) any(isinterior(x, points(:,2), points(:,1))), nucsInView.nucBoundary);
                cellsInView = p.IFboundaries.getAllCellBoundariesInRect(rect);
                cellIdx = arrayfun(@(x) any(isinterior(x, points(:,2), points(:,1))), cellsInView.cellBoundary);
                cellIDsToDelete = union(nucsInView.cellID(nucIdx), cellsInView.cellID(cellIdx));
            end
            %Delete cells entirely.
            p.IFboundaries.nucBoundaries2(ismember(p.IFboundaries.nucBoundaries2.cellID, cellIDsToDelete), :) = [];
            p.IFboundaries.cellBoundaries2(ismember(p.IFboundaries.cellBoundaries2.cellID, cellIDsToDelete), :) = [];
            p.IFquant(ismember(p.IFquant.cellID, cellIDsToDelete),:) = [];
        end
        
        function p = requantCell(p, cellID)
            nucPoly = p.IFboundaries.nucBoundaries2.nucBoundary(p.IFboundaries.nucBoundaries2.cellID == cellID);
            cellPoly = p.IFboundaries.cellBoundaries2.cellBoundary(p.IFboundaries.cellBoundaries2.cellID == cellID);
            tmpBB = d2utils.polyshapeBoundingBox(union(nucPoly, cellPoly)); %Union in case nuclei boundary extends beyond cell boundary
            tmpDapiMask = d2utils.polyvect2mask(tmpBB(3:4), regions(nucPoly), tmpBB(1:2), 'flip', true);
            tmpCellMask = d2utils.polyvect2mask(tmpBB(3:4), regions(cellPoly), tmpBB(1:2), 'flip', true);
            tmpCytoMask = tmpCellMask & ~tmpDapiMask;
            %Update IF quant table
            cellIdx = p.IFquant.cellID == cellID;
            for i = 1:numel(p.channels)
                [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, p.channels{i});
                cellChannelIdx = cellIdx & p.IFquant{:,p.channels{i}};
                p.IFquant{cellChannelIdx,{'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}} = [meanNuc, meanCyto, sumNuc, sumCyto];
            end
        end
        
        function p = requantCellChannel(p, channel, cellIDs)
            for i = 1:numel(cellIDs)
                nucPoly = p.IFboundaries.nucBoundaries2.nucBoundary(p.IFboundaries.nucBoundaries2.cellID == cellID);
                cellPoly = p.IFboundaries.cellBoundaries2.cellBoundary(p.IFboundaries.cellBoundaries2.cellID == cellID);
                tmpBB = d2utils.polyshapeBoundingBox(union(nucPoly, cellPoly)); %Union in case nuclei boundary extends beyond cell boundary
                tmpDapiMask = d2utils.polyvect2mask(tmpBB(3:4), regions(nucPoly), tmpBB(1:2), 'flip', true);
                tmpCellMask = d2utils.polyvect2mask(tmpBB(3:4), regions(cellPoly), tmpBB(1:2), 'flip', true);
                tmpCytoMask = tmpCellMask & ~tmpDapiMask;
                %Update IF quant table
                cellChannelIdx = p.IFquant.cellID == cellID & p.IFquant{:,channel};
                [meanNuc, meanCyto, sumNuc, sumCyto] = p.quantCellInChannel(tmpDapiMask, tmpCytoMask, tmpBB, channel);
                p.IFquant{cellChannelIdx,{'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}} = [meanNuc, meanCyto, sumNuc, sumCyto];
            end
        end
        
        function [meanNuc, meanCyto, sumNuc, sumCyto] = quantCellInChannel(p, dapiMask, cytoMask, rect, channel)
            tmpImg = p.scanObj.getImageRect(channel, rect);
            tmpMasks = p.maskObj.getChannelMasksInRect(rect, channel);
            if ~isempty(tmpMasks)
                maskIDs = unique(tmpMasks.maskID);
                for i = 1:numel(maskIDs)
                    bwmask = poly2mask(tmpMasks{tmpMasks.maskID == maskIDs(i), 'x'}, tmpMasks{tmpMasks.maskID == maskIDs(i), 'y'}, size(tmpDapiMask, 1), size(tmpDapiMask, 2));
                    dapiMask(bwmask) = false;
                    cytoMask(bwmask) = false;
                end
            end
            meanNuc = mean(tmpImg(dapiMask));
            meanCyto = mean(tmpImg(cytoMask));
            sumNuc = sum(tmpImg(dapiMask));
            sumCyto = sum(tmpImg(cytoMask));
        end
        
        function p = makeCentroidList(p, var)
            p.centroidLists = cell(0, numel(p.channels));
            for i = 1:numel(p.channels)
                p.centroidLists{i} = p.sortChannel(p.channels{i}, var);
            end
        end
        
        function p = updateCentroidList(p, channel, var)
            channelIdx = ismember(p.channels, channel);
            p.centroidLists{channelIdx} = p.sortChannel(p.channels{channelIdx}, var);
        end
        
        function p = updateExpressionColors(p, var)
            for i = 1:numel(p.channels)
                p.centroidLists{i}.expression_color = d2utils.expressionToColors(p.centroidLists{i}{:,var}, p.expressionColorPal{p.paletteIdx});
            end
        end
        
        function outTable = sortedChannelIndex(p, channel)
            outTable = varfun(@(x) d2utils.sortIndex(x, 'descend'), p.IFquant(p.IFquant{:,channel}, {'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}));
            outTable.Properties.VariableNames =  {'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'};
        end
        
        function outTable = sortChannel(p, channel, var)
            outTable = sortrows(p.IFquant(p.IFquant{:,channel}, {'cellID', 'x', 'y', var}), var, 'descend', 'MissingPlacement', 'last');
            outTable{:,var} = round(outTable{:,var});
            outTable.expression_color =  d2utils.expressionToColors(outTable{:,var}, p.expressionColorPal{p.paletteIdx});
        end
        
        function outTable = centroidTableInRect(p, channelIdx, rect)
            centroidIdx = p.centroidLists{channelIdx}.x >= rect(1)...
                & p.centroidLists{channelIdx}.x < rect(1) + rect(3) ...
                & p.centroidLists{channelIdx}.y >= rect(2)...
                & p.centroidLists{channelIdx}.y < rect(2) + rect(4);
            outTable = p.centroidLists{channelIdx}(centroidIdx,:);
        end
        
        
        function p = updateAllMasks(p) %This will overwrite previous maskIDs in nuclei. 
            
            maskTable = p.maskObj.masks(p.maskObj.masks.dapi, :);
            maskIDs = unique(maskTable.maskID);
            maskIDs(maskIDs == 0) = [];
            p.nuclei = p.nuclei(~(p.nuclei.nucID == 0),:);
            p.nuclei.maskID(:) = single(0);
            p.nuclei.status(:) = true;
            for i = 1:numel(maskIDs)
                idx = inpolygon(p.nuclei.x, p.nuclei.y,...
                    maskTable{maskTable.maskID == maskIDs(i), 'x'}, maskTable{maskTable.maskID == maskIDs(i), 'y'}) & p.nuclei.status;
                p.nuclei.maskID(idx) = maskIDs(i);
                p.nuclei.status(idx) = false;
            end
            p.addEmptyRows(1000);
        end
        
        function p = addNewMask(p)
            maxCellMask = max(p.maskObj.masksBB{p.maskObj.masksBB.dapi,'maskID'});
            maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxCellMask,{'BB'}}; %Only query nuclei within mask bouding box
            %polyRect = d2utils.boundingCorners2Rect(maskBB);
            
            nucIdx = p.getNucleiInRectIdx(maskBB);
            polyIdx = inpolygon(p.nuclei.x(nucIdx), p.nuclei.y(nucIdx),...
                p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'y'});
            
            nucIdx(nucIdx) = polyIdx; %Only nuclei in polygon remain true
            p.nuclei.maskID(nucIdx) = maxCellMask;
            p.nuclei.status(nucIdx) = false;
            
        end
       
        function p = updateMasksInRect(p, localRect) 
            %Use for removing masks or if we want to change multiple masks before updating nuclei
            %Probably could use some speeding up. 
            masksInRect = p.maskObj.getChannelMasksInRect(localRect, 'dapi');
            maskIDsinRect = unique(masksInRect.maskID);
            [tmpNuclei, nucIdx] = p.getNucleiInRect(localRect); %It might be faster to not subset the nuclei and just run inpolygon on entire nuclei table. 
            
            %Resest status for nuclei
            tmpNuclei.maskID(:) = single(0);
            tmpNuclei.status(:) = true;
            tmpNuclei.status(tmpNuclei.nucID == 0) = false;
            for i = 1:numel(maskIDsinRect)
                idx = inpolygon(tmpNuclei.x, tmpNuclei.y,...
                    masksInRect{masksInRect.maskID == maskIDsinRect(i), 'x'}, masksInRect{masksInRect.maskID == maskIDsinRect(i), 'y'}) & tmpNuclei.status;
                tmpNuclei.maskID(idx) = maskIDsinRect(i);
                tmpNuclei.status(idx) = false;
                
            end
            p.nuclei.maskID(nucIdx) = tmpNuclei.maskID;
            p.nuclei.status(nucIdx) = tmpNuclei.status;
 
        end
        
        function p = removeMasks(p) 
            
            masksToRemove = setdiff(p.nuclei.maskID, p.maskObj.masksBB.maskID(p.maskObj.masksBB.dapi));
            masksToRemove(masksToRemove == 0) = [];
            if ~isempty(masksToRemove)
                nucIdx = ismember(p.nuclei.maskID, masksToRemove);
                p.nuclei.maskID(nucIdx) = single(0);
                p.nuclei.status(nucIdx) = true;
                p.nucleiChanged = true;
            end
            
        end
        
    end
end
