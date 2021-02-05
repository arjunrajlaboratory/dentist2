classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        masks
        masksBB
        thresholds
        radius = 300; 
        spotChannels
        maxDistance = 100; %Should this be channel specific? 
        theFilter
        percentileToKeep = 98;
        
    end
    
    methods
        
        function p = spotTable(scanObject, varargin)
            n = inputParser;
            n.addRequired('scanObject');
            n.addParameter('spotsFile', '', @ischar); 
            n.addParameter('masksFile', '', @ischar); 
            n.addParameter('thresholdsFile', '', @ischar);

            n.parse(scanObject, varargin{:});
            
            channels = n.Results.scanObject.channels;
            p.spotChannels = channels(~ismember(channels,["dapi","trans"]));
            
            if isempty(n.Results.spotsFile)
                fprintf('New spots table\n');
                p.spots = cell2table(cell(0,9),...
                    'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'status', 'maskID', 'channel', 'distanceToCell'});
            else
                fprintf('Loading spot table');
                p.spots = readtable(n.Results.spotsFile,'TextType','string');
                p.spots{:,{'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'maskID', 'distanceToCell'}} =...
                    single(p.spots{:,{'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'maskID', 'distanceToCell'}});
                p.spots.status = logical(p.spots.status);
            end
            
            if isempty(n.Results.masksFile)
                fprintf('New masks table\n');
                p.masks = cell2table(cell(0,4), 'VariableNames', {'maskID', 'channel', 'x', 'y'});
                p.masksBB = cell2table(cell(0,4), 'VariableNames', {'maskID', 'channel', 'x', 'y'});
            else
                fprintf('Loading masks table');
                p.masks = readtable(n.Results.masksFile,'TextType','string');
                p.masks{:,{'maskID', 'x', 'y'}} =...
                    single(p.masks{:,{'maskID', 'x', 'y'}});
                p.allMasks2Corners();
            end
            
            if isempty(n.Results.thresholdsFile)
                fprintf('New thresholds table\n');
                p.thresholds = cell2table(cell(0,numel(p.spotChannels)), 'VariableNames', p.spotChannels);
            else
                fprintf('Loading thresholds table');
                p.thresholds = readtable(n.Results.thresholdsFile,'TextType','string');
            end
        end
        
        
        function p = maskAllPoints(p) 
            for i = 1:height(p.masks)
                poly = p.masks.polygon(i);
                poly = poly{1};
                channel = p.masks.channel(i);
                bb = d2utils.polygonBoundingBox(poly);
                spotsIdx = p.getSpotsInRectIndex(channel,bb);
                p.spots.maskID(spotsIdx) = p.masks.maskID(i);
            end
        end
        
        function p = maskAllPoints2(p)
            %create logical array with all spots
            %Create mask array
            %multiple
            %Update spots table

         end
        
        function p = addMask(p,maskPoly,localRect,channel)
            % Mask table should be:
            % maskID
            % handle
            % channel (can also be "all" to apply to all)
            % DO NOT NEED SINCE HAVE HANDLE isValid % 0 when masks are deleted
            % polygon in global coordinates
            % Should have a bounding box function to get the bounding box
            % of the polygon to test for intersection with current viewRect
            if isempty(p.masks)
                tempMaskID = single(1);
            else
                tempMaskID = single(max(p.masks.maskID)+1);
            end
            
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            %maskPoly = [x y];
            corners = d2utils.polygon2BoundingCorners([x,y]);
            
            p.masks = [p.masks; table(repmat(tempMaskID, length(x), 1),repmat(channel, length(x), 1),single(x),single(y), 'VariableNames', {'maskID', 'channel', 'x', 'y'})];
            p.masksBB = [p.masksBB; table(repmat(tempMaskID, 4, 1),repmat(channel, 4, 1),corners(:,1),corners(:,2), 'VariableNames', {'maskID', 'channel', 'x', 'y'})];
                        %UPDATE CELLS - inpolygon seems to fast enough for many millions of cells. 

            %UPDATE SPOTS - inpolygon seems to fast enough for many millions of spots (scales O(n)).
            %If it gets too slow, could try selecting cells in localRect first.
%             tempSpots = p.getValidNonmaskedSpotsInRect(channel,localRect);
%             tempIdx = inpolygon(tempSpots.x,tempSpots.y,x,y);
%             idx = ismember(p.spots.spotID, tempIdx.spotID);
            idx = p.spots.status & p.spots.channel == channel ...
                & inpolygon(p.spots.x,p.spots.y,x,y);
            p.spots.maskID(idx) = tempMaskID;
            p.spots.status(idx) = false;
        end
        
        function p = removeMasks(p,maskIDs, channel)
             p.masks(ismember(p.masks.maskID,maskIDs),:) = [];
             p.masksBB(ismember(p.masksBB.maskID,maskIDs),:) = [];
             
             p.spots(ismember(p.spots.maskID,maskIDs), 'maskID') = single(0);
             %Need to check threshold then update status.
             p.updateSpotStatus(channel)
        end
        
        function p = updateMaskPoly(p,channel,maskID,maskPoly,localRect)
            p.removeMasks(maskID);
            p.addMask(maskPoly,localRect,channel)
        end
        
        function outMasks = getMasksInRect(p,channel, rect)
            
            %tempMasksBB = p.masksBB(p.masksBB.channel == channel,:);

            idx = p.masksBB.channel == channel ...
                & p.masksBB.x >= rect(1) & p.masksBB.x < rect(1) + rect(3) ...
                & p.masksBB.y >= rect(2) & p.masksBB.y < rect(2) + rect(4);
            
            maskIDtoKeep = unique(tempMasksBB.maskID(idx));
            
            %outMasks = p.masks(p.masks.channel == channel,:);
            outMasks = p.masks(p.masks.channel == channel & ismember(p.masks.maskID, maskIDtoKeep) ,:);
        end
        
        function p = allMasks2Corners(p)
            p.masksBB = cell2table(cell(0,4), 'VariableNames', {'maskID', 'channel', 'x', 'y'});
            uniqueMaskIDs = unique(p.masks.maskID);
            for i = 1:numel(uniqueMaskIDs)
                tempPoly = p.masks{p.masks.maskID == uniqueMaskIDs(i), ['x', 'y']};
                tempChannel = p.masks{p.masks.maskID == uniqueMaskIDs(i), 'channel'};
                tempCorners = d2utils.polygon2BoundingCorners(tempPoly);
                p.masksBB = [p.masksBB; {repmat(uniqueMaskIDs(i), 4, 1),repmat(tempChannel, 4, 1),tempCorners(:,1),tempCorners(:,2)}];
            end
        end
        
%         function p = assignSpotsToCells2(p,cellTable) %Version in case we
%         want to have unique cellIDs for different channels
%             
%             channels = p.spotChannels;
%             for i = 1:numel(channels)
%                 tempSpots = p.spots(p.spots.channel == channel,:);
%                 tempCells = cellTable(cellTable.channel == channel,:);
%                 [idx, dist] = knnsearch([tempCells.x tempCells.y], [tempSpots.x tempSpots.y], 'K', 1, 'Distance', 'euclidean');
%                 
%                 p.spots{ismember(p.spots.spotID, tempSpots.spotID), 'nearestCellID'} = single(tempCells.cellID(idx));
%                 p.spots{ismember(p.spots.spotID, tempSpots.spotID), 'distanceToCell'} = single(dist);
%                 %farCells = tempSpots.spotID(dist > p.maxDistance);
%                 %p.spots{ismember(p.spots.spotID, farCells), 'nearestCellID'} = single(0);
%                 
%             end
%         end
        
        function p = assignSpotsToCells2(p,cellTable)
            
            [idx, dist] = knnsearch([cellTable.x cellTable.y], [p.spots.x p.spots.y], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestCellID = cellTable.cellID(idx);
            p.spots.distanceToCell = single(dist);
            
        end
        
        function p = assignSpotsInRect(p, channel, rect, cellObject) %Maybe useful for reassigning spots after add/removing cells 
            
            spotIdx = p.getValidSpotsInRectIndex(channel,rect);
            cellsNearRect = cellObject.getCellsNearRect(rect, p.maxDistance);
            [cellIdx, dist] = knnsearch([cellsNearRect.x cellsNearRect.y], [p.spots.x(spotIdx) p.spots.y(spotIdx)], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestCellID(spotIdx) = cellsNearRect.cellID(cellIdx);
            p.spots.distanceToCell(spotIdx) = single(dist);
            
            %spotsInRect.nearestCellID = single(cellsInRect.cellID(cellIdx));
            %spotsInRect.distanceToCell = single(dist);
            
            %spotsInRect{dist > p.maxDistance,'nearestCellID'} = single(0);
            
            %p.spots(ismemeber(p.spots.spotID, spotsInRect.spotID), :) = spotsInRect;

        end
        
        function intensities = getIntensities(p,channel)
            intensities = p.spots{p.spots.channel == channel,'intensity'};
        end
        
        function intensityThresh = getIntensityThreshold(p,channel)
            intensityThresh = p.thresholds{p.spots.channel == channel};
        end
        
        function [outSpots,idx] = getAllSpotsInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        function outSpots = getValidNonmaskedSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]
            %To speed things up, creating one index instead of multiple assignments.  
%             outSpots = p.spots(p.spots.channel == channel,:);
%             outSpots = outSpots(outSpots.status,:);
%             outSpots = outSpots(outSpots.maskID == 0,:);

            idx = p.spots.channel == channel & p.spots.status ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        function outSpots = getSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]
            %outSpots = p.spots(p.spots.channel == channel,:);
            %outSpots = outSpots(outSpots.status,:);

            idx = p.spots.channel == channel ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        function idx = getValidSpotsInRectIndex(p,channel,rect) 
            
            idx = p.spots.channel == channel & p.spots.status ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
        end
       
        function p = updateSpotStatus(p,channel)
            threshold = p.thresholds{p.spotChannels == channel};
            spotIdx = p.spots.channel == channel & p.spots.intensity >= threshold ...
                & p.spots.distanceToCell <= p.maxDistance & p.spots.maskID == 0;
            
            p.spots.status(spotIdx) = true;
            p.spots.status(~spotIdx) = false;
                    
        end
        
        function p = updateAllSpotStatus(p)
            for i = 1:numel(p.spotChannels)
                p.updateSpotStatus(p.spotChannels{i})
            end
        end
        
%         function applyIntensityThreshold(p,channel,threshold)
%             
%             
%             spotIdx = p.spots.channel == channel & p.spots.intensity < threshold;
%             
%             p.spots.status(spotIdx) = false;
%             p.spots.status(~spotIdx) = true;
%             
%             %UPDATE p.thresholds
%             
%         end
%         
%         function applyDistThreshold(p,channel,threshold) %Should this be channel specific? 
%             
%             
%             spotIdx = p.spots.channel == channel & p.spots.distanceToCell < threshold;
%             
%             p.spots.status(spotIdx) = false;
%             p.spots.status(~spotIdx) = true;
%             
%             %UPDATE p.maxDistance
%             
%         end
        
        
        function p = findSpots(p, scanObject) %BE should we make a faster version that uses blockproc on stitched data? 

            if isempty(p.theFilter) == 2
                filt = -fspecial('log',20,2);
                p.theFilter = filt;
            else
                filt = p.theFilter;
            end
            tilesTable = scanObject.tilesTable;
            scanMatrix = scanObject.scanMatrix;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            
            x = [];
            y = [];
            intensity = [];
            channel = [];
            reader = bfGetReader(scanObject.scanFile);
            
            for i = 1:numel(p.spotChannels)
                currChannel = p.spotChannels(i);
                fprintf('Finding %s spots\n',currChannel);
                iPlane = reader.getIndex(0, i - 1, 0) + 1;
                channelCount = 0;
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    [tempX, tempY, tempIntensity] = d2utils.findSpotsInImage(tmpPlane, p.percentileToKeep, 'filter', filt);
                    %Adjust local to global coordinates
                    tempX = tempX + tilesTable.left(tiles(ii)) - 1;
                    tempY = tempY + tilesTable.top(tiles(ii)) - 1; 
                    
                    x = [x ; tempX];
                    y = [y ; tempY];
                    intensity = [intensity ; tempIntensity];
                    channelCount = channelCount + length(tempX);
                end
                channel = [channel ; repmat(currChannel,channelCount,1)];
            end
            reader.close()
            
            spotID = single((1:length(x)))';
            nearestCellID = single(zeros(length(x),1));
            maskID = single(zeros(length(x),1));
            status = true(length(x),1);
            dist =  single(zeros(length(x),1));
            p.spots = table(spotID, single(x), single(y), intensity, nearestCellID, status, maskID, channel, dist,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'status', 'maskID', 'channel', 'distanceToCell'});
        end
        
        function outTable = getSpotsPerCell(p)
            %INSERT CODE
        end
        
        %Set functions
        function set.theFilter(p,filter)
            p.theFilter = filter;
        end
        
        function set.percentileToKeep(p,n)
            p.percentileToKeep = n;
        end
        
        function [] = saveTables(p, varargin)
            n = inputParser;
            n.addParameter('spotsFile', 'spots.csv', @(x) assert((ischar(x) & endsWith(x, '.csv')),...
                'Specify .csv filename to save spots. e.g "spots.csv"')); 
            n.addParameter('masksFile', 'masks.csv', @(x) assert((ischar(x) & endsWith(x, '.csv')),...
                'Specify .csv filename to save masks. e.g "masks.csv"')); 
            n.addParameter('thresholdsFile', 'thresholds.csv', @(x) assert((ischar(x) & endsWith(x, '.csv')),...
                'Specify .csv filename to save thresholds. e.g "thresholds.csv"'));

            n.parse(varargin{:});
            
            writetable(p.spots, n.Results.spotsFile);
           % writetable(p.masks, n.Results.masksFile);
           % writetable(p.thresholds, n.Results.thresholdsFile);
        end
        
        
    end
    
end 