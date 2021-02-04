classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        masks
        masksBB
        thresholds
        radius = 300; 
        spotChannels
        maxDistance = 100;
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
                p.spots = cell2table(cell(0,8),...
                    'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'status', 'maskID', 'channel'});
            else
                fprintf('Loading spot table');
                p.spots = readtable(n.Results.spotsFile,'TextType','string');
                p.spots{:,{'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'maskID'}} =...
                    single(p.spots{:,{'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'maskID'}});
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
                bb = d2utils.polygonBoundingBox(poly); %We could replace this with matlab builting func
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
                tempMaskID = 1;
            else
                tempMaskID = max(p.masks.maskID)+1;
            end
            
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            %maskPoly = [x y];
            corners = d2utils.polygon2BoundingCorners([x,y]);
            
            p.masks = [p.masks; table(repmat(single(tempMaskID), length(x), 1),repmat(channel, length(x), 1),single(x),single(y), 'VariableNames', {'maskID', 'channel', 'x', 'y'})];
            p.masksBB = [p.masksBB; table(repmat(single(tempMaskID), 4, 1),repmat(channel, 4, 1),corners(:,1),corners(:,2), 'VariableNames', {'maskID', 'channel', 'x', 'y'})];
            
            %UPDATE SPOTS
            tempSpots = p.getValidNonmaskedSpotsInRect(channel,localRect);
            idx = inpolygon(tempSpots.x,tempSpots.y,x,y);
            tempSpots{idx,'maskID'} = tempMaskID;
            tempSpots{idx,'status'} = false;
            p.spots(ismember(p.spots.spotID, tempSpots.spotID), :) = tempSpots;
        end
        
        function p = removeMasks(p,maskIDs)
            for i = 1:numel(maskIDs)
                p.masks(ismember(p.masks.maskID,maskIDs(i)),:) = [];
                p.masksBB(ismember(p.masksBB.maskID,maskIDs(i)),:) = [];
            
                %UPDATE SPOTS 
                spotIdx = ismember(p.spots.maskID,maskIDs(i));
                p.spots{spotIdx, 'status'} = true;
                p.spots{spotIdx, 'maskID'} = single(0);
            
            end
            
        end
        
        function p = updateMaskPoly(p,channel,maskID,maskPoly,localRect)
            p.removeMasks(maskID);
            p.addMask(maskPoly,localRect,channel)
        end
        
        function outMasks = getMasksInRect(p,channel, rect)
            
            tempMasksBB = p.masksBB(p.masksBB.channel == channel,:);

            idx = tempMasksBB.x >= rect(1) & tempMasksBB.x < rect(1) + rect(3) ...
                & tempMasksBB.y >= rect(2) & tempMasksBB.y < rect(2) + rect(4);
            
            maskIDtoKeep = unique(tempMasksBB.maskID(idx));
            
            outMasks = p.masks(p.masks.channel == channel,:);
            outMasks = outMasks(ismember(outMasks.maskID, maskIDtoKeep),:);
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
        
        function p = assignSpotsToCells2(p,cellTable)
            
            
            [idx, dist] = knnsearch([cellTable.x cellTable.y], [p.spots.x p.spots.y], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestCellID = single(cellTable.cellID(idx));
            
            p.spots{dist > p.maxDistance,'nearestCellID'} = single(0);

        end
        
        function p = assignSpotsInRect(p, channel, cellObject, rect) %Maybe useful for reassigning spots after add/removing cells 
            
            spotsInRect = getSpotsInRect(p.spots, channel,rect);
            cellsInRect = cellObject.getCellsInRect(rect);
            [idx_cell, dist] = knnsearch([cellsInRect.x cellsInRect.y], [spotsInRect.x spotsInRect.y], 'K', 1, 'Distance', 'euclidean');
            
            spotsInRect.nearestCellID = single(cellsInRect.cellID(idx_cell));
            
            spotsInRect{dist > p.maxDistance,'nearestCellID'} = single(0);
            
            p.spots(ismemeber(p.spots.spotID, spotsInRect.spotID), :) = spotsInRect;

        end
        
%         function p = updateSpotsForChannel(p, channel) 
            
           %add code

%         end
        
%         function p = updateSpotsInRect(p, channel, rect) 
            
           %add code

%         end
        
        function intensities = getIntensities(p,channel)
            tempTable = p.spots(p.spots.channel == channel,:);
            intensities = tempTable.intensity;
        end
        
        function [outSpots,idx] = getAllSpotsInRect(p,rect) %rect specified as [x y nrows ncols]
            outSpots = p.spots;

            idx = outSpots.x >= rect(1) & outSpots.x < rect(1) + rect(3) ...
                & outSpots.y >= rect(2) & outSpots.y < rect(2) + rect(4);
            
            outSpots = outSpots(idx,:);
        end
        
        function outSpots = getValidNonmaskedSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]
            outSpots = p.spots(p.spots.channel == channel,:);
            outSpots = outSpots(outSpots.status,:);
            outSpots = outSpots(outSpots.maskID == 0,:);

            idx = outSpots.x >= rect(1) & outSpots.x < rect(1) + rect(3) ...
                & outSpots.y >= rect(2) & outSpots.y < rect(2) + rect(4);
            
            outSpots = outSpots(idx,:);
        end
        
        function outSpots = getSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]
            outSpots = p.spots(p.spots.channel == channel,:);
            %outSpots = outSpots(outSpots.status,:);

            idx = outSpots.x >= rect(1) & outSpots.x < rect(1) + rect(3) ...
                & outSpots.y >= rect(2) & outSpots.y < rect(2) + rect(4);
            
            outSpots = outSpots(idx,:);
        end
        
        function idx = getSpotsInRectIndex(p,channel,rect) %BE - I don't know what this is for. 
            chanIdx = p.spots.channel == channel;
            outSpots = p.spots(p.spots.channel == channel,:);
            %outSpots = outSpots(outSpots.status,:);

            idx = outSpots.x >= rect(1) & outSpots.x < rect(1) + rect(3) ...
                & outSpots.y >= rect(2) & outSpots.y < rect(2) + rect(4);
            
            chanIdx(find(chanIdx)) = idx;
            
            idx = chanIdx;
            
        end
        
        function applyThreshold(p,channel,threshold)
            
            
            spotIdx = p.spots.channel == channel & p.spots.intensity < threshold;
            
            p.spots.status(spotIdx) = false;
            p.spots.status(~spotIdx) = true;

%             idx = p.spots.channel == channel;
%             
%             
%             
%             intensities = p.spots.intensity(idx);
%             spotIdx = intensities > threshold;
            
        end
        
        
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
                fprintf('Channel: %s\n',currChannel);
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
            
            p.spots = table(spotID, single(x), single(y), intensity, nearestCellID, status, maskID, channel,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'status', 'maskID', 'channel'});
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