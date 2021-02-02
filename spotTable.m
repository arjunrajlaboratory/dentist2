classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        masks
        thresholds
%         radius = 300; 
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
            end
            
            if isempty(n.Results.masksFile)
                fprintf('New masks table\n');
                p.masks = cell2table(cell(0,3), 'VariableNames', {'maskID', 'channel', 'polygon'});
            else
                fprintf('Loading masks table');
                p.spots = readtable(n.Results.masksFile,'TextType','string');
            end
            
            if isempty(n.Results.thresholdsFile)
                fprintf('New thresholds table\n');
                p.thresholds = cell2table(cell(0,numel(p.spotChannels)), 'VariableNames', p.spotChannels);
            else
                fprintf('Loading thresholds table');
                p.spots = readtable(n.Results.thresholdsFile,'TextType','string');
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
            maskPoly = [x y];
            
            p.masks = [p.masks; {tempMaskID,channel,{maskPoly},}];
            
            % INSERT CODE HERE TO UPDATE SPOTS - BE: Keep seperate? 
            
        end
        
        function p = removeMasks(p,maskIDs)
            idx = ismember(p.masks.maskID,maskIDs);
            % NEED TO save which ones got deleted for updating.
            p.masks(idx,:) = [];
            
            % INSERT CODE HERE TO UPDATE - BE: Keep seperate? 
        end
        
        function p = updateMaskPoly(p,maskID,maskPoly,localRect)
            idx = find(p.masks.maskID == maskID);
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            maskPoly = [x y];
            p.masks.polygon(idx) = {maskPoly};
            % INSERT CODE HERE TO UPDATE
        end
        
%         function p = assignSpotsToCells(p,cellTable)
%             
%             rect = zeros(1,4);
%             rect(3) = 2*p.radius+1; % Size of rectangle will always be the same
%             rect(4) = 2*p.radius+1;
%             for i = 1:height(p.spots)
%                 if mod(i,100) == 1
%                     fprintf('%d of %d spots\n',i,height(p.spots));
%                 end
%                 rect(1) = p.spots.x(i)-p.radius;
%                 rect(2) = p.spots.y(i)-p.radius;
%                 cells = cellTable.getCellsInRect(rect);
%                 if isempty(cells)
%                     p.spots.nearestCellID(i) = 0; % Should we set this to nan?
%                 else
% %                     if height(cells) == 1 % If there's just one nearby cell, just handle it
% %                         p.spots.nearestCellID(i) = cells.cellID;
% %                     else
%                         distMatrix = [(cells.x-p.spots.x(i))';(cells.y-p.spots.y(i))'];
%                         nm = vecnorm(distMatrix);
%                         [~,idx] = sort(nm);
%                         if nm(idx(1))<p.radius
%                             p.spots.nearestCellID(i) = cells.cellID(idx(1));
%                         end
% %                     end
%                 end
%             end
%                 
%         end
        
        function p = assignSpotsToCells2(p,cellTable)
            
            
            [idx, dist] = dsearchn([cellTable.x cellTable.y], [p.spots.x p.spots.y]);
            
            p.spots.nearestCellID = cellTable.cellID(idx);
            
            p.spots{dist > p.maxDistance,'nearestCellID'} = NaN;

        end
        
        function p = assignSpotsInRect(p, channel, cellTable, rect) %Maybe useful for reassigning spots after add/removing cells 
            
            spotsInRect = getSpotsInRect(p.spots, channel,rect);
            cellsInRect = getCellsInRect(cellTable,rect);
            [idx, dist] = dsearchn([cellsInRect.x cellsInRect.y], [spotsInRect.x spotsInRect.y]);
            
            spotsInRect.nearestCellID = cellsInRect.cellID(idx);
            
            spotsInRect{dist > p.maxDistance,'nearestCellID'} = NaN;
            
            p.spots(ismember(p.spots.spotID, spotsInRect.spotID), :) = spotsInRect;

        end
        
        function p = updateSpotsForChannel(p, channel) 
            
           %add code

        end
        
        function p = updateSpotsInRect(p, channel, rect) 
            
           %add code

        end
        
        function intensities = getIntensities(p,channel)
            tempTable = p.spots(p.spots.channel == channel,:);
            intensities = tempTable.intensity;
        end
        
        function [outSpots,idx] = getAllSpotsInRect(p,rect)
            outSpots = p.spots;

            idx = outSpots.x >= rect(1) & outSpots.x < rect(1) + rect(3) ...
                & outSpots.y >= rect(2) & outSpots.y < rect(2) + rect(4);
            
            outSpots = outSpots(idx,:);
        end
        
        function outSpots = getValidNonmaskedSpotsInRect(p,channel,rect)
            outSpots = p.spots(p.spots.channel == channel,:);
            outSpots = outSpots(outSpots.status,:);
            outSpots = outSpots(outSpots.maskID == 0,:);

            idx = outSpots.x >= rect(1) & outSpots.x < rect(1) + rect(3) ...
                & outSpots.y >= rect(2) & outSpots.y < rect(2) + rect(4);
            
            outSpots = outSpots(idx,:);
            
        end
        
        function outSpots = getSpotsInRect(p,channel,rect)
            outSpots = p.spots(p.spots.channel == channel,:);
            %outSpots = outSpots(outSpots.status,:);

            idx = outSpots.x >= rect(1) & outSpots.x < rect(1) + rect(3) ...
                & outSpots.y >= rect(2) & outSpots.y < rect(2) + rect(4);
            
            outSpots = outSpots(idx,:);
        end
        
        function idx = getSpotsInRectIndex(p,channel,rect)
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
        
        
        function p = findSpots(p, scanObject)

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
            
            spotID = (1:length(x))';
            nearestCellID = nan(length(x),1);
            maskID = zeros(length(x),1);
            status = true(length(x),1);
            
            p.spots = table(spotID, x, y, intensity, nearestCellID, status, maskID, channel,...
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
            writetable(p.masks, n.Results.masksFile);
            writetable(p.thresholds, n.Results.thresholdsFile);
        end
        
        
    end
    
end 