classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        centroidLists %Not sure if this should go into the View
        
        thresholds
        spotChannels
        maxDistance = 100;  
        theFilter
        percentileToKeep = 98;
        intensitiesToPlot
        
        scanObj
        maskObj 
        nucleiObj 
        
    end
    
    methods
        
        function p = spotTable(scanObject, maskObject, nucleiObject, varargin)
            n = inputParser;
            n.addRequired('scanObject');
            n.addRequired('maskObject');
            n.addRequired('nucleiObject');
            n.addParameter('spotsFile', '', @ischar); 
            n.addParameter('thresholdsFile', '', @ischar);

            n.parse(scanObject, maskObject, nucleiObject, varargin{:});
            
            channels = n.Results.scanObject.channels;
            p.spotChannels = channels(~ismember(channels,{'dapi','trans'}));
            p.scanObj = scanObject;
            p.maskObj = n.Results.maskObject;
            p.nucleiObj = n.Results.nucleiObject;
            
            if isempty(n.Results.spotsFile)
                fprintf('New spots table\n');
                p.spots = cell2table(cell(0,9),...
                    'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'});
            else
                fprintf('Loading spot table\n');
                tmpSpots = readtable(n.Results.spotsFile,'TextType','string');
                p.spots = convertvars(tmpSpots, {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'maskID', 'distanceToNuc'}, 'single');
                p.spots.status = logical(p.spots.status);
            end
                        
            if isempty(n.Results.thresholdsFile)
                fprintf('New thresholds table\n');
                p.thresholds = cell(0,numel(p.spotChannels));
            else
                fprintf('Loading thresholds table\n'); %Need to check this
                tmpThreshold = readtable(n.Results.thresholdsFile,'TextType','string');
                p.thresholds = table2cell(tmpThreshold);
            end
        end
        
        
        function p = assignSpotsToNuclei(p)
            validNucleiIdx = p.nucleiObj.nuclei.status;
            
            [idx, dist] = knnsearch([p.nucleiObj.nuclei.x(validNucleiIdx), p.nucleiObj.nuclei.y(validNucleiIdx)], [p.spots.x, p.spots.y], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestNucID = p.nucleiObj.nuclei.nucID(idx);
            p.spots.distanceToNuc = single(dist);
            p.spots.colors = p.nucleiObj.nuclei.colors(idx, :);
            
        end
        
        function p = assignSpotsInRect(p, rect) %Maybe useful for reassigning spots after add/removing cells 
            
            [spotsInRect, spotIdx] = p.getAllSpotsInRect(rect);
            nucleiNearRect = p.nucleiObj.getNucleiNearRect(rect, p.maxDistance);
            [nucIdx, dist] = knnsearch([nucleiNearRect.x nucleiNearRect.y], [spotsInRect.x spotsInRect.y], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestNucID(spotIdx) = nucleiNearRect.nucID(nucIdx);
            p.spots.distanceToNuc(spotIdx) = single(dist);
            p.spots.colors(spotIdx, :) = nucleiNearRect.colors(nucIdx, :);
        end
        
        function [outSpots,idx] = getAllSpotsInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        function [outSpots, idx] = getValidSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]

            idx = p.spots.channel == channel & p.spots.status ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        function [outSpots,idx] = getSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]
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
        
        function idx = getSpotsInRectIndex(p,channel,rect)

            idx = p.spots.channel == channel ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
        end
                
        function p = updateAllMasks(p)  %Note, this will overwrite previous maskIDs in spotTable. 
            
            for i = 1:numel(p.spotChannels)
                channelIdx = p.maskObj.masks(:, p.spotChannels(i));
                maskTable = p.maskObj.masks(channelIdx,:);
                maskIDs = unique(maskTable.maskID);
                maskIDs(maskIDs == 0) = [];
                spotIdx = p.spots.channel == channel;
                p.spots.maskID(spotIdx) =  single(0);
                p.updateSpotStatus(p.spotChannels(i))
                
                validSpotIdx = p.spots.channel == channel & p.spots.status;
                spotTable = p.spots(validSpotIdx,:);
                
                for ii = 1:numel(maskIDs)
                    idx = inpolygon(spotTable.x, spotTable.y,...
                        maskTable{maskTable.maskID == maskIDs(ii), 'x'}, maskTable{maskTable.maskID == maskIDs(ii), 'y'}) & spotTable.status;
                    spotTable.maskID(idx) = maskIDs(ii); 
                    spotTable.status(idx) = false; 
                end
                p.spots.maskID(validSpotIdx) = spotTable.maskID;
                p.spots.status(validSpotIdx) = spotTable.status;
            end 
            
        end
        
        function p = addNewMask(p, channel)
            
            maxSpotMask = max(p.maskObj.masks{p.maskObj.masks{:,channel},'maskID'});
            maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxSpotMask,{'BB'}}; %Only query spots within mask bouding box
            %polyRect = d2utils.boundingCorners2Rect(maskBB);
            spotIdx = p.getSpotsInRectIndex(channel, maskBB);
            idx = inpolygon(p.spots.x(spotIdx), p.spots.y(spotIdx),...
                p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'y'});
            
            spotIdx(spotIdx) = idx; %Only spots in polygon remain true
                
            p.spots.maskID(spotIdx) = maxSpotMask;
            p.spots.status(spotIdx) = false;
            
        end
        
        function p = updateMasksInRect(p, channel, localRect) 
            %Use for removing masks or if we want to change multiple masks
            %before updating spots
            %Probably could use some speeding up. 
            masksInRect = p.maskObj.getChannelMasksInRect(localRect, channel);
            maskIDsinRect = unique(masksInRect.maskID);
            [tmpSpots, spotIdx] = p.getSpotsInRect(channel, localRect);            
            
            %Resest status for nuclei
            tmpSpots.maskID(:) = single(0);
            tmpSpots.status(:) = true;
            for i = 1:numel(maskIDsinRect)
                idxPoly = inpolygon(tmpSpots.x, tmpSpots.y,...
                    masksInRect{masksInRect.maskID == maskIDsinRect(i), 'x'}, masksInRect{masksInRect.maskID == maskIDsinRect(i), 'y'}) & tmpSpots.status;
                
                tmpSpots.maskID(idxPoly) = maskIDsinRect(i);
                tmpSpots.status(idxPoly) = false;
                
            end
            
            p.spots.maskID(spotIdx) = tmpSpots.maskID;
            p.updateSpotStatus(channel);
        end
        
        function p = removeMasks(p, channel) 
            masksToRemove = setdiff(p.spots.maskID, p.maskObj.masks.maskID(p.maskObj.masks{:,channel}));
            p.spots.maskID(ismember(p.spots.maskID, masksToRemove)) = single(0);
            p.updateSpotStatus(channel);
        end
        
        function p = removeMasks2(p, channel, rect)
            %Not sure if it'll be faster to first get spots and masks in
            %rect. 
            [spotsInRect, spotIdx] = p.getSpotsInRect(channel, rect);
            maskIDsInRect = p.maskObj.getChannelMaskIDsInRect(rect, channel);
            
            goodLabelIdx = ismember(spotsInRect.maskID, [maskIDsInRect; 0]);  %add zero in case maskTable is full
            spotIdx(goodLabelIdx) = false;
            p.spots.maskID(spotIdx) = single(0);
            p.updateSpotStatus(channel);
        end
        
        function p = setThreshold(p, channel, value)
            p.thresholds{ismember(p.spotChannels, channel)} = value;
            %Update spot status
            p.updateSpotStatus(channel);
            p.updateCentroidList(channel);
        end
       
        function p = updateSpotStatus(p,channel)
            threshold = p.thresholds{ismember(p.spotChannels, channel)};
            channelIdx = ismember(p.spots.channel, channel);
            spotIdx = p.spots.intensity(channelIdx) >= threshold ...
                & p.spots.distanceToNuc(channelIdx) <= p.maxDistance & p.spots.maskID(channelIdx) == 0 ...
                & ismember(p.spots.nearestNucID(channelIdx), p.nucleiObj.nuclei.nucID(p.nucleiObj.nuclei.status));%spots near masked cells will be set to false. 
            
            p.spots.status(channelIdx) = spotIdx;
            %p.spots.status(~spotIdx) = false;
                    
        end
        
        function p = updateAllSpotStatus(p)
            for i = 1:numel(p.spotChannels)
                p.updateSpotStatus(p.spotChannels{i});
            end
        end
       
        function p = findSpots(p) %BE should we make a faster version that uses blockproc on stitched data? 

            if isempty(p.theFilter) == 2
                filt = -fspecial('log',20,2);
                p.theFilter = filt;
            else
                filt = p.theFilter;
            end
            tilesTable = p.scanObj.tilesTable;
            scanMatrix = p.scanObj.scanMatrix;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            
            x = [];
            y = [];
            intensity = [];
            channel = [];
            reader = bfGetReader(p.scanObj.scanFile);
            
            for i = 1:numel(p.spotChannels)
                currChannel = p.spotChannels{i};
                fprintf('Finding %s spots\n',currChannel);
                iPlane = reader.getIndex(0, i - 1, 0) + 1;
                channelCount = 0;
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    [tempX, tempY, tempIntensity] = d2utils.findSpotsInImage(tmpPlane, p.percentileToKeep, 'filter', filt);
                    %Adjust local to global coordinates
                    tempX = tempX + tilesTable.left(tiles(ii));
                    tempY = tempY + tilesTable.top(tiles(ii)); 
                    
                    x = [x ; tempX];
                    y = [y ; tempY];
                    intensity = [intensity ; tempIntensity];
                    channelCount = channelCount + length(tempX);
                end
                channel = [channel ; repmat(string(currChannel),channelCount,1)]; %somewhat less memory with string array vs cell array
            end
            reader.close()
            
            spotID = single((1:length(x)))';
            nearestNucID = single(zeros(length(x),1));
            maskID = single(zeros(length(x),1));
            status = true(length(x),1);
            dist =  single(zeros(length(x),1));
            p.spots = table(spotID, single(x), single(y), intensity, nearestNucID, status, maskID, channel, dist,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'});
        end
        
        function p = findSpots2(p) 
            %Consider method that first stitches the filtered images then
            %finds regional maxima to avoid artifacts of spots in tile
            %overlap. c
        end
        
        function intensities = getIntensities(p, channel)
            intensities = p.spots{p.spots.channel == channel & p.spots.distanceToNuc <= p.maxDistance, {'intensity'}}; %May want to exclude spots not assigned to nuclei
        end
        
        function p = makeIntensitiesToPlot(p)
            p.intensitiesToPlot = cell(0, numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.intensitiesToPlot{i} = sort(uint16(p.getIntensities(p.spotChannels{i})));
            end
        end
        
        function p = defaultThresholds(p) %Need to update this with some better heuristic
            for i = 1:numel(p.spotChannels)
                p.thresholds{i} = round(mean(p.getIntensities(p.spotChannels{i})));
            end
            p.updateAllSpotStatus();
            p.makeCentroidList();
        end
        
        function p = makeCentroidList(p)
            p.centroidLists = cell(0, numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.centroidLists{i} = sortrows(p.tabulateChannel(p.spotChannels{i}), 'GroupCount', 'descend');
                p.centroidLists{i}.expression_color = d2utils.expressionToColors2(p.centroidLists{i}.GroupCount);
            end
        end
        
        function p = updateCentroidList(p, channel)
            channelIdx = ismember(p.spotChannels, channel);
            p.centroidLists{channelIdx}...
                = sortrows(p.tabulateChannel(channel), 'GroupCount', 'descend');
            p.centroidLists{channelIdx}.expression_color = d2utils.expressionToColors2(p.centroidLists{channelIdx}.GroupCount);
        end
        
        function outTable = centroidTableInRect(p, channelIdx, rect)
%             channelIdx = ismember(p.spotChannels, channel);
            centroidIdx = p.centroidLists{channelIdx}.x >= rect(1)...
                & p.centroidLists{channelIdx}.x < rect(1) + rect(3) ...
                & p.centroidLists{channelIdx}.y >= rect(2)...
                & p.centroidLists{channelIdx}.y < rect(2) + rect(4);
            outTable = p.centroidLists{channelIdx}(centroidIdx,:);
        end
        
        function outTable = tabulateChannel(p, channel)
            idx = ismember(p.spots.channel, channel) & p.spots.status;
            tmpSpots = groupsummary(p.spots(idx, :), 'nearestNucID');
            outTable = outerjoin(p.nucleiObj.nuclei(p.nucleiObj.nuclei.status,{'nucID', 'x', 'y', 'colors'}), tmpSpots, 'Type', 'left', 'LeftKeys', 'nucID', 'RightKeys', 'nearestNucID');
            outTable{isnan(outTable.GroupCount), 'GroupCount'} = single(0);
        end
        
        function outTable = tabulateAllChannels(p)
            tmpSpots = groupsummary(p.spots, {'channel', 'nearestNucID'});
            outTable = outerjoin(p.nucleiObj.nuclei(p.nucleiObj.nuclei.status,{'nucID', 'x', 'y'}), tmpSpots, 'Type', 'left', 'LeftKeys', 'nucID', 'RightKeys', 'nearestNucID');
            outTable.nearestNucID = [];
            outTable{isnan(outTable.GroupCount), 'GroupCount'} = single(0);
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
            n.addParameter('thresholdsFile', 'thresholds.csv', @(x) assert((ischar(x) & endsWith(x, '.csv')),...
                'Specify .csv filename to save thresholds. e.g "thresholds.csv"'));

            n.parse(varargin{:});
            
            writetable(p.spots, n.Results.spotsFile);
           % writetable(p.thresholds, n.Results.thresholdsFile);
        end
        
        
    end
    
end 