classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        centroidLists %Not sure if this should go into the View
        
        thresholds
        spotChannels
        maxDistance = 100;  
        theFilter
        percentileToKeep = 98;
        spotsIntensitiesWithMasked
        spotsIntensitiesNoMasked
        expressionColorPal = {'BuYlRd', 'YlOrRd', 'GrBu', 'BuGnYlRd'}
        paletteIdx = 1;
        
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
            %At the moment, assigning spots to nearest nucleus, even if that nucleus is masked (status = false)
            %If we want to assign spots to only valid nucleu, then change
            %below to validNuclei = p.nucleiObj.nuclei(p.nucleiObj.nuclei.status, :);

            validNuclei = p.nucleiObj.nuclei(~p.nucleiObj.nuclei.nucID == 0,:); %There may be a more efficient way fo doing this using only index arrays and not creating a new table. 
            
%             validNucleiIdx = p.nucleiObj.nuclei.status;
            
            [idx, dist] = knnsearch([validNuclei.x, validNuclei.y], [p.spots.x, p.spots.y], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestNucID = validNuclei.nucID(idx);
            p.spots.distanceToNuc = single(dist);
            p.spots.colors = validNuclei.colors(idx, :);
            
        end
        
        function p = assignSpotsInRect(p, rect) %Maybe useful for reassigning spots after add/removing cells 
            
            [spotsInRect, spotIdx] = p.getAllSpotsInRect(rect);
            nucleiNearRect = p.nucleiObj.getNucleiNearRect(rect, p.maxDistance);
            [nucIdx, dist] = knnsearch([nucleiNearRect.x nucleiNearRect.y], [spotsInRect.x spotsInRect.y], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestNucID(spotIdx) = nucleiNearRect.nucID(nucIdx);
            p.spots.distanceToNuc(spotIdx) = single(dist);
            p.spots.colors(spotIdx, :) = nucleiNearRect.colors(nucIdx, :);
            p.allIntensities(); %Because there could be new valid spots to count.
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
                channel = p.spotChannels{i};
                channelIdx = p.maskObj.masks{:, channel};
                maskTable = p.maskObj.masks(channelIdx,:);
                maskIDs = unique(maskTable.maskID);
                maskIDs(maskIDs == 0) = [];
                spotIdx = p.spots.channel == channel;
                p.spots.maskID(spotIdx) =  single(0);
                p.updateSpotStatus(p.spotChannels{i});
                
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
                %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
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
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
            
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
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = removeMasks(p, channel) 
            masksToRemove = setdiff(p.spots.maskID, p.maskObj.masks.maskID(p.maskObj.masks{:,channel}));
            p.spots.maskID(ismember(p.spots.maskID, masksToRemove)) = single(0);
            p.updateSpotStatus(channel);
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = removeMasks2(p, channel, rect)
            %Not sure if it'll be faster to first get spots and masks in
            %rect. 
            [spotsInRect, spotIdx] = p.getSpotsInRect(channel, rect);
            maskIDsInRect = p.maskObj.getChannelMaskIDsInRect(rect, channel);
            maskIDsInRect(maskIDsInRect == 0) = [];
            goodSpotIdx = ~ismember(spotsInRect.maskID, maskIDsInRect);
            spotIdx(spotIdx) = goodSpotIdx;
            p.spots.maskID(spotIdx) = single(0);
            p.updateSpotStatus(channel);
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
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
                & ismember(p.spots.nearestNucID(channelIdx), p.nucleiObj.nuclei.nucID(p.nucleiObj.nuclei.status));%spots assigned to masked cells will be set to false. 
            
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
            %Try using aTrous function for finding spots
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
                    [tempX, tempY, tempIntensity] = d2utils.findSpotsaTrous(tmpPlane);
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
        
        function intensities = getIntensities(p, channel)
            intensities = p.spots{p.spots.channel == channel & p.spots.distanceToNuc <= p.maxDistance, {'intensity'}}; 
        end
        
        function intensities = getIntensitiesNoMasked(p, channel)
            intensities = p.spots{p.spots.channel == channel & p.spots.distanceToNuc <= p.maxDistance & p.spots.maskID == 0, {'intensity'}}; 
        end
        
        function p = allIntensities(p)
            p.spotsIntensitiesWithMasked = cell(0, numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.spotsIntensitiesWithMasked{i} = sort(uint16(p.getIntensities(p.spotChannels{i})));
            end
        end
        
        function p = allIntensitiesNoMasked(p)
            p.spotsIntensitiesNoMasked = cell(0, numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.spotsIntensitiesNoMasked{i} = sort(uint16(p.getIntensitiesNoMasked(p.spotChannels{i})));
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
                p.centroidLists{i}.expression_color = d2utils.expressionToColors(p.centroidLists{i}.GroupCount, p.expressionColorPal{p.paletteIdx});
            end
        end
        
        function p = updateCentroidList(p, channel)
            channelIdx = ismember(p.spotChannels, channel);
            p.centroidLists{channelIdx}...
                = sortrows(p.tabulateChannel(channel), 'GroupCount', 'descend');
            p.centroidLists{channelIdx}.expression_color = d2utils.expressionToColors(p.centroidLists{channelIdx}.GroupCount, p.expressionColorPal{p.paletteIdx});
        end
        
        function p = updateExpressionColors(p)
            for i = 1:numel(p.spotChannels)
                p.centroidLists{i}.expression_color = d2utils.expressionToColors(p.centroidLists{i}.GroupCount, p.expressionColorPal{p.paletteIdx});
            end
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
            tmpSpots = groupsummary(p.spots(p.spots.status,:), {'channel', 'nearestNucID'}, 'IncludeEmptyGroups', true);
            spotLessCells = setdiff(p.nucleiObj.nuclei.nucID(p.nucleiObj.nuclei.status), tmpSpots.nearestNucID);
            spotLessTable = table(repmat(p.spotChannels', numel(spotLessCells), 1), repelem(spotLessCells, 3), zeros(numel(spotLessCells)*3, 1),...
                'VariableNames', tmpSpots.Properties.VariableNames);
            tmpSpots = [tmpSpots; spotLessTable];
            outTable = outerjoin(p.nucleiObj.nuclei(p.nucleiObj.nuclei.status,{'nucID', 'x', 'y'}), tmpSpots, 'Type', 'left', 'LeftKeys', 'nucID', 'RightKeys', 'nearestNucID');
            outTable.nearestNucID = [];
        end
        
        function updateScanSummary(p, varargin)
            if nargin == 1
                inFileName = 'scanSummary.txt';
            elseif nargin == 2
                inFileName = varargin{1};
            end
            
            inFileObj = fopen(inFileName);
            scanSummaryArray = textscan(inFileObj, '%s%q', 'Delimiter', '\t');
            fclose(inFileObj);
            
            if any(ismember(scanSummaryArray{1}, 'spotChannels'))
                idx = ismember(scanSummaryArray{1}, 'spotChannels');
                scanSummaryArray{2}{idx} = strjoin(p.spotChannels);
            else
                newIdx = height(scanSummaryArray{1})+1;
                scanSummaryArray{1}{newIdx} = 'spotChannels';
                scanSummaryArray{2}{newIdx} = strjoin(p.spotChannels);
            end
            
            if any(ismember(scanSummaryArray{1}, 'thresholds'))
                idx = ismember(scanSummaryArray{1}, 'thresholds');
                scanSummaryArray{2}{idx} = num2str([p.thresholds{:}]);
            else
                newIdx = height(scanSummaryArray{1})+1;
                scanSummaryArray{1}{newIdx} = 'thresholds';
                scanSummaryArray{2}{newIdx} = num2str([p.thresholds{:}]);
            end
            
            scanSummaryTable = cell2table(scanSummaryArray{2}, 'RowNames', scanSummaryArray{1}');
            writetable(scanSummaryTable, inFileName, 'WriteRowNames', true, 'WriteVariableNames', false, 'QuoteStrings', false, 'Delimiter', '\t')  
        end
        
        function saveSpotsTable(p, varargin)
            if ~isempty(p.spots)
                if nargin == 1
                    writetable(p.spots(:,{'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'}), 'spots.csv');
                elseif nargin == 2
                    writetable(p.spots(:, {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'}), varargin{1});
                end
            else
                fprintf("spots table is empty. Run findSpots and try again")
            end
        end
        
        function exportSpotsSummary(p, varargin)
            if isempty(p.spots(p.spots.status,:)) || all(p.spots.nearestNucID == 0)
                fprintf("There are no valid spots or spots are not assigned to cells. Try running findSpots and assignSpotsToNuclei")
            else
                outTable = p.tabulateAllChannels;
                if nargin == 1
                    writetable(outTable, 'spotsSummary.csv');
                elseif nargin == 2
                    writetable(outTable, varargin{1});
                end
            end
        end
    end
    
end 