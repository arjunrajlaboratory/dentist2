classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        
        thresholds
        spotChannels
        maxDistance = 100;  
        theFilter
        percentileToKeep = 98;
        
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
                fprintf('Loading spot table');
                tmpSpots = readtable(n.Results.spotsFile,'TextType','string');
                p.spots = convertvars(tmpSpots, {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'maskID', 'distanceToNuc'}, 'single');
                p.spots.status = logical(p.spots.status);
            end
                        
            if isempty(n.Results.thresholdsFile)
                fprintf('New thresholds table\n');
                p.thresholds = cell2table(cell(0,numel(p.spotChannels)), 'VariableNames', p.spotChannels);
            else
                fprintf('Loading thresholds table');
                p.thresholds = readtable(n.Results.thresholdsFile,'TextType','string');
            end
        end
        
        
        function p = assignSpotsToNuclei(p)
            validNuclei = p.nucleiObj.nuclei(p.nucleiObj.nuclei.status);
            
            [idx, dist] = knnsearch([validNuclei.x, validNuclei.y], [p.spots.x, p.spots.y], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestNucID = validNuclei.nucID(idx);
            p.spots.distanceToNuc = single(dist);
            
        end
        
        function p = assignSpotsInRect(p, channel, rect) %Maybe useful for reassigning spots after add/removing cells 
            
            spotIdx = p.getValidSpotsInRectIndex(channel,rect);
            nucleiNearRect = p.nucleiObj.getNucleiNearRect(rect, p.maxDistance);
            [nucIdx, dist] = knnsearch([nucleiNearRect.x nucleiNearRect.y], [p.spots.x(spotIdx) p.spots.y(spotIdx)], 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestNucID(spotIdx) = nucleiNearRect.nucID(nucIdx);
            p.spots.distanceToNuc(spotIdx) = single(dist);
           
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
        
        function p = updateAllMasks(p)  %Note, this will overwrite previous maskIDs in spotTable. 
            
            for i = 1:numel(p.spotChannels)
                channelIdx = p.maskObj.masks(:, p.spotChannels(i));
                maskTable = p.maskObj.masks(channelIdx,:);
                maskIDs = unique(maskTable.maskID);
                spotIdx = p.spots.channel == channel;
                spotTable = p.spots(spotIdx,:);
                spotTable.maskID = single(0);
                for ii = 1:numel(maskIDs)
                    idx = inpolygon(spotTable.x, spotTable.y,...
                        maskTable{maskTable.maskID == maskIDs(ii), 'x'}, maskTable{maskTable.maskID == maskIDs(ii), 'y'}) & spotTable.maskID == 0;
                    spotTable.maskID(idx) = maskIDs(ii); 
                    spotTable.status(idx) = false; 
                end
                p.spots.maskID(spotIdx) = spotTable.maskID;
                p.spots.status(spotIdx) = spotTable.status;
            end 
            
        end
        
        function p = addNewMaskInRect(p, channel)
            
            maxSpotMask = max(p.maskObj.masks{p.maskObj.masks{:,channel},'maskID'});
            maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxSpotMask,{'x', 'y'}}; %Only query spots within mask bouding box
            polyRect = d2utils.boundingCorners2Rect(maskBB);
            spotIdx = getValidSpotsInRectIndex(channel, polyRect);
            idx = inpolygon(p.spots.x(spotIdx), p.spots.y(spotIdx),...
                p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'y'});
            
            spotIdx(spotIdx) = idx; %Only spots in polygon remain true
                
            p.spots.maskID(spotIdx) = maxSpotMask;
            p.spots.status(spotIdx) = false;
            
        end
        
%         function p = updateMasksInRect(p, channel, localRect) %%NEEDS TO
%         BE FIXED
%             %Use for removing masks or if we want to change multiple masks
%             %before updating spots
%             %Probably could use some speeding up. 
%             masksInRect = p.maskObj.getChannelMasksInRect(localRect, channel);
%             maskIDsinRect = unique(masksInRect.maskID);
%             [tmpSpots, spotIdx] = getValidSpotsInRect(channel, localRect);            
%             
%             %Resest status for nuclei
%             p.spots.maskID(spotIdx) = single(0);
%             tmpSpots.status = true;
%             for i = 1:numel(maskIDsinRect)
%                 idxPoly = inpolygon(tmpSpots.x, tmpSpots.y,...
%                     masksInRect{masksInRect.maskID == maskIDsinRect(i), 'x'}, masksInRect{masksInRect.maskID == maskIDsinRect(i), 'y'}) & tmpSpots.status;
%                 
%                 tmpSpots.maskID(idxPoly) = maskIDsinRect(i);
%                 tmpSpots.status(idxPoly) = false;
%                 
%             end
%             
%             p.spots.maskID(spotIdx) = tmpSpots.maskID;
%             p.updateSpotStatus();
%             
%             p.spots.status(spotIdx) = tmpSpots.status;
%  
%         end
        
        function p = removeMasks(p, channel) 
            %If spot falls within multiple masks and only 1 is removed,
            %this function may incorrectly set status to true. Can instead
            %use updateMasksInRect
            
            masksToRemove = setdiff(p.nuclei.maskID, p.maskObj.masks.maskID(p.maskObj.masks{:,channel}));
            p.spots.maskID(ismember(p.spots.maskID, masksToRemove)) = single(0);
            p.spots.status(ismember(p.spots.maskID, masksToRemove)) = true;
        end
       
        function p = updateSpotStatus(p,channel)
            threshold = p.thresholds{p.spotChannels == channel};
            spotIdx = p.spots.channel == channel & p.spots.intensity >= threshold ...
                & p.spots.distanceToNuc <= p.maxDistance & p.spots.maskID == 0;
            
            p.spots.status(spotIdx) = true;
            p.spots.status(~spotIdx) = false;
                    
        end
        
        function p = updateAllSpotStatus(p)
            for i = 1:numel(p.spotChannels)
                p.updateSpotStatus(p.spotChannels{i})
            end
        end
        
        function intensities = getIntensities(p,channel)
            intensities = p.spots{p.spots.channel == channel,'intensity'};
        end
        
        function intensityThresh = getIntensityThreshold(p,channel)
            intensityThresh = p.thresholds{p.spots.channel == channel};
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
            nearestNucID = single(zeros(length(x),1));
            maskID = single(zeros(length(x),1));
            status = true(length(x),1);
            dist =  single(zeros(length(x),1));
            p.spots = table(spotID, single(x), single(y), intensity, nearestNucID, status, maskID, channel, dist,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'});
        end
        
        function outTable = getSpotsPerNuc(p)
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
            n.addParameter('thresholdsFile', 'thresholds.csv', @(x) assert((ischar(x) & endsWith(x, '.csv')),...
                'Specify .csv filename to save thresholds. e.g "thresholds.csv"'));

            n.parse(varargin{:});
            
            writetable(p.spots, n.Results.spotsFile);
           % writetable(p.thresholds, n.Results.thresholdsFile);
        end
        
        
    end
    
end 