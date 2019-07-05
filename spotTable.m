classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        masks % per channel
        thresholds % per channel
        radius = 300;
        spotChannels
        
        theFilter
        percentileToKeep = 98;
        
    end
    
    methods
        
        function p = spotTable(fileTable)
            p.spots = [];
            p.masks = [];
            p.thresholds = [];
            channels = fileTable.channels;
            p.spotChannels = channels(~ismember(channels,["dapi","trans"]));
            
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
                p.masks = cell2table(cell(0,3), 'VariableNames', {'maskID', 'channel', 'polygon'});
                tempMaskID = 1;
            else
                tempMaskID = max(p.masks.maskID)+1;
            end
            
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            maskPoly = [x y];
            
            p.masks = [p.masks; table(tempMaskID,channel,{maskPoly},'VariableNames', {'maskID', 'channel', 'polygon'})];
            
            % INSERT CODE HERE TO UPDATE SPOTS
            
        end
        
        function p = removeMasks(p,maskIDs)
            idx = ismember(p.masks.maskID,maskIDs);
            % NEED TO save which ones got deleted for updating.
            p.masks(idx,:) = [];
            
            % INSERT CODE HERE TO UPDATE
        end
        
        function p = updateMaskPoly(p,maskID,maskPoly,localRect)
            idx = find(p.masks.maskID == maskID);
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            maskPoly = [x y];
            p.masks.polygon(idx) = {maskPoly};
            % INSERT CODE HERE TO UPDATE
        end
        
        function p = assignSpotsToCells(p,cellTable)
            
            rect = zeros(1,4);
            rect(3) = 2*p.radius+1; % Size of rectangle will always be the same
            rect(4) = 2*p.radius+1;
            for i = 1:height(p.spots)
                if mod(i,100) == 1
                    fprintf('%d of %d spots\n',i,height(p.spots));
                end
                rect(1) = p.spots.x(i)-p.radius;
                rect(2) = p.spots.y(i)-p.radius;
                cells = cellTable.getCellsInRect(rect);
                if isempty(cells)
                    p.spots.nearestCellID(i) = 0; % Should we set this to nan?
                else
%                     if height(cells) == 1 % If there's just one nearby cell, just handle it
%                         p.spots.nearestCellID(i) = cells.cellID;
%                     else
                        distMatrix = [(cells.x-p.spots.x(i))';(cells.y-p.spots.y(i))'];
                        nm = vecnorm(distMatrix);
                        [~,idx] = sort(nm);
                        if nm(idx(1))<p.radius
                            p.spots.nearestCellID(i) = cells.cellID(idx(1));
                        end
%                     end
                end
            end
                
            
%             for i = 1:height(cellTable.cells)
%                 x = cellTable.cells.x(i);
%                 y = cellTable.cells.y(i);
%                 
%                 rect = [x-p.radius,y-p.radius,2*p.radius+1,2*p.radius+1];
%                 [outSpots,idx] = p.getAllSpotsInRect(rect);
%                 dists = norm
%                 p.spots.nearestCellID(idx) = cellTable.cells.cellID(i);
%                 
%             end
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
        
        function findSpots(p,fileTable)
            
            %p.spots = cell2table(cell(0,8), 'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'status', 'channel', 'fileID'});
            
            %spotID = [];
            x = [];
            y = [];
            intensity = [];
            %nearestCellID = [];
            %status = [];
            channel = [];
            fileID = [];
            
            p.theFilter = -fspecial('log',20,2);
            for i = 1:length(p.spotChannels)
                currChannel = p.spotChannels(i);
                fprintf('Channel: %s\n',currChannel);
                files = fileTable.files(fileTable.files.channel == currChannel,:);
                for j = 1:size(files,1)
                    fprintf('File: %s\n',files.fileName(j));
                    im = imread(files.fileName(j));
                    % Convert to single to "smooth" things so that
                    % imregional max will isolate a particular pixel
                    filt = imfilter(im2single(im),p.theFilter,'replicate');
                    irm = imregionalmax(filt);
                    tempSpots = filt(irm)';
                    thresh = prctile(tempSpots,p.percentileToKeep);
                    %idx = tempSpots > thresh;
                    filt = filt.*single(irm);
                    filt(filt < thresh) = 0;
                    filt = im2uint16(filt);
                    idx = filt > 0;
                    intensities = filt(idx);
                    [tempX,tempY] = ind2sub(size(filt),find(idx));
                    
                    tempX = tempX + files.top(j) - 1;
                    tempY = tempY + files.left(j) - 1; 
                    
                    x = [x ; tempX];
                    y = [y ; tempY];
                    intensity = [intensity ; intensities];
                    channel = [channel ; repmat(currChannel,length(tempX),1)];
                    fileID = [fileID ; repmat(files.fileName(j),length(tempX),1)];
                end
            end
            
            spotID = (1:length(x))';
            nearestCellID = nan(length(x),1);
            maskID = zeros(length(x),1);
            status = true(length(x),1);
            
            p.spots = table(spotID, x, y, intensity, nearestCellID, status, maskID, channel, fileID);
            
%             p.spots.spotID = spotID;
%             p.spots.x = x;
%             p.spots.y = y;
%             p.spots.intensity = intensity;
%             p.spots.status = status;
%             p.spots.nearestCellID = nearestCellID;
%             p.spots.fileID = fileID;
%             p.spots.channel = channel;
%             
%             
        end

        
    end
    
end 