classdef nucleiTable < handle
    
    properties (Access = public)
        
        nuclei
        minNucleusSize = 1000; %Update set method to call findNuclei and updateAllMasks
        dapiMask
        
        scanObj
        maskObj  
       
    end
    
    methods
        
        function p = nucleiTable(scanObject, maskObj, varargin)
            p.scanObj = scanObject;
            p.maskObj = maskObj;
            if nargin == 2
                fprintf('New Table\n');
                p.nuclei = cell2table(cell(0,6), 'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea'}); 
            elseif nargin == 3
                fprintf('Loading Table\n');
                p.nuclei = readtable(varargin{1},'TextType','string');
            end
        end
        
        function p = stitchDAPImask(p, varargin)  
            tileTable = p.scanObj.tilesTable;
            tilesTmp = transpose(p.scanObj.scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = p.scanObj.tileSize();
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
            
            p.dapiMask = tmpStitch;
        end
        
        function p = findNuclei(p)
            
            CC = bwconncomp(p.dapiMask);
            rp = regionprops(CC);
            area = [rp.Area];
            idx = area >= p.minNucleusSize;
            rp = rp(idx);
            centroids = [rp.Centroid];
            centroids = round(reshape(centroids,2,[])');
            centroids = single(centroids);
            
            area = single([rp.Area]);
            
            status = true(height(centroids), 1);
            maskID = single(zeros(height(centroids), 1));

            p.nuclei = table((1:height(centroids))',centroids(:,2),centroids(:,1), status, maskID, area',...
                'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea'});
        end
        
        function [outNuclei, idx] = getNucleiInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.nuclei.x >= rect(1) & p.nuclei.x < rect(1) + rect(3) ...
                & p.nuclei.y >= rect(2) & p.nuclei.y < rect(2) + rect(4);
            
            outNuclei = p.nuclei(idx,:);
        end
        
        
        function outNuclei = getValidNucleiInRect(p, rect)
            
            idx = p.nuclei.status ... 
                & p.nuclei.x >= rect(1) & p.nuclei.x < rect(1) + rect(3) ...
                & p.nuclei.y >= rect(2) & p.nuclei.y < rect(2) + rect(4);
            
            outNuclei = p.nuclei(idx,:);
            
        end
        
        function idx = getNucleiInRectIdx(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.nuclei.x >= rect(1) & p.nuclei.x < rect(1) + rect(3) ...
                & p.nuclei.y >= rect(2) & p.nuclei.y < rect(2) + rect(4) ...
                & p.nuclei.status;
        end
        
        function outNuclei = getNucleiNearRect(p,rect, radius) %rect specified as [x y nrows ncols]
            
            idx = p.nuclei.x >= rect(1) - radius & p.nuclei.x < rect(1) + rect(3) + radius ...
                & p.nuclei.y >= rect(2) - radius & p.nuclei.y < rect(2) + rect(4) + radius ...
                & p.nuclei.status;
            
            outNuclei = p.nuclei(idx,:);
        end
        
        function p = addCell(p, x, y)  
            if ~isempty(p.nuclei)
                maxID = max(p.nucID);
                newNuc = {maxID+1, single(x), single(y), true, single(0), single(0)}; %Should area be NaN?
                p.nuclei = [p.nuclei; newNuc];
                
            else
                p.nuclei(1,:) = {single(1), single(x), single(y), true, single(0), single(area)};
            end
        end
        
        function p = removeCell(p, x, y)
            if ~isempty(p.nuclei)
                [idx, dist] = knnsearch(p.nuclei{:,{'x','y'}}, [x, y], 'K', 1, 'Distance', 'euclidean');
                if dist < 40 %make sure that the query point is near a nucleus
                    p.nuclei(idx,:) = [];
                end
                
            else
                disp("NucleiTable is empty.")
            end
        end
        
        function p = updateAllMasks(p) %This will overwrite previous maskIDs in nuclei. 
            
            maskTable = p.maskObj.masks(p.maskObj.masks.dapi, :);
            maskIDs = unique(maskTable.maskID);
            p.nuclei.maskID(:) = single(0);
            p.nuclei.status(:) = true;
            for i = 1:numel(maskIDs)
                idx = inpolygon(p.nuclei.x, p.nuclei.y,...
                    maskTable{maskTable.maskID == maskIDs(i), 'x'}, maskTable{maskTable.maskID == maskIDs(i), 'y'}) & p.nuclei.status;
                p.nuclei.maskID(idx) = maskIDs(i);
                p.nuclei.status(idx) = false;
            end
            
        end
        
        function p = addNewMask(p)
            
            maxCellMask = max(p.maskObj.masks{p.maskObj.masks.dapi,'maskID'});
            maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxCellMask,{'x', 'y'}}; %Only query nuclei within mask bouding box
            polyRect = d2utils.boundingCorners2Rect(maskBB);
            
            nucIdx = p.getNucleiInRectIdx(polyRect);
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
            %If nucleus falls within multiple masks and only 1 is removed,
            %this function may incorrectly set status to true. Can instead
            %use updateMasksInRect
            
            masksToRemove = setdiff(p.nuclei.maskID, p.maskObj.masks.maskID(p.maskObj.masks.dapi));
            p.nuclei.maskID(ismember(p.nuclei.maskID, masksToRemove)) = single(0);
            p.nuclei.status(ismember(p.nuclei.maskID, masksToRemove)) = true;
        end
        
        
        function [] = saveNucleiTable(p, varargin)
            if ~isempty(p.nuclei)
                if nargin == 1
                    writetable(p.nuclei, 'nuclei.csv')
                elseif nargin ==2 
                    writetable(p.nuclei, varargin{1})
                end
            else
                fprintf("nuclei is empty. Run findNuclei and try again")
            end
        end
        
    end
    
end 