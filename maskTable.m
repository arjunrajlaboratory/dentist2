classdef maskTable < handle
    
    properties (Access = public)
        
        masks
        masksBB
    end
    
    methods
        
        function p = maskTable(scanObject, varargin)
            channels = convertStringsToChars(scanObject.channels);
            channels = channels(~ismember(channels, 'trans'));
            if nargin == 1
                fprintf('New Mask Table\n');
                p.masks = array2table(false(0,numel(channels) + 3), 'VariableNames', [{'maskID','x','y'} , channels]);
                p.masksBB = array2table(false(0,numel(channels) + 3), 'VariableNames', [{'maskID','x','y'} , channels]);
            elseif nargin == 2 % Otherwise, load the specified table
                fprintf('Loading Table\n');
                tmpMasks = readtable(varargin{1},'TextType','string');
                tmpMasks = convertvars(tmpMasks,{'maskID', 'x', 'y'},'single'); %NEED TO FIX: even though converting to single, x and y not == to orignal data
                p.masks = convertvars(tmpMasks,4:width(tmpMasks),'logical');
                p.allMasks2Corners();
            end
        end
        
        function p = addMask(p,maskPoly,localRect,channel)
            
            if isempty(p.masks)
                tempMaskID = single(1);
            else
                tempMaskID = single(max(p.masks.maskID)+1);
            end
            
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            corners = d2utils.polygon2BoundingCorners([x,y]);
            channelIdx = ismember(p.masks.Properties.VariableNames, channel);
            tmpMaskTable = array2table(repmat(channelIdx, length(x), 1), 'VariableNames', p.masks.Properties.VariableNames);
            tmpMaskTable.maskID = repmat(tempMaskID,  length(x), 1);
            tmpMaskTable.x = single(x);
            tmpMaskTable.y = single(y);
            p.masks = [p.masks; tmpMaskTable];
            
            tmpMaskBBTable = array2table(repmat(channelIdx, 4, 1), 'VariableNames', p.masksBB.Properties.VariableNames);
            tmpMaskBBTable.maskID = repmat(tempMaskID,  4, 1);
            tmpMaskBBTable.x =  single(corners(:,1));
            tmpMaskBBTable.y =  single(corners(:,2));
            p.masksBB = [p.masksBB; tmpMaskBBTable];

        end
        
        function p = removeMasks(p,maskIDs)
             p.masks(ismember(p.masks.maskID,maskIDs),:) = [];
             p.masksBB(ismember(p.masksBB.maskID,maskIDs),:) = [];
             
        end
        
        function p = updateMaskPoly(p,channel,maskID,maskPoly,localRect)
            p.removeMasks(maskID);
            p.addMask(maskPoly,localRect,channel)
        end
        
        function outMasks = getAllMasksInRect(p, rect)
            
            idx = p.masksBB.x >= rect(1) & p.masksBB.x < rect(1) + rect(3) ...
                & p.masksBB.y >= rect(2) & p.masksBB.y < rect(2) + rect(4);
            
            maskIDtoKeep = unique(p.masksBB.maskID(idx));
            
            outMasks = p.masks(ismember(p.masks.maskID, maskIDtoKeep) ,:);
        end
        
        function outMasks = getChannelMasksInRect(p, rect, channel)
            channelIdx = ismember(p.masks.Properties.VariableNames, channel);
            
            idx = p.masksBB{:,channelIdx} ...
                & p.masksBB.x >= rect(1) & p.masksBB.x < rect(1) + rect(3) ...
                & p.masksBB.y >= rect(2) & p.masksBB.y < rect(2) + rect(4);
            
            maskIDtoKeep = unique(p.masksBB.maskID(idx));
            
            outMasks = p.masks(ismember(p.masks.maskID, maskIDtoKeep) ,:);
        end
        
        function p = allMasks2Corners(p)
            p.masksBB = array2table(false(0,width(p.masks)), 'VariableNames', p.masks.Properties.VariableNames);
            uniqueMaskIDs = unique(p.masks.maskID);
            for i = 1:numel(uniqueMaskIDs)
                tempMask = p.masks(p.masks.maskID == uniqueMaskIDs(i), :);
                tempMaskBB = repmat(tempMask(1,:), 4, 1);
                tempCorners = d2utils.polygon2BoundingCorners(tempMask{:,{'x', 'y'}});
                tempMaskBB{:,{'x', 'y'}} =  single(tempCorners);
                p.masksBB = [p.masksBB; tempMaskBB];
            end
        end
       
        function [] = saveTables(p, varargin)
           if nargin == 1
               writetable(p.masks, 'masks.csv');
           elseif nargin == 2
               writetable(p.masks, varargin{1});
               
           end
        end
        
    end
    
end 