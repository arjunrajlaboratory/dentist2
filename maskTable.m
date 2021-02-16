classdef maskTable < handle
    
    properties (Access = public)
        
        masks
        masksBB
        channels
    end
    
    methods
        
        function p = maskTable(scanObject, varargin)
            channels = convertStringsToChars(scanObject.channels);
            p.channels = channels(~ismember(channels, 'trans'));
            if nargin == 1
                fprintf('New Mask Table\n');
                p.masks = table('Size', [0, numel(p.channels) + 3],...
                    'VariableNames', [{'maskID','x','y'} , p.channels],...
                    'VariableTypes', [repmat({'single'}, 1, 3) , repmat({'logical'}, 1, numel(p.channels))]);
                p.masksBB = table('Size', [0, numel(p.channels) + 2],...
                    'VariableNames', [{'maskID','BB'} , p.channels],...
                    'VariableTypes', [repmat({'single'}, 1, 2) , repmat({'logical'}, 1, numel(p.channels))]);
%                 p.masks = table(false(0,numel(p.channels) + 3), 'VariableNames', [{'maskID','x','y'} , p.channels]);
%                 p.masksBB = array2table(false(0,numel(p.channels) + 5), 'VariableNames', [{'maskID','x','y', 'h', 'w'} , p.channels]);
            elseif nargin == 2 % Otherwise, load the specified table
                fprintf('Loading Table\n');
                tmpMasks = readtable(varargin{1},'TextType','string');
                tmpMasks = convertvars(tmpMasks,{'maskID', 'x', 'y'},'single');
                p.masks = convertvars(tmpMasks,4:width(tmpMasks),'logical');
                p.allMasks2BB();
            end
        end
        
        function p = addMask(p,maskPoly,localRect,channel)
            
            if isempty(p.masks)
                tempMaskID = single(1);
            else
                tempMaskID = single(max(p.masks.maskID)+1);
            end
            
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            channelIdx = ismember(p.channels, channel);
            
            tmpCoords = table(repmat(tempMaskID,  length(x), 1), single(x), single(y), 'VariableNames',{'maskID', 'x', 'y'});
            tmpChannelTable = array2table(repmat(channelIdx, length(x), 1), 'VariableNames', p.channels);

            p.masks = [p.masks; [tmpCoords, tmpChannelTable]];

            BB = d2utils.polygonBoundingBox([x,y]);
            tmpCoords = cell2table({tempMaskID, single(BB)}, 'VariableNames',{'maskID', 'BB'});
            p.masksBB = [p.masksBB; [tmpCoords, tmpChannelTable(1,:)]];

        end
        
        function p = addMaskLocalCoords(p,maskPoly,channel)
            
            if isempty(p.masks)
                tempMaskID = single(1);
            else
                tempMaskID = single(max(p.masks.maskID)+1);
            end
            
%             [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            channelIdx = ismember(p.channels, channel);
            
            tmpCoords = table(repmat(tempMaskID,  height(maskPoly), 1), single(maskPoly(:,2)), single(maskPoly(:,1)), 'VariableNames',{'maskID', 'x', 'y'});
            tmpChannelTable = array2table(repmat(channelIdx, height(maskPoly), 1), 'VariableNames', p.channels);

            p.masks = [p.masks; [tmpCoords, tmpChannelTable]];

            BB = d2utils.polygonBoundingBox(fliplr(maskPoly));
            tmpCoords = cell2table({tempMaskID, single(BB)}, 'VariableNames',{'maskID', 'BB'});
            p.masksBB = [p.masksBB; [tmpCoords, tmpChannelTable(1,:)]];

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
            
            idx = rectint(p.masksBB.BB,rect) > 0;
                        
            outMasks = p.masks(ismember(p.masks.maskID, p.masksBB.maskID(idx)) ,:);
        end
        
        function outMasks = getChannelMasksInRect(p, rect, channel)
            
            idx = rectint(p.masksBB.BB, rect) > 0 & p.masksBB{:, channel};
           
            outMasks = p.masks(ismember(p.masks.maskID, p.masksBB.maskID(idx)) ,:);
        end
        
        function p = removeMasksByPoints(p, points, rect, channel)
            masksInRect = p.getChannelMasksInRect(rect, channel);
            maskIDs = unique(masksInRect.maskID);
            [x,y] = d2utils.localToGlobalCoords(rect,points(:,2),points(:,1));
            for i = 1:numel(maskIDs)
                if any(inpolygon(x, y, masksInRect{masksInRect.maskID == maskIDs(i), 'x'}, masksInRect{masksInRect.maskID == maskIDs(i), 'y'}))
                    p.removeMasks(maskIDs(i));
                end
            end
        end
        
        function p = removeMasksByLocalPoints(p, points, rect)
            masksInRect = p.getAllMasksInRect(rect);
            maskIDs = unique(masksInRect.maskID);
            for i = 1:numel(maskIDs)
                if any(inpolygon(points(:,2), points(:,1), masksInRect{masksInRect.maskID == maskIDs(i), 'x'}, masksInRect{masksInRect.maskID == maskIDs(i), 'y'}))
                    p.removeMasks(maskIDs(i));
                end
            end
        end
        
        function p = allMasks2BB(p)
            p.masksBB = table('Size', [0, numel(p.channels) + 2],...
                    'VariableNames', [{'maskID','BB'} , p.channels],...
                    'VariableTypes', [repmat({'single'}, 1, 2) , repmat({'logical'}, 1, numel(p.channels))]);
            uniqueMaskIDs = unique(p.masks.maskID);
            for i = 1:numel(uniqueMaskIDs)
                tempMask = p.masks(p.masks.maskID == uniqueMaskIDs(i), :);
                BB = d2utils.polygonBoundingBox(tempMask{:,{'x', 'y'}});
                tmpCoords = cell2table({uniqueMaskIDs(i), single(BB)}, 'VariableNames',{'maskID', 'BB'});
                p.masksBB = [p.masksBB; [tmpCoords, tempMask(1,4:end)]];
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