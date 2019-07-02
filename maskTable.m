classdef maskTable < handle
    
    properties (Access = public)
        
        masks
        
    end
    
    methods
        
        function p = maskTable(varargin)
            if nargin == 0
                fprintf('New Table\n');
                p.allPoints = cell2table(cell(0,6), 'VariableNames', {'spotID','x','y','intensity','channel','nearestCellID','isGood','fileID'});
            else % Otherwise, load the specified table
                fprintf('Loading Table\n');
                p.allPoints = readtable(varargin{1},'TextType','string');
            end
        end
        
        function p = changeAnnotation(p, pointIDs, newAnnotation)
            idx = ismember(p.allPoints.pointID,pointIDs);
            p.allPoints.annotation(idx) = newAnnotation;
        end
        
        function [p, newPoints] = addRawPoints(p, frame, centroids)
            sz = size(centroids);
            
            mxID = max(p.allPoints.pointID); % Will assign pointIDs to the end of existing list.
            if isempty(mxID)
                mxID = 0;
            end
            
            T = table( ((mxID+1):(mxID+sz(1)))', frame*ones(sz(1),1), centroids(:,1), centroids(:,2), nan(sz(1),1));
            T.Properties.VariableNames = {'pointID','frameNumber','xCoord','yCoord','parentID'};
            T.annotation = repmat("none",sz(1),1);
            p.allPoints = [p.allPoints; T];
            newPoints = T; % This allows you to capture just the new points that you just added, but in the table format
        end
        
        function p = removePoints(p,removeIDs)
            idx = ismember(p.allPoints.pointID,removeIDs);
            p.allPoints(idx,:) = [];
            idx = ismember(p.allPoints.parentID,removeIDs);
            p.allPoints.parentID(idx) = nan;
        end
        
        function p = guessParents(p)
            
            endFrame = max(unique(p.allPoints.frameNumber));
            
            for i = 2:endFrame
                prevFrame = i-1;
                prevList = p.allPoints(p.allPoints.frameNumber == prevFrame,:);
                currList = p.allPoints(p.allPoints.frameNumber == i, :);
                
                xy1 = [prevList.xCoord prevList.yCoord];
                xy2 = [currList.xCoord currList.yCoord];
                
                D = pdist2(xy1,xy2);
                [~,idx] = min(D);
                
                parentIDs = prevList.pointID(idx);
                
                p.allPoints(p.allPoints.frameNumber == i,:).parentID = parentIDs;
            end
        end
        
        function p = assignParent(p, childID, parentID)
            id = p.allPoints.pointID == childID;
            p.allPoints.parentID(id) = parentID;
        end
        
        function p = setPointCoordinates(p, pointID, xCoord, yCoord)
            p.allPoints(p.allPoints.pointID == pointID,:).xCoord = xCoord;
            p.allPoints(p.allPoints.pointID == pointID,:).yCoord = yCoord;
        end
        
        function p = setPointParent(p, pointID, parentID)
            p.allPoints(p.allPoints.pointID == pointID,:).parentID = parentID;
        end
        
        function outTab = getAllPointsInFrame(p,frameNumber)
            outTab = p.allPoints(p.allPoints.frameNumber == frameNumber,:);
        end
        
        % Probably need something to set an entire frame's worth of points
        % UNTESTED!! This also might be done better in a more abstract way
        % Also, not clear this is actually needed.
        function p = updateEntireFrame(p,frameNumber,newTable)
            idx = p.allPoints.frameNumber == frameNumber;
            p.allPoints(idx,:) = [];
            p.allPoints = [p.allPoints newTable];
        end

        
    end
    
end 