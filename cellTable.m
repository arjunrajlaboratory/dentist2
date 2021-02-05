classdef cellTable < handle
    
    properties (Access = public)
        
        cells
        cellMasks %Do we want to mask cells? 
        cellMasksBB
        minNucleusSize = 1000; % Set method should probably call findCells

    end
    
    methods
        
        function p = cellTable(varargin)
            if nargin == 0
                fprintf('New Table\n');
                p.cells = cell2table(cell(0,6), 'VariableNames', {'cellID', 'x', 'y', 'status', 'maskID', 'nucleusArea'}); 
            elseif nargin == 1
                fprintf('Loading Table\n');
                p.cells = readtable(varargin{1},'TextType','string');
            end
        end
        
        
        function p = findCells(p,dapiMask)
            
            CC = bwconncomp(dapiMask);
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

            p.cells = table((1:height(centroids))',centroids(:,2),centroids(:,1), status, maskID, area',...
                'VariableNames', {'cellID', 'x', 'y', 'status', 'maskID', 'nucleusArea'});
        end
        
        function outCells = getCellsInRect(p,rect) %rect specified as [x y nrows ncols]
            
            %outCells = p.cells(p.cells.status==1,:);

            idx = p.cells.x >= rect(1) & p.cells.x < rect(1) + rect(3) ...
                & p.cells.y >= rect(2) & p.cells.y < rect(2) + rect(4) ...
                & p.cells.status;
            
            outCells = p.cells(idx,:);
        end
        
        function outCells = getCellsNearRect(p,rect, radius) %rect specified as [x y nrows ncols]
            
            %outCells = p.cells(p.cells.status==1,:);

            idx = p.cells.x >= rect(1) - radius & p.cells.x < rect(1) + rect(3) + radius ...
                & p.cells.y >= rect(2) - radius & p.cells.y < rect(2) + rect(4) + radius ...
                & p.cells.status;
            
            outCells = p.cells(idx,:);
        end
                
        function p = setMinNucleusSize(p,area, dapiMask)
            p.minNucleusSize = area;
            %UPDATE CELLS. Shoudl we just change status instead? 
            p.findCells(dapiMask)
        end
        
        function p = addCell(p, x, y, status, area)
            if ~isempty(p.cells)
                tempCells = p.cells;
                maxID = max(p.cellID);
                newCell = {maxID+1, x, y, status, 0, area};
                tempCells = [tempCells;newCell];
                p.cells = tempCells;
                % UPDATE SPOTS?
            else
                p.cells(1,:) = {1, x, y, status, 0, area};
            end
        end
        
        function p = removeCell(p, x, y)
            if ~isempty(p.cells)
                tempCells = p.cells;
                [idx, dist] = knnsearch(p.cells{:,{'x','y'}}, [x, y], 'K', 1, 'Distance', 'euclidean');
                if dist < 40 %make sure that the query point is near a cell
                    p.cells(idx,:) = [];
                end
                % UPDATE SPOTS
            else
                disp("cells is empty.")
            end
        end
        
        function outCells = getValidCellsInRect(p, rect)
            
            idx = p.cells.status ... 
                & p.cells.x >= rect(1) & p.cells.x < rect(1) + rect(3) ...
                & p.cells.y >= rect(2) & p.cells.y < rect(2) + rect(4);
            
            outCells = p.cells(idx,:);
            
        end
        
        function p = addMask(p,maskPoly,localRect)
            if isempty(p.cellMasks)
                tempMaskID = single(1);
            else
                tempMaskID = single(max(p.cellMasks.maskID)+1);
            end
            
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            %maskPoly = [x y];
            corners = d2utils.polygon2BoundingCorners([x,y]);
            
            p.cellMasks = [p.cellMasks; table(repmat(tempMaskID, length(x), 1),single(x),single(y), 'VariableNames', {'maskID', 'x', 'y'})];
            p.cellMasksBB = [p.cellMasksBB; table(repmat(tempMaskID, 4, 1),corners(:,1),corners(:,2), 'VariableNames', {'maskID', 'x', 'y'})];
            
            %UPDATE CELLS - inpolygon seems to fast enough for many millions of cells. If it gets too slow, could try selecting cells in localRect first.
%             tempcells = p.getValidCellsInRect(localRect);
%             tempIdx = inpolygon(tempcells.x,tempcells.y,x,y);
            idx = inpolygon(p.cells.x, p.cells.y, x, y) & p.cells.status;
            p.cells.maskID(idx) = tempMaskID;
            p.cells.status(idx) = false;
        end
        
        function p = removeMasks(p,maskIDs)
             
            p.masks(ismember(p.masks.maskID,maskIDs),:) = [];
            p.masksBB(ismember(p.masksBB.maskID,maskIDs),:) = [];
            
            p.cells(ismember(p.cells.maskID,maskIDs), 'maskID') = single(0);
            p.cells(ismember(p.cells.maskID,maskIDs), 'status') = true; %This assumes that status is determined only by masks
            
        end
        
        function [] = saveCellTable(p, varargin)
            if ~isempty(p.cells)
                if nargin == 1
                    writetable(p.cells, 'cells.csv')
                elseif nargin ==2 
                    writetable(p.cells, varargin{1})
                end
            else
                fprintf("cells is empty. Run findCells and try again")
            end
        end
        
    end
    
end 