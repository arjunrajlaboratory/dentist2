classdef cellTable < handle
    
    properties (Access = public)
        
        cells
        masks
        
        minNucleusSize = 1000; % Set method should probably call findCells

    end
    
    methods
        
        function p = cellTable(varargin)
            if nargin == 0
                fprintf('New Table\n');
                p.cells = cell2table(cell(0,5), 'VariableNames', {'cellID', 'x', 'y', 'status', 'area'}); 
            else
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
            
            area = [rp.Area];
            status = ones(size(centroids, 1),1);
            
            
            p.cells = table((1:size(centroids, 1))',centroids(:,2),centroids(:,1),status, area', 'VariableNames', {'cellID', 'x', 'y', 'status', 'area'});
        end
        
        function outCells = getCellsInRect(p,rect) %rect specified as [x y width height]
            
            outCells = p.cells(p.cells.status==1,:);

            idx = outCells.x >= rect(1) & outCells.x < rect(1) + rect(3) ...
                & outCells.y >= rect(2) & outCells.y < rect(2) + rect(4);
            
            outCells = outCells(idx,:);
        end
                
        function set.minNucleusSize(p,area)
            p.minNucleusSize = area;
        end
        
        function p = addCell(p, x, y, status, area)
            if ~isempty(p.cells)
                tempCells = p.cells;
                maxID = max(p.cellID);
                newCell = {maxID+1, x, y, status, area};
                tempCells = [tempCells;newCell];
                p.cells = tempCells;
            else
                p.cells(1,:) = {1, x, y, status, area};
            end
        end
        
        function p = removeCell(p, x, y)
            if ~isempty(p.cells)
                tempCells = p.cells;
                [idx, dist] = dsearchn(tempCells{:,{'x','y'}}, [x, y]);
                if dist < 40 %make sure that the query point is near a cell
                    tempCells(idx,:) = [];
                    p.cells = tempCells;
                end
            else
                disp("cells is empty.")
            end
        end

        
    end
    
end 