classdef cellTable < handle
    
    properties (Access = public)
        
        cells
        masks %Do we want to mask cells? 
        
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
            
            outCells = p.cells(p.cells.status==1,:);

            idx = outCells.x >= rect(1) & outCells.x < rect(1) + rect(3) ...
                & outCells.y >= rect(2) & outCells.y < rect(2) + rect(4);
            
            outCells = outCells(idx,:);
        end
                
        function set.minNucleusSize(p,area)
            p.minNucleusSize = area;
            %UPDATE CELLS
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
                [idx, dist] = knnsearch(tempCells{:,{'x','y'}}, [x, y], 'K', 1, 'Distance', 'euclidean');
                if dist < 40 %make sure that the query point is near a cell
                    tempCells(idx,:) = [];
                    p.cells = tempCells;
                end
                % UPDATE SPOTS
            else
                disp("cells is empty.")
            end
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