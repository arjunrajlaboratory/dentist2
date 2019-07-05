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
                p.cells = cell2table(cell(0,5), 'VariableNames', {'cellID', 'x', 'y', 'status', 'fileID'});
            else
                fprintf('Loading Table\n');
                p.cells = readtable(varargin{1},'TextType','string');
            end
        end
        
        
        function p = findCells(p,fileTable)
            files = fileTable.files(fileTable.files.channel=="dapi",:);
            
            cellID = [];
            x = [];
            y = [];
            status = [];
            fileID = [];
            
            for i = 1:size(files,1)
                fprintf('Finding nuclei in %s\n',files.fileName(i));
                %dapiIm = imread(files.fileName(i));
                [dapiIm, outRect] = fileTable.getPaddedImage(files.fileID(i),100);
                T = adaptthresh(dapiIm,0.3,'ForegroundPolarity','bright');
                dp = imbinarize(dapiIm,T);
                
                CC = bwconncomp(dp);
                rp = regionprops(CC);
                area = [rp.Area];
                
                idx = area > p.minNucleusSize; % Get rid of small stuff
                
                rp = rp(idx);
                centroids = [rp.Centroid];
                centroids = round(reshape(centroids,2,[])');
                % Second column is the rows, first column is the column
                centroids = fliplr(centroids);
                % first column is the rows, second column is the column
                
                % First, let's adjust to global coordinates
                centroidRow = centroids(:,1) + outRect(1) - 1;
                centroidCol = centroids(:,2) + outRect(2) - 1;
                
                % Now crop
                top = files.top(i);
                left = files.left(i);
                height = files.height(i);
                width = files.width(i);
                
                idx = centroidRow >= top & centroidRow < top + height & centroidCol >= left & centroidCol < left + width;
                
                centroidRow = centroidRow(idx);
                centroidCol = centroidCol(idx);
                
                x = [x ; centroidRow];
                y = [y ; centroidCol];
                status = [status ; ones(numel(centroidRow),1)];
                fileID = [fileID ; repmat(files.fileID(i),numel(centroidRow),1)];
            end
            
            p.cells = table((1:length(x))',x,y,status,fileID,'VariableNames', {'cellID', 'x', 'y', 'status', 'fileID'});
        end
        
        function outCells = getCellsInRect(p,rect)
            %outCells = p.cells(p..channel == channel,:);
            outCells = p.cells(p.cells.status==1,:);

            idx = outCells.x >= rect(1) & outCells.x < rect(1) + rect(3) ...
                & outCells.y >= rect(2) & outCells.y < rect(2) + rect(4);
            
            outCells = outCells(idx,:);
        end

        
    end
    
end 