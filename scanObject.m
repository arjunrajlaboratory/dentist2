classdef scanObject < handle
    
    properties (Access = public)
        
        scanFile

        tilesTable
        channels
        dapiMasks
        %Not sure if we want to keep the rest as properties or make
        %variables instead. 
        scanMatrix
        rowTransformCoords 
        columnTransformCoords
        
    end
    
    methods
        
        function p = scanObject(scanFile, varargin) % 
            
            if nargin == 0
                files = dir('*.nd2');
                files = {files.name};
                p.scanFile = files{1};
                p.channels = d2utils.readND2Channels(p.scanFile);
                fprintf('New Table\n');
                p.tilesTable = cell2table(cell(0,5), 'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            elseif nargin == 1
                p.scanFile = scanFile;
                p.channels = d2utils.readND2Channels(p.scanFile);
                fprintf('New Table\n');
                p.tilesTable = cell2table(cell(0,5), 'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            elseif nargin == 2
                p.scanFile = scanFile;
                p.channels = d2utils.readND2Channels(p.scanFile);
                fprintf('Loading Table\n');
                p.tilesTable = readtable(varargin{1},'TextType','string');
            end
            
        end
        
        function p = loadTiles(p, scanMatrix, rowTransformCoords, columnTransformCoords)
            
            [height, width] = tileSize(p, varargin);
            
            for i = 1:numel(scanMatrix)
                [row,col] = find(scanMatrix == i);
                topCoords(i)  = col*columnTransformCoords(1) + row*rowTransformCoords(1);
                leftCoords(i) = row*rowTransformCoords(2) + col*columnTransformCoords(2);
            end

            topCoords = topCoords - min(topCoords) + 1;
            leftCoords = leftCoords - min(leftCoords) + 1;
            
            p.tilesTable = table((1:numel(scanMatrix))', topCoords', leftCoords', repmat(height, numel(scanMatrix),1), repmat(width, numel(scanMatrix),1), ...
                'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            
        end
        
         function [height, width] = tileSize(p)
             
             reader = bfGetReader(p.scanFile);
             omeMeta = reader.getMetadataStore();
             width = omeMeta.getPixelsSizeX(0).getValue();
             height = omeMeta.getPixelsSizeY(0).getValue();

         end
         
        %Set rowShift and columnShift
        function set.rowTransformCoords(p,coords)
            p.rowTransformCoords = coords;
        end
        
        function set.columnTransformCoords(p,coords)
            p.columnTransformCoords = coords;
        end
         
         
         %Mask DAPI
         
         
    end
         
        
end