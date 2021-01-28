classdef scanObject < handle
    
    properties (Access = public)
        
        scanFile

        tilesTable
        channels
        %Not sure if we want to keep the rest as properties or make
        %variables in other scripts instead. 
        scanMatrix
        rowTransformCoords 
        columnTransformCoords
        dapiMask
        
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
            
            [height, width] = tileSize(p);
            
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
         
        %Set scanMatrix, rowShift and columnShift
        function set.scanMatrix(p,matrix)
            p.scanMatrix = matrix;
        end
        
        function set.rowTransformCoords(p,coords)
            p.rowTransformCoords = coords;
        end
        
        function set.columnTransformCoords(p,coords)
            p.columnTransformCoords = coords;
        end
         
         %Mask DAPI
        function p = stitchDAPImask(p, scanMatrix)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = tileSize(p);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'logical');
            channel = find(p.channels == 'dapi');
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            for i = 1:numel(tiles)
                
                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);
                
                tileMask = d2utils.makeDAPImask(scale(tmpPlane));
                
                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tileMask;
            end
            reader.close()
            p.dapiMask = tmpStitch;
        end
        
        function tmpStitch = stitchChannel(p, scanMatrix, c)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = tileSize(p);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            channel = find(p.channels == c);
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            for i = 1:numel(tiles)
                
                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);
                                
                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = im2uint16(tmpPlane);
            end
            reader.close()
        end
         
         
    end
         
        
end