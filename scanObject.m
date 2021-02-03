classdef scanObject < handle
    
    properties (Access = public)
        
        scanFile

        tilesTable
        channels
        %Not sure if we want to keep the rest as properties or make
        %variables in other scripts instead.
        scanMatrix
        dapiStitch
        stitchedChannels
        rowTransformCoords 
        columnTransformCoords
        dapiMask
        dapiMask2
        
    end
    
    methods
        
        function p = scanObject(scanFile, varargin) % 
            
            if nargin == 1
                files = dir('*.nd2');
                files = {files.name};
                p.scanFile = files{1};
                p.channels = d2utils.readND2Channels(p.scanFile);
                fprintf('New Table\n');
                p.tilesTable = cell2table(cell(0,5), 'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            elseif nargin == 2
                p.scanFile = scanFile;
                p.channels = d2utils.readND2Channels(p.scanFile);
                fprintf('New Table\n');
                p.tilesTable = cell2table(cell(0,5), 'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            elseif nargin == 3
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
            
            p.scanMatrix = scanMatrix;
            p.tilesTable = table((1:numel(scanMatrix))', topCoords', leftCoords', repmat(height, numel(scanMatrix),1), repmat(width, numel(scanMatrix),1), ...
                'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            
        end
        
%         function p = loadTiles2(p, scanMatrix)
%             
%             reader = bfGetReader(p.scanFile);
%             omeMeta = reader.getMetadataStore();
%             
%             width = omeMeta.getPixelsSizeX(0).getValue();
%             height = omeMeta.getPixelsSizeY(0).getValue();
%             
%              topCoords(1) = 1;
%              leftCoords(1) = 1;
%             pixelSize =  omeMeta.getPixelsPhysicalSizeX(0).value.doubleValue;
%             
%             tilesTmp = transpose(scanMatrix);
%             tiles = tilesTmp(:);
%             for i = 1:numel(tiles)
%                 
%                 %xpos_previous = omeMeta.getPlanePositionX(tiles(i-1)-1,0).value.doubleValue;
%                 xpos = omeMeta.getPlanePositionX(tiles(i)-1,0).value.doubleValue;
%                 %deltaX = (xpos-xpos_previous)/pixelSize
%                 topCoords(i) = xpos;
%                 
%                 %ypos_previous = omeMeta.getPlanePositionY(tiles(i-1)-1,0).value.doubleValue;
%                 ypos = omeMeta.getPlanePositionY(tiles(i)-1,0).value.doubleValue;
%                 %deltaY = (ypos-ypos_previous)/pixelSize
%                 leftCoords(i) = ypos;
%             end
%             reader.close()
%             
%             topCoords = topCoords/pixelSize;
%             topCoords = round(abs(topCoords - topCoords(1)));
%             topCoords = topCoords - min(topCoords) + 1;
%             
%             leftCoords = leftCoords/pixelSize;
%             leftCoords = round(abs(leftCoords - leftCoords(1)));
%             leftCoords = leftCoords - min(leftCoords) + 1;
%             
%             p.tilesTable2 = table((1:numel(scanMatrix))', topCoords', leftCoords', repmat(height, numel(tiles),1), repmat(width, numel(tiles),1), ...
%                 'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
%             
%         end
        
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
         
        %Stitch DAPI
        function p = stitchDAPI(p, scanMatrix)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = tileSize(p);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            channel = find(p.channels == 'dapi');
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            for i = 1:numel(tiles)
                
                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);
                                
                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = im2uint16(tmpPlane);
            end
            reader.close()
            p.dapiStitch = tmpStitch;
        end
        
        function p = maskDAPI(p, varargin) %vargargin option to specify sensitivity value (between 0 and 1) for adaptthresh.
            if ~isempty(p.dapiStitch)
                if nargin == 1
                    function_binarize = @(block_struct)...
                        imbinarize(scale(block_struct.data), adaptthresh(scale(block_struct.data), 0.05, 'ForegroundPolarity','bright'));
                    binary = blockproc(p.dapiStitch, [5000 5000], function_binarize, 'BorderSize', [200 200]);
                    p.dapiMask = binary;
                elseif nargin == 2
                    function_binarize = @(block_struct)...
                        imbinarize(scale(block_struct.data), adaptthresh(scale(block_struct.data), varargin{1}, 'ForegroundPolarity','bright'));
                    binary = blockproc(p.dapiStitch, [5000 5000], function_binarize, 'BorderSize', [200 200]);
                    p.dapiMask = binary;
                end
            else
                disp("dapiStitch is empty. Run stitchDAPI method and try again.")
            end
        end
        
        function p = stitchDAPImask(p, scanMatrix, varargin)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = tileSize(p);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'logical');
            channel = find(p.channels == 'dapi');
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            
            if nargin == 2
                s = 0.1;
            elseif nargin == 3
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
            
            p.dapiMask2 = tmpStitch;
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
        
        function p = stitchAllChannels(p, scanMatrix)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = tileSize(p);
            nChannels = numel(p.channels);
            stitches = cell(1,nChannels);
            reader = bfGetReader(p.scanFile);
            for i = 1:nChannels
                tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
                
                iPlane = reader.getIndex(0, i - 1, 0) + 1;
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    
                    tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                        tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = im2uint16(tmpPlane);
                end
                stitches{i} = tmpStitch;
            
            end
            reader.close()
            p.stitchedChannels = stitches;
            
        end
         
         
    end
         
        
end