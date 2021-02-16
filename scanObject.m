classdef scanObject < handle
    
    properties (Access = public)
        
        scanFile

        tilesTable
        channels
        scanMatrix
        dapiStitch
        smallDapiStitch
        stitchedScans
        smallStitchedScans
        resizeFactor = 4
%         dapiMask
%         dapiMask2
        
    end
    
    methods
        
        function p = scanObject(varargin) % 
            n = inputParser;
            n.addParameter('scanFile', '', @ischar); 
            n.addParameter('tilesTable', '', @ischar); 

            n.parse(varargin{:});
            
            if isempty(n.Results.scanFile)
                files = dir('*.nd2');
                files = {files.name};
                p.scanFile = files{1};
                p.channels = d2utils.readND2Channels(p.scanFile);
            else
                p.scanFile = n.Results.scanFile;
                p.channels = d2utils.readND2Channels(p.scanFile);
            end
            
            if isempty(n.Results.tilesTable)
                fprintf('New Table\n');
                p.tilesTable = cell2table(cell(0,5), 'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            else
                fprintf('Loading Table\n');
                p.tilesTable = readtable(n.Results.tilesTable,'TextType','string');
            end
            
        end
        
        function p = loadTiles(p, scanMatrix, rowTransformCoords, columnTransformCoords)
            
            [height, width] = p.tileSize();
            
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
          
         function [height, width] = tileSize(p)
             reader = bfGetReader(p.scanFile);
             omeMeta = reader.getMetadataStore();
             width = omeMeta.getPixelsSizeX(0).getValue();
             height = omeMeta.getPixelsSizeY(0).getValue();
         end
         
        %Stitch DAPI
        function p = stitchDAPI(p, scanMatrix)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = p.tileSize();
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            channel = find(ismember(p.channels, 'dapi'));
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
        
        function p = contrastDAPIstitch(p)
            function_scale =  @(block_struct) im2uint16(scale(block_struct.data));
            
            p.dapiStitch = blockproc(p.dapiStitch, [2000 2000], function_scale, 'BorderSize', [200 200]);
        end
        
        function p = resizeStitchedScans(p)
            %Can use blockproc as below but set 'BorderSize' to [0 0]. 
%             function_resize =  @(block_struct)... 
%             imresize(block_struct.data, 1/p.resizeFactor); 
            p.smallStitchedScans.labels = p.stitchedScans.labels;
            p.smallStitchedScans.stitches = cell(0, numel(p.stitchedScans.stitches));
            for i = 1:numel(p.stitchedScans.stitches)
                p.smallStitchedScans.stitches{i} = imresize(p.stitchedScans.stitches{i}, 1/p.resizeFactor);
            end
            p.smallDapiStitch = imresize(p.dapiStitch, 1/p.resizeFactor);
%             p.smallDapiStitch = blockproc(p.dapiStitch, [5000 5000], function_resize, 'BorderSize', [0 0]);
        end
        
%         function p = maskDAPI(p, varargin) %vargargin option to specify sensitivity value (between 0 and 1) for adaptthresh.
%             if ~isempty(p.dapiStitch)
%                 if nargin == 1
%                     function_binarize = @(block_struct)...
%                         imbinarize(scale(block_struct.data), adaptthresh(scale(block_struct.data), 0.05, 'ForegroundPolarity','bright'));
%                     binary = blockproc(p.dapiStitch, [5000 5000], function_binarize, 'BorderSize', [200 200]);
%                     p.dapiMask = binary;
%                 elseif nargin == 2
%                     function_binarize = @(block_struct)...
%                         imbinarize(scale(block_struct.data), adaptthresh(scale(block_struct.data), varargin{1}, 'ForegroundPolarity','bright'));
%                     binary = blockproc(p.dapiStitch, [5000 5000], function_binarize, 'BorderSize', [200 200]);
%                     p.dapiMask = binary;
%                 end
%             else
%                 disp("dapiStitch is empty. Run stitchDAPI method and try again.")
%             end
%         end
        
%         function p = stitchDAPImask(p, scanMatrix, varargin)  
%             tileTable = p.tilesTable;
%             tilesTmp = transpose(scanMatrix);
%             tiles = tilesTmp(:);
%             [height, width] = p.tileSize();
%             tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'logical');
%             channel = find(p.channels == 'dapi');
%             reader = bfGetReader(p.scanFile);
%             iPlane = reader.getIndex(0, channel - 1, 0) + 1;
%             
%             if nargin == 2
%                 s = 0.1;
%             elseif nargin == 3
%                 s = varargin{1};
%             end
%                 
%             for i = 1:numel(tiles)
%                 
%                 reader.setSeries(tiles(i)-1);
%                 tmpPlane  = bfGetPlane(reader, iPlane);
%                 
%                 tileMask = d2utils.makeDAPImask(scale(tmpPlane), 'sensitivity', s);
%                 
%                 tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
%                 tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tileMask;
%             end
%             reader.close()
%             
%             p.dapiMask2 = tmpStitch;
%         end
                
        function tmpStitch = stitchChannel(p, scanMatrix, channel)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = p.tileSize();
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            c = find(ismember(p.channels,channel));
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, c - 1, 0) + 1;
            for i = 1:numel(tiles)
                
                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);
                                
                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = im2uint16(tmpPlane);
            end
            reader.close()
        end
        

        function p = stitchChannels(p, channels)
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = p.tileSize();
            stitches = cell(1,numel(channels));
            channelIdx = find(ismember(p.channels, channels));
            reader = bfGetReader(p.scanFile);
            for i = 1:numel(channels)
                tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
                
                iPlane = reader.getIndex(0, channelIdx(i) - 1, 0) + 1;
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    
                    tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                        tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = im2uint16(tmpPlane);
                end
                stitches{i} = tmpStitch;
            
            end
            reader.close()
            p.stitchedScans.labels = channels;
            p.stitchedScans.stitches = stitches;
        end
        
        function outIm = getImageRect(p, channel, rect) %rect specified as [x y nrows ncols]
            channelIdx = ismember(p.stitchedScans.labels, channel);
            outIm = p.stitchedScans.stitches{channelIdx};
            outIm = outIm(rect(1):rect(1)+rect(3), rect(2):rect(2)+rect(4)); %Kinda ugly, would prefer using imcrop
        end
        
        function outIm = getSmallImageRect(p, channel, rect) %rect specified as [x y nrows ncols]
            channelIdx = ismember(p.smallStitchedScans.labels, channel);
            smallRect = round(rect/p.resizeFactor);
            smallRect(1:2) = max([1, 1], smallRect(1:2));
            outIm = p.stitchedScans.stitches{channelIdx};
            outIm = outIm(smallRect(1):smallRect(1)+smallRect(3), smallRect(2):smallRect(2)+smallRect(4)); %Kinda ugly, would prefer using imcrop
        end
        
        function outIm = getDapiImage(p, rect) %rect specified as [x y nrows ncols]
            outIm = p.dapiStitch(rect(1):rect(1)+rect(3), rect(2):rect(2)+rect(4));
        end
        
        function outIm = getSmallDapiImage(p, rect) %rect specified as [x y nrows ncols]
            smallRect = round(rect/p.resizeFactor);
            smallRect(1:2) = max([1, 1], smallRect(1:2));
            outIm = p.smallDapiStitch(smallRect(1):smallRect(1)+smallRect(3), smallRect(2):smallRect(2)+smallRect(4));
        end
        
        function [] = savetilesTable(p, varargin)
            if ~isempty(p.tilesTable)
                if nargin == 1
                    writetable(p.tilesTable, 'tilesTable.csv')
                elseif nargin ==2 
                    writetable(p.tilesTable, varargin{1})
                end
            else
                fprintf("tilesTable is empty. Run loadTiles and try again")
            end
        end
        
       function [] = saveStitches(p)
            if ~isempty(p.stitchedScans)
                temp = {p.stitchedScans};
                save('stitchedScans.mat', 'temp', '-v7.3')
            else
                fprintf("stitchedScans is empty. Run stitchChannels and try again")
            end
       end
       
       function p = loadStitches(p)
           tempStitches = load('stitchedScans.mat');
           p.stitchedScans.labels = tempStitches.temp{1}.labels;
           p.stitchedScans.stitches = tempStitches.temp{1}.stitches;
       end
        
    end
         
end