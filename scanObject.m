classdef scanObject < handle
    
    properties (Access = public)
        
        scanFile

        tilesTable
        tileSize
        channels
        scanMatrix
        scanDim
        snake
        startPos
        direction
        stitchDim
        dapiStitch
        smallDapiStitch
        stitchedScans
        smallStitchedScans
        resizeFactor = 4
        rowTransformCoords
        columnTransformCoords
        imRotation = 0;
%         dapiMask
%         dapiMask2
        
    end
    
    methods
        
        function p = scanObject(scanDim, varargin) % 
            n = inputParser;
            n.addRequired('scanDim', @(x)validateattributes(x,{'numeric'},{'size',[1 2]}));
            n.addParameter('scanFile', '', @ischar); 
            n.addParameter('tilesTable', '', @ischar); 
            n.addParameter('scanSummary', '', @ischar);
            
            n.parse(scanDim, varargin{:});
            
            p.scanDim = n.Results.scanDim;
            
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
                p.tilesTable = cell2table(cell(0,5), 'VariableNames', {'tileID', 'top', 'left', 'height', 'width'}); %Unncecessary? loadTiles creates a new table 
            else
                fprintf('Loading Table\n');
                opts = detectImportOptions(n.Results.tilesTable);
                opts = setvartype(opts, 'single'); %Probably unnecessary given the typical scan size. 
                p.tilesTable = readtable(n.Results.tilesTable, opts);
            end
            
            if isempty(n.Results.scanSummary)
                fprintf('New scan matrix\n');
                p.defaultScanParameters();
            else
                fprintf('Loading scan summary\n');
                p.parseScanSummary(n.Results.scanSummary);
            end
            
        end
        
        function defaultScanParameters(p)
            p.getTileSize();
            p.snake = true;
            p.direction = 'horizontal';
            p.startPos = 'top left';
            p.scanMatrix = d2utils.makeScanMatrix(p.scanDim);
            p.rowTransformCoords = [0 round(0.1*p.tileSize(1))]; %set default overlap to 10%
            p.columnTransformCoords = [round(0.1*p.tileSize(2)), 0];
        end
        
        function p = loadTiles(p) %should remove scanMatrix and TransformCoords arguments and leave them as properties that need to be set.    
            
            height = p.tileSize(1);
            width = p.tileSize(2);
            nTiles = numel(p.scanMatrix);
            topCoords = zeros(nTiles,1);
            leftCoords = zeros(nTiles,1);
            for i = 1:numel(nTiles)
                [row,col] = find(p.scanMatrix == i);
                topCoords(i)  = col*p.columnTransformCoords(1) + row*p.rowTransformCoords(1);
                leftCoords(i) = row*p.rowTransformCoords(2) + col*p.columnTransformCoords(2);
            end
            
            topCoords = topCoords - min(topCoords) + 1;
            leftCoords = leftCoords - min(leftCoords) + 1;
            
            p.tilesTable = table((1:nTiles)', topCoords', leftCoords', repmat(height, nTiles,1), repmat(width, nTiles,1), ...
                'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            
        end
          
         function p = getTileSize(p)
             reader = bfGetReader(p.scanFile);
             omeMeta = reader.getMetadataStore();
             p.tileSize = [omeMeta.getPixelsSizeY(0).getValue(), omeMeta.getPixelsSizeX(0).getValue()];
%              width = omeMeta.getPixelsSizeX(0).getValue();
%              height = omeMeta.getPixelsSizeY(0).getValue();
         end
         
         function outIm = getTileFromScan(p, tile, channel)
             reader = bfGetReader(p.scanFile);
             reader.setSeries(tile-1)
             channelIdx = find(ismember(p.channels, channel));  
             iPlane = reader.getIndex(0, channelIdx-1, 0) + 1;
             outIm  = bfGetPlane(reader, iPlane);
         end
         
        %Stitch DAPI
        function p = stitchDAPI(p)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
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
            p.stitchDim = size(tmpStitch);
        end
        
        function p = contrastDAPIstitch(p)
            function_scale =  @(block_struct) im2uint16(scale(block_struct.data));
            
            p.dapiStitch = blockproc(p.dapiStitch, [5000 5000], function_scale, 'BorderSize', [0 0]);
        end
        
        function p = contrastStitchedScans(p, percentiles, scaleFactor)  %Can modify this function for different default contrast
            function_contrast =  @(block_struct) im2uint16(d2utils.percentileScaleImage(block_struct.data, percentiles, scaleFactor));
            for i = 1:numel(p.stitchedScans.stitches)
                p.stitchedScans.stitches{i} = blockproc(p.stitchedScans.stitches{i}, [5000 5000], function_contrast, 'BorderSize', [0 0]);
            end
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
                
        function tmpStitch = stitchChannel(p, channel)  
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            c = find(ismember(p.channels,channel));
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, c - 1, 0) + 1;
            for i = 1:numel(tiles)

                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);

                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tmpPlane;
            end
            reader.close()
        end
        

        function p = stitchChannels(p, channels)
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
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
        
        function outIm = stitchTiles(p, rows, cols, channel, varargin) %Option to specify transformCoords - used with stitchingGUI
            if nargin == 4
                rowTransform = p.rowTransformCoords;
                colTransform = p.columnTransformCoords;
            elseif nargin > 4
                rowTransform = varargin{1};
                colTransform = varargin{2};
            end
            
            %Height and width may need to be switched for rotated images.
            %I'm not sure how tilesize changes when the image from Elements is rotated.  
            height = p.tileSize(1);
            width = p.tileSize(2);
            localScanMatrix = p.scanMatrix(rows(1):rows(2), cols(1):cols(2));
            tiles = transpose(localScanMatrix);
            tiles = tiles(:);
            topCoords = zeros(numel(tiles),1);
            leftCoords = zeros(numel(tiles),1);
            for i = 1:numel(tiles) 
                [r,c] = find(localScanMatrix == tiles(i));
                topCoords(i)  = c*colTransform(1) + r*rowTransform(1);
                leftCoords(i) = r*rowTransform(2) + c*colTransform(2);
            end
            topCoords = topCoords - min(topCoords) + 1;
            leftCoords = leftCoords - min(leftCoords) + 1;
            outIm = zeros(max(leftCoords)+height-1,max(topCoords)+width-1, 'uint16');
            reader = bfGetReader(p.scanFile);
            channelIdx = find(ismember(p.channels, channel));  
            iPlane = reader.getIndex(0, channelIdx - 1, 0) + 1;
            
            if logical(p.imRotation)
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                
                    outIm(leftCoords(ii):leftCoords(ii)+height-1, ...  
                        topCoords(ii):topCoords(ii)+width-1) = rot90(tmpPlane, p.imRotation);
                end
            else
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    outIm(leftCoords(ii):leftCoords(ii)+height-1, ...
                        topCoords(ii):topCoords(ii)+width-1) = tmpPlane;
                end
            end
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
        
        function savetilesTable(p, varargin)
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
        
        function parseScanSummary(p, varargin)
            if nargin == 1
                inFileName = 'scanSummary.txt';
            elseif nargin == 2
                inFileName = varargin{1};
            end
            
            inFileObj = fopen(inFileName);
            scanSummaryArray = textscan(inFileObj, '%s%q', 'Delimiter', '\t');
            fclose(inFileObj);
            scanSummaryTable = cell2table(scanSummaryArray{2}, 'RowNames', scanSummaryArray{1}');
            %scanSummaryTable = convertvars(scanSummaryTable, 'Variable', 'string');
            p.tileSize = str2num(cell2mat(scanSummaryTable{'imageSize',1}));
            p.snake = strcmp(scanSummaryTable{'snake',1}{:}, 'true');
            p.direction = scanSummaryTable{'scanDirection',1}{:};
            p.startPos = scanSummaryTable{'startPosition',1}{:};
            p.rowTransformCoords = str2num(cell2mat(scanSummaryTable{'rowTransform',1}));
            p.columnTransformCoords = str2num(cell2mat(scanSummaryTable{'columnTransform',1}));
            p.scanMatrix = d2utils.makeScanMatrix(p.scanDim, 'start', p.startPos, 'snake', p.snake,'direction', p.direction); 
        end
        
        function saveScanSummary(p, varargin)
            %NOTE! Avoid tabs ('\t') in your file name. Commas are OK. 
            if nargin == 1
                outFileName = 'scanSummary.txt';
            elseif nargin == 2
                outFileName = varargin{1};
            end
            
            outTableArray = {p.scanFile; num2str(p.scanDim); num2str(p.tileSize);...
                p.startPos; p.direction; string(p.snake); strjoin(p.channels);...
                num2str(p.rowTransformCoords); num2str(p.columnTransformCoords)};
            outTable = cell2table(outTableArray,'RowNames', {'scanFileName', 'scanDimensions', 'imageSize',...
                'startPosition', 'scanDirection', 'snake', 'channels', 'rowTransform', 'columnTransform'});
            writetable(outTable, outFileName, 'WriteRowNames', true, 'WriteVariableNames', false, 'QuoteStrings', false, 'Delimiter', '\t')  
        end
        
        function saveStitches(p)
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