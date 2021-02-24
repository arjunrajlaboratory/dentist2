classdef scanObject < handle
    
    properties (Access = public)
        
        scanFile
        tilesTableName

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
        %Including imRotation in case the image acquired by Elements is flipped or somehow
        %rotated. Consider deleting. Also, height and width 
        %may need to be switched for rotated images. I'm not sure if 
        %tilesize changes when the image from Elements is rotated.
        imRotation = 0;  
    end
    
    methods
        
        function p = scanObject(varargin) % 
            n = inputParser;
            n.addParameter('scanSummary', '', @ischar);
            n.addParameter('scanFile', '', @ischar); 
            n.addParameter('scanDim', [0,0], @(x)validateattributes(x,{'numeric'},{'size',[1 2]}));
            
            n.parse(varargin{:});

            if ~isempty(n.Results.scanSummary)
                    fprintf('Loading %s\n', n.Results.scanSummary);
                    p.loadScanSummary(n.Results.scanSummary);
            else
                if ~isempty(n.Results.scanFile) && all(n.Results.scanDim)
                    fprintf('Loading file %s\n', n.Results.scanFile)
                    p.scanFile = n.Results.scanFile;
                    p.channels = d2utils.readND2Channels(p.scanFile);
                    p.scanDim = n.Results.scanDim;
                    p.defaultScanParameters()
                else
                    fprintf('Unable to create the scan object.\nPlease specify a scan summary file (e.g. scanSummary.txt) or the scan file name and scan dimensions.\n')
                    return
                    
                end
            end
        end
        
        function defaultScanParameters(p)
            p.getTileSize();
            p.snake = true;
            p.direction = 'horizontal';
            p.startPos = 'top left';
            p.scanMatrix = d2utils.makeScanMatrix(p.scanDim);
%             p.rowTransformCoords = [0 round(0.1*p.tileSize(1))]; %set default overlap to 10%
%             p.columnTransformCoords = [round(0.1*p.tileSize(2)), 0];
        end
        
        function p = loadTiles(p)  
            
            height = p.tileSize(1);
            width = p.tileSize(2);
            nTiles = numel(p.scanMatrix);
            topCoords = zeros(nTiles,1);
            leftCoords = zeros(nTiles,1);
            for i = 1:nTiles
                [row,col] = find(p.scanMatrix == i);
                topCoords(i)  = col*p.columnTransformCoords(1) + row*p.rowTransformCoords(1);
                leftCoords(i) = row*p.rowTransformCoords(2) + col*p.columnTransformCoords(2);
            end
            
            topCoords = topCoords - min(topCoords) + 1;
            leftCoords = leftCoords - min(leftCoords) + 1;
            p.tilesTable = table((1:nTiles)', topCoords, leftCoords, repmat(height, nTiles,1), repmat(width, nTiles,1), ...
                'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            p.stitchDim = [max(leftCoords)+p.tileSize(1)-1, max(topCoords)+p.tileSize(2)-1];
        end
          
         function p = getTileSize(p)
             reader = bfGetReader(p.scanFile);
             omeMeta = reader.getMetadataStore();
             p.tileSize = [omeMeta.getPixelsSizeY(0).getValue(), omeMeta.getPixelsSizeX(0).getValue()];
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
            
            if logical(p.imRotation) %In case the image acquired by Elements is flipped or somehow rotated. Consider deleting. 
                for i = 1:numel(tiles)
                    reader.setSeries(tiles(i)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                        tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = rot90(tmpPlane, p.imRotation);
                end
            else
                for i = 1:numel(tiles)
                    reader.setSeries(tiles(i)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                        tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tmpPlane;
                end
            end
            
            reader.close()
            p.dapiStitch = tmpStitch;
%             p.stitchDim = size(tmpStitch);
        end
        
        function p = contrastDAPIstitch(p)
            function_scale =  @(block_struct) im2uint16(scale(block_struct.data));
            
            p.dapiStitch = blockproc(p.dapiStitch, [5000 5000], function_scale, 'BorderSize', [0 0], 'UseParallel', true);
        end
        
        function p = contrastStitchedScans(p, percentiles, scaleFactor)  %Can modify this function for different default contrast
            function_contrast =  @(block_struct) im2uint16(d2utils.percentileScaleImage(block_struct.data, percentiles, scaleFactor));
            for i = 1:numel(p.stitchedScans.stitches)
                p.stitchedScans.stitches{i} = blockproc(p.stitchedScans.stitches{i}, [5000 5000], function_contrast, 'BorderSize', [0 0], 'UseParallel', true);
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
        
        function tmpStitch = stitchChannel(p, channel) %Maybe unnecessary
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            c = find(ismember(p.channels,channel));
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, c - 1, 0) + 1;
            if logical(p.imRotation) %In case the image acquired by Elements is flipped or somehow rotated. Consider deleting. 
                for i = 1:numel(tiles)
                    reader.setSeries(tiles(i)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);

                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                    tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = rot90(tmpPlane, p.imRotation);
                end
            else
                for i = 1:numel(tiles)
                    reader.setSeries(tiles(i)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);

                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                    tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tmpPlane;
                end
            end
            reader.close()
        end

        function p = stitchChannels(p, varargin)
            if nargin == 1
                channelsToStitch = p.channels(~ismember(p.channels,{'dapi','trans'}));
            elseif nargin == 2
                channelsToStitch = varargin{1};
            end
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
            stitches = cell(1,numel(channelsToStitch));
            channelIdx = find(ismember(p.channels, channelsToStitch));
            reader = bfGetReader(p.scanFile);
            for i = 1:numel(channelsToStitch)
                tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
                iPlane = reader.getIndex(0, channelIdx(i) - 1, 0) + 1;
                if logical(p.imRotation) %In case the image acquired by Elements is flipped or somehow rotated. Consider deleting.
                    for ii = 1:numel(tiles)
                        reader.setSeries(tiles(ii)-1);
                        tmpPlane  = bfGetPlane(reader, iPlane);
                    
                        tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                            tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = rot90(tmpPlane, p.imRotation);
                    end
                    stitches{i} = tmpStitch;
                else
                    for ii = 1:numel(tiles)
                        reader.setSeries(tiles(ii)-1);
                        tmpPlane  = bfGetPlane(reader, iPlane);
                    
                        tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                            tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = tmpPlane;
                    end
                    stitches{i} = tmpStitch;
                end
            end
            reader.close()
            p.stitchedScans.labels = channelsToStitch;
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
            outIm = p.smallStitchedScans.stitches{channelIdx};
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
        
        function saveTilesTable(p, varargin)
            if ~isempty(p.tilesTable) %Not sure we need the option to specify alternative filename
                if nargin == 1
                    p.tilesTableName = 'tilesTable.csv';
                elseif nargin ==2 
                    p.tilesTableName = varargin{1};
                end
                writetable(p.tilesTable, p.tilesTableName)
            else
                fprintf("tilesTable is empty. Run loadTiles and try again")
            end
        end
        
        function loadTilesTable(p, varargin)
            if nargin == 1
                inFileName = 'tilesTable.csv';
            elseif nargin == 2
                inFileName = varargin{1};
            end
            fprintf('Loading %s\n', inFileName);
            opts = detectImportOptions(inFileName);
            opts = setvartype(opts, 'single'); %Probably unnecessary given the typical scan size.
            p.tilesTable = readtable(inFileName, opts);
        end
        
        function loadScanSummary(p, varargin)
            if nargin == 1
                inFileName = 'scanSummary.txt';
            elseif nargin == 2
                inFileName = varargin{1};
            end
            
            scanSummaryTable = d2utils.parseScanSummary(inFileName);
            p.scanFile = scanSummaryTable{'scanFileName',1}{:};%Could possibly check that the scanFile tilesTable exist
            p.tilesTableName = scanSummaryTable{'tilesTableName',1}{:};
            p.scanDim = str2double(split(scanSummaryTable{'scanDimensions',1})');
            p.tileSize = str2double(split(scanSummaryTable{'imageSize',1})');
            p.stitchDim = str2double(split(scanSummaryTable{'stitchDimensions',1})');
            p.snake = strcmp(scanSummaryTable{'snake',1}{:}, 'true');
            p.direction = scanSummaryTable{'scanDirection',1}{:};
            p.startPos = scanSummaryTable{'startPosition',1}{:};
            p.rowTransformCoords = str2double(split(scanSummaryTable{'rowTransform',1})');
            p.columnTransformCoords = str2double(split(scanSummaryTable{'columnTransform',1})');
            
            p.scanMatrix = d2utils.makeScanMatrix(p.scanDim, 'start', p.startPos, 'snake', p.snake,'direction', p.direction); 
            p.channels = d2utils.readND2Channels(p.scanFile);
            %Load tiles table.
            if isfile(p.tilesTableName)
                p.loadTilesTable(p.tilesTableName);
            end 
        end
        
        function saveScanSummary(p, varargin) 
            %NOTE! Avoid tabs ('\t') in your file name. Commas are OK. 
            %Should maybe check that all the necessary properties are not
            %empty. 
            if nargin == 1 %Not sure we need the option to specify alternative filename
                outFileName = 'scanSummary.txt';
            elseif nargin == 2
                outFileName = varargin{1};
            end
            
            outTableArray = {p.scanFile; p.tilesTableName; num2str(p.scanDim); num2str(p.tileSize); num2str(p.stitchDim);...
                p.startPos; p.direction; string(p.snake); strjoin(p.channels);...
                num2str(p.rowTransformCoords); num2str(p.columnTransformCoords)};
            outTable = cell2table(outTableArray,'RowNames', {'scanFileName', 'tilesTableName', 'scanDimensions', 'imageSize', 'stitchDimensions',...
                'startPosition', 'scanDirection', 'snake', 'channels', 'rowTransform', 'columnTransform'});
            writetable(outTable, outFileName, 'WriteRowNames', true, 'WriteVariableNames', false, 'QuoteStrings', false, 'Delimiter', '\t')  
        end
                
        function saveStitches(p, varargin) 
            %For the time being, saving to .mat files. 
            %At some point it may be worth modifying this method to save
            %bigTiff files for compatibility with other software. 
            %Will need to update loadStitches as well. 
            if nargin == 1 %Not sure we need the option to specify alternative filename
                outFileName = 'stitchedScans.mat';
            elseif nargin == 2
                outFileName = sprintf('%s.mat', varargin{1});
            end
            
            if ~isempty(p.dapiStitch) && ~isempty(p.stitchedScans) 
                d2StitchedScans = {p.dapiStitch, p.stitchedScans};
                save(outFileName, 'd2StitchedScans', '-v7.3')
            else
                if isempty(p.dapiStitch)
                    fprintf("dapiStitch is empty. Run stitchDAPI and try again\n")
                end
                if isempty(p.stitchedScans)
                    fprintf("stitchedScans is empty. Run stitchChannels and try again\n")
                end
            end
       end
       
       function p = loadStitches(p)
           load('stitchedScans.mat', 'd2StitchedScans');
           p.dapiStitch = d2StitchedScans{1};
           p.stitchedScans = d2StitchedScans{2};
%            clear d2StitchedScans
       end
        
    end
         
end