function [] = splitStitchedScan(varargin)
    p = inputParser;
    p.addParameter('scanIm', zeros(2), @(x)validateattributes(x,{'numeric'}, {'2d'}));
    p.addParameter('scanFile', '', @ischar);
    p.addParameter('channel', 'dapi', @ischar);
    p.addParameter('tileSize', [1000 1000], @(x)validateattributes(x,{'numeric'}, {'size', [1, 2]}));    
    p.addParameter('overlap', [50 50], @(x)validateattributes(x,{'numeric'}, {'size', [1, 2]}));
    p.addParameter('resize', 1, @(x)validateattributes(x,{'uint8'}, {'scalar', '>', 0}));
    p.addParameter('contrastPercentiles', [0 0], @(x)validateattributes(x,{'numeric'}, {'size', [1, 2], '>', 0, '<=', 100}));
    p.addParameter('contrastScaleFactor', [1 1], @(x)validateattributes(x,{'numeric'}, {'size', [1, 2], '>', 0}));
    p.addParameter('outDir', 'splitTiles', @ischar);    
   
    p.parse(varargin{:});
    scanIm = p.Results.scanIm;
    scanFile = p.Results.scanFile;
    tileSize = p.Results.tileSize;
    overlap = p.Results.overlap;
    resize = p.Results.resize;
    contrastPercentiles = p.Results.contrastPercentiles;
    contrastScaleFactor = p.Results.contrastPercentiles;
    outDir = p.Results.outDir;
    
    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end
    
    if isempty(scanFile) && ~any(scanIm, 'all')
        disp('Please input an image or specify an image .nd2 file. Do not input both.')
        return
    elseif ~isempty(scanFile) && any(scanIm, 'all')
        disp('Please input an image or specify an image .nd2 file. Do not input both.')
        return
    elseif ~isempty(scanFile) && ~any(scanIm, 'all')
        channels = d2utils.readND2Channels(scanFile);
        reader = bfGetReader(scanFile);
        eader.setSeries(0); %Will only load the first scan 
        channelIdx = find(ismember(channels, p.Results.channel));
        iPlane = reader.getIndex(0, channelIdx - 1, 0) + 1; %Will only read the first z-plane
        scanIm = bfGetPlane(reader, iPlane);        
    end
    
    rowOffset = tileSize(1) - overlap(1);
    colOffset = tileSize(2) - overlap(2);
    
    scanDim = size(scanIm);
    rowStart = 1:rowOffset:scanDim(1)-1
    rowEnd = rowStart+tileSize(1)-1
    rowEnd(end) = scanDim(1)
    
    colStart = 1:colOffset:scanDim(2)-1
    colEnd = colStart+tileSize(2)-1
    colEnd(end) = scanDim(2)
    if all(contrastPercentiles) %If you're contrasting the images. 
        for i = 1:numel(rowStart)
            for ii = 1:numel(colStart)
                imwrite(imresize(d2utils.percentileAdjustImage(scanIm(rowStart(i):rowEnd(i), colStart(ii):colEnd(ii)), contrastPercentiles, contrastScaleFactor), 1/resize), sprintf('%s/tile%d_%d.png', outDir, i, ii))
            end
        end
    else
        for i = 1:numel(rowStart)
            for ii = 1:numel(colStart)
                imwrite(imresize(scanIm(rowStart(i):rowEnd(i), colStart(ii):colEnd(ii)), 1/resize), sprintf('%s/tile%d_%d.png', outDir, i, ii))
            end
        end
    end
    
end