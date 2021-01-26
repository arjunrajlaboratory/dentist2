function tmpImage = ND2tileReader(inFile, tile, varargin)
    p = inputParser;

    p.addRequired('inFile', @ischar);
    p.addRequired('tile', @isnumeric);
    
    p.addParameter('channel', 1, @isnumeric);
    p.addParameter('Z', 1, @isnumeric);
    
    p.parse(inFile, tile, varargin{:});
    
    reader = bfGetReader(p.Results.inFile);
    
    reader.setSeries(tile-1);
    iPlane = reader.getIndex(p.Results.Z - 1, p.Results.channel - 1, 0) + 1;
    
    tmpImage  = bfGetPlane(reader, iPlane);
end
