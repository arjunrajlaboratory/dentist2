function outMask = polyshape2mask(polyIn, shift, sz, varargin)
    p = inputParser;
    p.addRequired('polyIn', @(x)validateattributes(x,{'polyshape'}, {'scalar'}));
    p.addRequired('shift', @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}));
    p.addRequired('sz', @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}));
    p.addParameter('flip', false, @islogical);
    p.parse(polyIn, shift, sz, varargin{:});
    
    if p.Results.flip
        tmpVertices = fliplr(p.Results.polyIn.Vertices - p.Results.shift);
    else
        tmpVertices = p.Results.polyIn.Vertices - p.Results.shift;
    end
    sz = p.Results.sz;
    outMask = poly2mask(tmpVertices(:,1), tmpVertices(:,2), sz(1),  sz(2));
end