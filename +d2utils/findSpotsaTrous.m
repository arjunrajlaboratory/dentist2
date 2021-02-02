function [X, Y, spotInt] = findSpotsaTrous(img, varargin)
    p = inputParser;
    p.addRequired('img', @(x)validateattributes(x,{'numeric'}, {'2d'}));
    p.addParameter('numLevels', 3, @(x)validateattributes(x,{'numeric'}, {'scaler', '>',0})); 
    p.addParameter('sigma', 2, @(x)validateattributes(x,{'numeric'}, {'scaler', '>',0})); 
    
    p.parse(img, varargin{:});
    
    img = p.Results.img;
    numLevels = p.Results.numLevels;
    sigma = p.Results.sigma;
    
    [aTrous, ~] = aTrousWaveletTransform(img,'numLevels',numLevels,'sigma',sigma);
    imgAT = sum(aTrous,3);
    bw = imregionalmax(imgAT);
    regionalMaxValues = imgAT(bw);
    regionalMaxIndices = find(bw);

    [regionalMaxValues,I] = sort(regionalMaxValues,'ascend');

    regionalMaxIndices = regionalMaxIndices(I);
    %Auto threshold
    [threshold] = imregmaxThresh(regionalMaxValues);
    if isempty(threshold)
        threshold = max(regionalMaxValues) + 1; %beyond max
    end

    spotInds = regionalMaxIndices(regionalMaxValues>(threshold/2));
    spotInt = regionalMaxValues(regionalMaxValues>(threshold/2));

    [X, Y] = ind2sub(size(bw),spotInds);  % convert 1D to 2D
end