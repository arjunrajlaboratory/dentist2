function outRect = polyshapeBoundingBox(polyin)
    [xlim, ylim] = boundingbox(polyin);
    if isempty(xlim) %Assuming that if xlim is not empty, then ylim is not empty
        outRect = zeros(1,4);
       
    else
        start = max([xlim(1), ylim(1)], [1, 1]);
        sz = [diff(xlim), diff(ylim)];
        outRect = [floor(start), ceil(sz)];
    end
    
end