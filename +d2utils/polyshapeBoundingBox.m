function outRect = polyshapeBoundingBox(polyin)
    [xlim, ylim] = boundingbox(polyin);
    start = max([xlim(1), ylim(1)], [1, 1]);
    sz = [diff(xlim), diff(ylim)];
    outRect = [floor(start), ceil(sz)];
end