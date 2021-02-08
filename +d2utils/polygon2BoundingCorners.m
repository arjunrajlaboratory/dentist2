function outRect = polygon2BoundingCorners(polygon)
    mins = min(polygon, [], 1);
    maxs = max(polygon, [], 1)+1;
    outRect = single(combvec([mins(1), maxs(1)], [mins(2), maxs(2)])');
end