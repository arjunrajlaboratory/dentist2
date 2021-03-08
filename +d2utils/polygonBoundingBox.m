function outRect = polygonBoundingBox(polygon)
    mins = min(polygon);
    maxs = max(polygon);
    outRect = [mins, maxs-mins+1];
end