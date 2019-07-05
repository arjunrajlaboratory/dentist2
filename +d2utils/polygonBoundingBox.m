function outRect = polygonBoundingBox(polygon)
    outRect(1) = min(polygon(:,1));
    outRect(2) = min(polygon(:,2));
    outRect(3) = max(polygon(:,1))-min(polygon(:,1))+1;
    outRect(4) = max(polygon(:,2))-min(polygon(:,2))+1;
end