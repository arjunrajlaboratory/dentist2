function outRect = polygonBoundingBox(polygon)
    mins = min(polygon);
    maxs = max(polygon);
    outRect = [mins, maxs-mins+1];
%     outRect(1) = min(polygon(:,1));
%     outRect(2) = min(polygon(:,2));
%     outRect(3) = max(polygon(:,1))-min(polygon(:,1))+1;
%     outRect(4) = max(polygon(:,2))-min(polygon(:,2))+1;
end