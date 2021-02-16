function outRect = expandView2x(inRect, maxDim)
    startPos = max([1, 1], [inRect(1) - (inRect(3)/2)+1, inRect(2) - (inRect(4)/2)+1]); %Avoid rect outside range of scan. 
    startPos = min(maxDim- inRect(3:4)+1,round(startPos));
    outRect = [startPos, 2 * inRect(3:4)];
end