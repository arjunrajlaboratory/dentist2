function outColors = expressionToColors(expressionVector)
    sampleColors = [0 0 1; 0 0.75 1; 0.165 1 0.189; 1 0.85 0; 1 0.66 0.04; 1 0 0];
    outColors = zeros(numel(expressionVector), 3);
    outColors(expressionVector== 0,:) = repmat([0.7 0.7 0.7], sum(expressionVector== 0),1); %set 0 to gray.
    positiveIntensities = expressionVector(expressionVector >0);
    outColors(expressionVector >0,:) = interp1(round(linspace(min(positiveIntensities), max(positiveIntensities), 6)), sampleColors, positiveIntensities);
    outColors = single(outColors);
end