function channels = readND2Channels(fileND2)

    % Dictionary of commonly used channels and their names in Elements (in
    % Raj and Shaffer labs)
    channelMap = containers.Map({'Brightfield', 'DAPI', 'YFP', 'GFP', 'CY3', 'Cy3', 'A594', 'CY5', 'A647', '700', 'CY7','NIR'},...
                            {'trans'      , 'dapi', 'gfp', 'gfp', 'tmr', 'tmr', 'alexa', 'cy', 'cy'  , 'nir', 'nir','nir'});

    reader = bfGetReader(fileND2);
    omeMeta = reader.getMetadataStore();
    
    channels = cell(1, omeMeta.getPixelsSizeC(0).getValue());
    
    for i = 1:numel(channels)
        c = omeMeta.getChannelName(0, i-1);
        channels{i} = channelMap(c.toCharArray');
    end
    channels = string(channels);
end