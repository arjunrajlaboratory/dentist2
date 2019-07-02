classdef spotTable < handle
    
    properties (Access = public)
        
        spots
        masks % per channel
        thresholds % per channel
        radius = 100;
        spotChannels
        
        theFilter
        percentileToKeep = 98;
        
    end
    
    methods
        
        function p = spotTable(fileTable)
            p.spots = [];
            p.masks = [];
            p.thresholds = [];
            channels = fileTable.channels;
            p.spotChannels = channels(~ismember(channels,["dapi","trans"]));
            
        end
        
        function findSpots(p,fileTable)
            
            p.spots = cell2table(cell(0,8), 'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestCellID', 'status', 'channel', 'fileID'});
            
            %spotID = [];
            x = [];
            y = [];
            intensity = [];
            %nearestCellID = [];
            %status = [];
            channel = [];
            fileID = [];
            
            p.theFilter = -fspecial('log',20,2);
            for i = 1:length(p.spotChannels)
                currChannel = p.spotChannels(i);
                fprintf('Channel: %s\n',currChannel);
                files = fileTable.files(fileTable.files.channel == currChannel,:);
                for j = 1:size(files,1)
                    fprintf('File: %s\n',files.fileName(j));
                    im = imread(files.fileName(j));
                    % Convert to single to "smooth" things so that
                    % imregional max will isolate a particular pixel
                    filt = imfilter(im2single(im),p.theFilter,'replicate');
                    irm = imregionalmax(filt);
                    tempSpots = filt(irm)';
                    thresh = prctile(tempSpots,p.percentileToKeep);
                    %idx = tempSpots > thresh;
                    filt = filt.*single(irm);
                    filt(filt < thresh) = 0;
                    filt = im2uint16(filt);
                    idx = filt > 0;
                    intensities = filt(idx);
                    [tempX,tempY] = ind2sub(size(filt),find(idx));
                    x = [x ; tempX];
                    y = [y ; tempY];
                    intensity = [intensity ; intensities];
                    channel = [channel ; repmat(currChannel,length(x),1)];
                    fileID = [fileID ; repmat(files.fileName(j),length(x),1)];
                end
            end
            
            4+4
            
            
            
        end

        
    end
    
end 