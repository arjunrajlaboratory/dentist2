classdef d2MainAxesController < handle
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        scanDim
        viewObj
        threshZoom = false
        viewRect %Image coords not axes coords
        panRect
        cellViewRadius = 500
        zoomROI
        zoomStart
        zoomRect
        ThreshAxisMin
        ThreshAxisMax
        fixedZoom = false
        panStart
        channelIdx
        zoomMode = true
        scatterH %Not sure if we need these handle. 
        imageH
        spotScatterH
        thumbPlotH
        thumbRectH
        zoomThreshStartH
        zoomThreshEndH
        zoomThreshRectH
        nucleiScatterH
        maskH
        paletteIdx = 1 
        imagesInView
        dapiInView
        
        resizedImg %Not sure if we want to save these or regenerate
        resizedDapi
        
    end
    
    methods
        function p = d2MainAxesController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
            p.scanDim = size(p.scanObj.dapiStitch);
            p.viewRect = [1, 1, min(p.scanDim, [25000 25000])];
            p.startup()
        end
        
        function startup(p)
            p.channelIdx = p.viewObj.channelPopup.Value;
            p.spotTable.makeIntensitiesToPlot();
            p.plotIntensityHistogram()
            p.plotIntensityThreshold()
            p.plotScatterMain();
            p.plotThumbnail();
            p.updateImageInView();
            %p.updateMainAxes()
        end
        
        function p = changeChannel(p, ~, ~)
            p.channelIdx = p.viewObj.channelPopup.Value;
            %Update centroid listbox
            p.updateCentroidListView();
            p.updateMainAxes();
            p.plotIntensityHistogram();
        end
        
        function changeColormap(p, ~, ~)
            p.paletteIdx = p.viewObj.colormapPopup.Value;
            p.spotTable.updateExpressionColors(p.paletteIdx);
            p.updateMainAxes();
        end
        
        function p = updateMainAxes(p, ~, ~)
            if isvalid(p.viewObj.mainAxes.Children)
                delete(get(p.viewObj.mainAxes, 'Children'));
            end
            
            if logical(p.viewObj.scatterCheckBox.Value)
                p.plotScatterMain();
            else
                p.showImage();
                p.overlaySpots();
                p.overlayNuclei();
                p.overlayMasks();
            end
        end
        
        function updateCentroidListView(p)
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
        end
        
        function p = centroidSelected(p, ~, ~)
            if strcmp(get(p.viewObj.figHandle, 'SelectionType'), 'open') %Respond to double mouse click
                cellIdx = get(p.viewObj.centroidList, 'Value');
                cellPos = p.spotTable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanDim);
                set(p.viewObj.scatterCheckBox, 'Value', 0)
                p.updateImageInView();
                p.updateMainAxes();
                p.overlayThumbnailRect();
            end
        end
        
        function plotIntensityHistogram(p)
            %channelIdx = p.viewObj.channelPopup.Value;
            if ~isempty(p.viewObj.histogramLineH)
                delete(p.viewObj.histogramLineH)
            end
            intensities = p.spotTable.intensitiesToPlot{p.channelIdx};
            logRank = log(numel(intensities):-1:1);
            p.viewObj.histogramLineH = line(p.viewObj.threshAxes, intensities, logRank, ...
                'HitTest', 'off', ...
                'Color', 'k');
            
            p.ThreshAxisMin = intensities(1);
            p.ThreshAxisMax = intensities(end) * 1.05;
            
            set(p.viewObj.threshAxes, 'XLim', [p.ThreshAxisMin p.ThreshAxisMax]);
            
            yaxismax = logRank(1)*1.1;
          
            set(p.viewObj.threshAxes, 'YLim', [0 yaxismax]);
        end
        
        function plotIntensityThreshold(p)
            %channelIdx = p.viewObj.channelPopup.Value;
            if ~isempty(p.viewObj.thresholdLineH) && ishandle(p.viewObj.thresholdLineH)
                delete(p.viewObj.thresholdLineH)
            end
            if isempty(p.spotTable.thresholds)
                p.spotTable.defaultThresholds();
            end
            threshold = p.spotTable.thresholds{p.channelIdx};
            yaxis = get(p.viewObj.threshAxes, 'Ylim');
            p.viewObj.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis,...
                'Color', 'b', 'HitTest', 'off');
            p.viewObj.threshValue.String = num2str(threshold);
        end
        
        function plotScatterMain(p)
            %Could probably save on time by adding colormap to
            %centroidTable in the spotTableObject.
            centroidsInView = p.spotTable.centroidTableInRect(p.channelIdx, p.viewRect);
            centroidsInView = flipud(centroidsInView); %In order to plot cells with 0 expression in the back
            set(p.viewObj.mainAxes, 'XLim', [p.viewRect(2),  p.viewRect(2)+p.viewRect(4)])
            set(p.viewObj.mainAxes, 'YLim', [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)])
            hold(p.viewObj.mainAxes, 'on')
            %Note that it is slightly slower specifying RGB colors with expression_color
            %instead of changing the colormap on the axes however, specifying RGB colors
            %ensures that the colors don't change with the view 
            p.scatterH = scatter(centroidsInView.y, centroidsInView.x,...
                20, centroidsInView.expression_color, 'filled',...
                'Parent', p.viewObj.mainAxes, 'HitTest','off');
            hold(p.viewObj.mainAxes, 'off')
            %colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
        end
        
        function plotThumbnail(p)
            set(p.viewObj.thumbAxes, 'XLim', [1,  p.scanDim(2)])
            set(p.viewObj.thumbAxes, 'YLim', [1,  p.scanDim(1)])
            hold(p.viewObj.thumbAxes, 'on')
            %Note that it is slightly slower specifying RGB colors with expression_color
            %instead of changing the colormap on the axes however, specifying RGB colors
            %ensures that the colors don't change with the view 
            p.thumbPlotH = binscatter(p.spotTable.centroidLists{p.channelIdx}.y, p.spotTable.centroidLists{p.channelIdx}.x, round(min(p.scanDim/1000)),...
                'Parent', p.viewObj.thumbAxes, 'HitTest','off');
            p.overlayThumbnailRect();
            hold(p.viewObj.thumbAxes, 'off')
            
        end
        
        function overlayThumbnailRect(p)
            if any(p.viewRect(3:4) < p.scanDim) && all(p.viewRect(3:4) > 4) %No overlay if too zoomed out or in
                if isempty(p.thumbRectH) || ~isvalid(p.thumbRectH)
                    p.thumbRectH = rectangle(...
                    'Position', d2utils.rotateRectROI(p.viewRect), ...
                    'EdgeColor', 'b', 'LineWidth', 2,...
                    'Parent', p.viewObj.thumbAxes,...
                    'Hittest', 'off');
                else
                    set(p.thumbRectH, 'Position', d2utils.rotateRectROI(p.viewRect))
                end
            else
                delete(p.thumbRectH)
            end
        end
        
        function zoomInPressed(p, ~, ~)
            p.zoomMode = true;
            %iptPointerManager(p.viewObj.figHandle, 'disable');
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
            %Below are old functions that draw rectangle using matlab
            %ROI instead of mainAxesZoom fcn. Probably slightly more efficient but doesn't allow using
            %the console until roi completed. 
            %set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDownZoomMode}) %need a zoomMode that doesn't use 'normal' clicks
%             while p.zoomMode
%                 if p.fixedZoom
%                     p.zoomROI = drawrectangle(p.viewObj.mainAxes, 'Color', 'r', 'FaceAlpha', 0, 'FixedAspectRatio', true);
%                 else 
%                     p.zoomROI = drawrectangle(p.viewObj.mainAxes, 'Color', 'r', 'FaceAlpha', 0);
%                 end
%                 
%                 if isempty(p.zoomROI.Position) %if escape pressed
%                     break
%                 elseif all(p.zoomROI.Position(3:4) > 4) %minimum view size
%                     p.viewRect  = d2utils.coordToPixelRect(round(p.zoomROI.Position));
%                     delete(p.zoomROI) %Possibly unnecessary since Children deleted in p.updateMainAxes(); 
%                     p.updateImageInView();
%                     p.updateMainAxes();
%                     p.fixedZoom = false;
%                 end
%             end
        end
    
        function figWindowDown(p, ~, ~)
            if d2utils.pointInSideViewRect(p.viewRect, get(p.viewObj.mainAxes, 'CurrentPoint'))
                switch(p.getSelectionType)
                    %'normal' click was used for drawing rectangle before using
                    %rectangle roi
                    case 'normal'
                        if p.zoomMode
                            p.zoomStart = get(p.viewObj.mainAxes, 'CurrentPoint');
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopDragFcn});
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.mainAxesZoom});
                        else
                            p.panStart =  get(p.viewObj.mainAxes, 'CurrentPoint');
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopPan})
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.panView});
                        end
                    case 'open'
                        p.viewRect = [1, 1, min(p.scanDim, [25000 25000])];
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.overlayThumbnailRect();
                    case 'alt'
                        p.viewRect = d2utils.expandView2x(p.viewRect, p.scanDim);
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.overlayThumbnailRect();
                end
            end
        end
        
        function mainAxesZoom(p, ~, ~)
            if ~isempty(p.viewObj.zoomH)
                delete(p.viewObj.zoomH)
            end
            
            if p.fixedZoom
                p.getSelectedRectangleCoordsFixed()
            else
                p.getSelectedRectangleCoords()
            end
            
            if all(p.zoomRect(3:4) > 4) %min view size
                p.viewObj.zoomH = rectangle(...
                    'Position', p.zoomRect, ...
                    'EdgeColor', 'r', ...
                    'Parent', p.viewObj.mainAxes,...
                    'Hittest', 'off');
            end

        end
        
        function stopDragFcn(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            if ~isempty(p.zoomRect) && all(p.zoomRect > 0)
                p.viewRect  = d2utils.coordToPixelRect(p.zoomRect); % Should notify event "viewChange"
                %Update imagesInView
                p.updateImageInView();
                %Update main axes
                p.updateMainAxes();
                p.overlayThumbnailRect();
                p.zoomRect = [];
                p.fixedZoom = false;
            end
            delete(p.viewObj.zoomH)
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
        
        function getSelectedRectangleCoords(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)),...
                abs(p.zoomStart(1,1:2) - currentPoint(1,1:2))]);
        end
        
        function getSelectedRectangleCoordsFixed(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            sz = max(abs(p.zoomStart(1,1:2) - currentPoint(1,1:2)));
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)), [sz sz]]);
        end
        
        function panViewPressed(p, ~, ~)
            if any(p.viewRect(3:4) < p.scanDim) %No need to pan if view is of entire scan
                p.zoomMode = false; 
                %iptPointerManager(p.viewObj.figHandle, 'enable');
                %iptSetPointerBehavior(p.viewObj.mainAxes, @(hfig, currentPoint) set(hfig, 'Pointer', 'hand'));
                set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
            else
                p.zoomMode = false; 
                set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            end
        end
                
        function stopPan(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '')
        end
        
        function panView(p, ~, ~)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            displacement = p.panStart(1,1:2) - currentPoint(1,1:2);
            p.viewRect = d2utils.updateViewPanning(p.viewRect, displacement, p.scanDim);
            p.updateImageInView;
            p.updateMainAxes;
            p.overlayThumbnailRect();
        end
        
        function thumbAxesButtonDown(p, ~, ~)
            if any(p.viewRect(3:4) < p.scanDim)
                currentPoint = get(p.viewObj.thumbAxes, 'CurrentPoint');
                switch(p.getSelectionType)
                    case 'open'
                        newRect = d2utils.getRectAroundPoint(currentPoint(1,2:-1:1), p.viewRect(3), p.viewRect(4), p.scanDim);
                        p.viewRect = newRect;
                        p.updateImageInView;
                        p.updateMainAxes;
                        p.overlayThumbnailRect();
                    case 'normal'
                        if d2utils.pointInSideViewRect(p.viewRect, currentPoint)
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopThumbDrag});
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.dragThumb});
                        end
                end
            end
        end
        
        function dragThumb(p, ~, ~)
            currentPoint = get(p.viewObj.thumbAxes, 'CurrentPoint');
            newRect = d2utils.getRectAroundPoint(currentPoint(1,2:-1:1), p.viewRect(3), p.viewRect(4), p.scanDim);
            p.viewRect = newRect;
            p.overlayThumbnailRect();
            p.updateImageInView;
            p.updateMainAxes;
        end
        
        function stopThumbDrag(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
        
        function threshZoomButtonDown(p, ~, ~)
            p.threshZoom = true;
        end
        
        function thresholdButtonDown(p, ~, ~)
            %set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopThresholdDrag});
            %set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.thresholdDrag});
            switch(p.getSelectionType)
                case 'normal'
                    currentPoint = get(p.viewObj.threshAxes, 'CurrentPoint');
                    if ~logical(p.threshZoom)
                        if abs(currentPoint(1,1) - p.viewObj.thresholdLineH.XData(1)) < 150
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopThreshDrag});
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.dragThresh});
                        end
                    else
                        yaxis = get(p.viewObj.threshAxes, 'Ylim');
                        p.zoomThreshRectH = patch('YData', [yaxis, fliplr(yaxis)], 'XData', [currentPoint(1:2,1)', currentPoint(1:2,1)'],...
                            'FaceColor', 'red', 'FaceAlpha', 0.2,...
                            'Parent', p.viewObj.threshAxes, 'Hittest', 'off');
                        set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopzoomThresh});
                        set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.zoomThresh});
                    end
                case 'open'
                    set(p.viewObj.threshAxes, 'XLim', [p.ThreshAxisMin p.ThreshAxisMax]);
                    p.threshZoom = false;
            end
        end
        
        function dragThresh(p, ~, ~)
            currentPoint = get(p.viewObj.threshAxes, 'CurrentPoint');
            newThresh = max(currentPoint(1,1), 1);
            newThresh = min(newThresh, p.viewObj.threshAxes.XLim(2));
            set(p.viewObj.thresholdLineH, 'XData', [newThresh, newThresh])
            p.viewObj.threshValue.String = num2str(round(newThresh));
        end
        
        function stopThreshDrag(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            %Should notify threshold change
            newThresh = str2double(p.viewObj.threshValue.String);
            p.spotTable.setThreshold(p.spotTable.spotChannels{p.channelIdx}, newThresh);
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
            p.updateMainAxes();
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
        
        function zoomThresh(p, ~, ~)
            currentPoint = get(p.viewObj.threshAxes, 'CurrentPoint');
            newThresh = max(currentPoint(1,1), 1);
            newThresh = min(newThresh, p.viewObj.threshAxes.XLim(2));
            %set(p.zoomThreshEndH, 'XData', [newThresh, newThresh])
            p.zoomThreshRectH.XData(3:4) = [newThresh; newThresh];
        end
        
        function stopzoomThresh(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            set(p.viewObj.threshAxes, 'XLim', [p.zoomThreshRectH.XData(1) p.zoomThreshRectH.XData(3)]);
            delete(p.zoomThreshRectH)
            p.threshZoom = false;
            %Should notify threshold change
            %newThresh = str2double(p.viewObj.threshValue.String);
            %p.spotTable.setThreshold(p.spotTable.spotChannels{p.channelIdx}, newThresh);
            %p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
            %p.updateMainAxes();
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
      
        function threshValueChange(p, ~, ~)
            if ~isnan(str2double(p.viewObj.threshValue.String)) && isreal(str2double(p.viewObj.threshValue.String))
                newThresh = round(str2double(p.viewObj.threshValue.String));
                set(p.viewObj.thresholdLineH, 'XData', [newThresh, newThresh])
                %Update threshold for spots
                p.spotTable.setThreshold(p.spotTable.spotChannels{p.channelIdx}, newThresh);
                p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
%                 p.updateColorMap(); %Specifying colors directly in scatter instead
                p.updateMainAxes();
            else
                %If not number, return to previous threshold or default
                %threshold
                if isvalid(p.viewObj.thresholdLineH)
                    oldThresh = get(p.viewObj.thresholdLineH, 'XData');
                    set(p.viewObj.threshValue, 'String', num2str(oldThresh(1)))
                else
                    threshold = p.spotTable.thresholds{p.channelIdx};
                    yaxis = get(p.viewObj.threshAxes, 'Ylim');
                    p.viewObj.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis,...
                        'Color', 'b', 'HitTest', 'off');
                    p.viewObj.threshValue.String = num2str(threshold);
                    p.spotTable.setThreshold(p.spotTable.spotChannels{p.channelIdx}, num2str(threshold));
                    p.updateMainAxes();
                end
            end
            
        end
%         function threshValueKey(p, ~, evt) %Could use this to toggle up and down with arrow keys
%             if strcmp(evt.Key, 'return')
%                 disp(p.viewObj.threshValue.String)
%             end
%         end
       
        function p = updateImageInView(p)
            p.imagesInView = cell(0, numel(p.spotTable.spotChannels));
            if p.viewRect(3) * p.viewRect(4) < 4000001
                for i = 1:numel(p.spotTable.spotChannels)
                    p.imagesInView{i} = p.scanObj.getImageRect(p.spotTable.spotChannels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getDapiImage(p.viewRect);
            elseif p.viewRect(3) * p.viewRect(4) < 64000001
                for i = 1:numel(p.spotTable.spotChannels)
                    p.imagesInView{i} = p.scanObj.getSmallImageRect(p.spotTable.spotChannels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getSmallDapiImage(p.viewRect);
            else
                disp('View is too large to display image. Plotting scatter')
                set(p.viewObj.scatterCheckBox, 'Value', 1)
            end
        end
        
        function showImage(p)
           %Adjust contrast
            contrastIn = [get(p.viewObj.lowerContrastSlider, 'Value'), get(p.viewObj.upperContrastSlider, 'Value')];
            tmpIm = imadjust(p.imagesInView{p.channelIdx}, contrastIn, []);
            %Decide if overlay with DAPI?
            tmpRGB = cat(3, tmpIm, tmpIm, tmpIm+p.dapiInView);
            xlimits = [p.viewRect(2),  p.viewRect(2)+p.viewRect(4)];
            ylimits = [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)];
            %axis(p.viewObj.mainAxes, [xlimits, ylimits], 'square')
            set(p.viewObj.mainAxes, 'XLim', xlimits)
            set(p.viewObj.mainAxes, 'YLim', ylimits)
            hold(p.viewObj.mainAxes, 'on')
            p.imageH = imshow(tmpRGB, 'XData', xlimits, 'YData', ylimits, 'Parent', p.viewObj.mainAxes);
            %axis fill
            %pbaspect auto
            set(p.viewObj.mainAxes, 'Visible', 'on')
            hold(p.viewObj.mainAxes, 'off')
        end
        
        function overlaySpots(p, ~, ~)
            if logical(p.viewObj.spotsCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                [spotsInView, spotIdx] = p.spotTable.getValidSpotsInRect(p.spotTable.spotChannels{p.channelIdx}, p.viewRect);
                hold(p.viewObj.mainAxes, 'on')
                p.spotScatterH = scatter(spotsInView.y, spotsInView.x, 20, spotsInView.colors,...
                    'Parent', p.viewObj.mainAxes, 'HitTest','off');
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.spotScatterH)
            end
        end
        
        function overlayNuclei(p, ~, ~)
            if logical(p.viewObj.centroidsCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                outTableTmp = p.spotTable.centroidTableInRect(p.channelIdx, p.viewRect);
                hold(p.viewObj.mainAxes, 'on')
                p.nucleiScatterH = scatter(outTableTmp.y, outTableTmp.x, 35, outTableTmp.colors, 'filled',...
                    'Parent', p.viewObj.mainAxes, 'HitTest','off');
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.nucleiScatterH)
            end
        end
        
        function overlayMasks(p, ~, ~)
            masks = findobj(p.viewObj.mainAxes, 'type', 'images.roi.freehand');
            delete(masks);
            if logical(p.viewObj.masksCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                maskTableTmp = p.maskObj.getChannelMasksInRect(p.viewRect, p.spotTable.spotChannels{p.channelIdx});
                maskIDs = unique(maskTableTmp.maskID);
                maskIDs(maskIDs == 0) = [];
                for i = 1:numel(maskIDs)
                    drawfreehand(p.viewObj.mainAxes, 'Position', maskTableTmp{maskTableTmp.maskID == maskIDs(i), {'y', 'x'}},...
                        'Color', 'red', 'InteractionsAllowed', 'none');
                end
                cellMasksTmp = p.maskObj.getChannelMasksInRect(p.viewRect, 'dapi');
                maskIDs = unique(cellMasksTmp.maskID);
                maskIDs(maskIDs == 0) = [];
                for i = 1:numel(maskIDs)
                    drawfreehand(p.viewObj.mainAxes, 'Position', cellMasksTmp{cellMasksTmp.maskID == maskIDs(i), {'y', 'x'}},...
                        'Color', 'blue', 'InteractionsAllowed', 'none');
                end
            end
        end
        
        function scatterCallback(p, ~, ~)
            if p.viewObj.scatterCheckBox.Value == 0
                p.updateImageInView();
            end
            p.updateMainAxes();
        end
        
        function addSpotMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            channel = p.spotTable.spotChannels{p.channelIdx};
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            %addlistener(p.maskH, 'DrawingFinished', @p.maskSpots);
            if ~isempty(p.maskH.Position) %Allows 'escape' from ROI
                p.maskSpots(p.maskH, channel)
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskSpots(p, roi, channel)
            %tmpPoly = roi.Position;
            set(p.viewObj.masksCheckBox, 'Value', true)
            p.maskObj.addMaskLocalCoords(roi.Position, channel);
            delete(roi)
            tic
            p.spotTable.addNewMask(channel);
            p.spotTable.updateSpotStatus(channel);
            toc
            p.spotTable.updateCentroidList(channel, p.paletteIdx);
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function addCellMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.maskH.Position) %Allows 'escape' from ROI
                p.maskCells(p.maskH)
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskCells(p, roi)
            %tmpPoly = roi.Position;
            channel = p.spotTable.spotChannels{p.channelIdx};
            set(p.viewObj.masksCheckBox, 'Value', true)
            p.maskObj.addMaskLocalCoords(roi.Position, 'dapi');
            delete(roi)
            p.nucleiObj.addNewMask();
            p.spotTable.updateSpotStatus(channel);
            p.spotTable.updateCentroidList(channel, p.paletteIdx);
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function deleteMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            channel = p.spotTable.spotChannels{p.channelIdx};
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.   
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.maskObj.removeMasksByLocalPoints(ptsInView, p.viewRect);
                tic
                p.nucleiObj.removeMasks();
                toc
                %p.nucleiObj.updateMasksInRect(p.viewRect);
                tic
                p.spotTable.removeMasks2(channel, p.viewRect);
                toc
                %p.spotTable.updateMasksInRect(channel, p.viewRect);
                p.spotTable.updateCentroidList(channel, p.paletteIdx);
                p.updateCentroidListView();
                p.updateMainAxes();
                set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            end
        end
        
        function addCells(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.   
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.nucleiObj.addCell(ptsInView(:,2), ptsInView(:,1));
                tic
                p.spotTable.assignSpotsInRect(p.viewRect)
                toc
                tic
                p.spotTable.updateAllSpotStatus();
                toc
                tic
                p.spotTable.makeCentroidList(p.paletteIdx);
                toc
                p.updateCentroidListView();
                p.updateMainAxes();
                set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            end
        end
        
        function deleteCells(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.nucleiObj.removeCell(ptsInView(:,2), ptsInView(:,1));
                tic
                p.spotTable.assignSpotsInRect(p.viewRect);
                toc
                tic
                p.spotTable.updateAllSpotStatus();
                toc
                tic
                p.spotTable.makeCentroidList(p.paletteIdx);
                toc
                p.updateCentroidListView();
                p.updateMainAxes();
                set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            end
        end
        
        function type = getSelectionType(p)
            type = get(p.viewObj.figHandle, 'SelectionType');
        end
                
        function keyPressFcns(p, ~, evt)
            keyPressed = evt.Key;
            modifierPressed = evt.Modifier;
            if isempty(modifierPressed)
                switch(keyPressed)
                    case 'z'
                        p.zoomInPressed();
                    case 'p'
                        p.panViewPressed();
                    case 's'
                        p.overlaySpots();
                    case 'n'
                        p.overlayNuclei();
                    case 'c'
                        p.scatterCallback();
                    case 'm'
                        disp('m')
                        p.addSpotMask();
                    case 'd'
                        p.deleteMask();
                    case 'f' 
                        p.fixedZoom = true;
                    case 'x'
                        set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
                        %iptPointerManager(p.viewObj.figHandle, 'disable');
                end
            elseif strcmp(modifierPressed{1}, 'shift')
                switch(keyPressed)
                    case 'm'
                        disp('shift M')
                        p.addCellMask();
                    case 's'
                    case 'e'
                        
                end
            end
        end
    end
end