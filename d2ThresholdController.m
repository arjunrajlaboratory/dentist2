classdef d2ThresholdController < handle
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        
        viewObj
        threshZoom
        viewRect %Image coords not axes coords
        cellViewRadius = 750
        zoomROI
        zoomStart
        zoomRect
        channelIdx
        scatterH %Not sure if we'll need these handle. 
        imageH
        
        imagesInView
        dapiInView
        
        resizedImg %Not sure if we want to save these or regenerate
        resizedDapi
        
    end
    
    methods
        function p = d2ThresholdController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
            
            p.viewRect = [1, 1, size(p.scanObj.dapiStitch)];
            p.startup()
        end
        
        function startup(p)
            p.channelIdx = p.viewObj.channelPopup.Value;
            p.spotTable.makeIntensitiesToPlot();
            p.plotIntensityHistogram()
            p.plotIntensityThreshold()
            p.plotScatterMain();
            p.updateImageInView();
            %p.updateMainAxes()
        end
        
        function p = changeChannel(p, ~, ~)
            p.channelIdx = p.viewObj.channelPopup.Value;
            %Update centroid listbox
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
            p.updateMainAxes();
            p.plotIntensityHistogram();
        end
        
        function p = updateMainAxes(p)
%             disp(p.viewObj.scatterCheckBox.Value)
%             disp(class(p.viewObj.scatterCheckBox.Value))
%             disp(isvalid(p.scatterH))
%             disp(isgraphics(p.scatterH))
%             disp(logical(p.viewObj.scatterCheckBox.Value))
            if isvalid(p.viewObj.mainAxes.Children)
                delete(get(p.viewObj.mainAxes, 'Children'));
            end
            
            if logical(p.viewObj.scatterCheckBox.Value)
                p.plotScatterMain();
            else
                p.showImage();
            end
        end
        
        
        function p = centroidSelected(p, src, ~)
            if strcmp(src.Parent.SelectionType, 'open') %Respond to double mouse click
                p.channelIdx = p.viewObj.channelPopup.Value;
                cellIdx = p.viewObj.centroidList.Value;
                disp(p.spotTable.centroidLists{p.channelIdx}(cellIdx, :));
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
            
            xAxisMin = intensities(1);
            xAxisMax = intensities(end) * 1.05;
            
            set(p.viewObj.threshAxes, 'XLim', [xAxisMin xAxisMax]);
            
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
            p.viewObj.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis, ...
                'Color', 'b', 'HitTest', 'off');
            p.viewObj.threshValue.String = num2str(threshold);
        end
        
        function plotScatterMain(p)
            %Could probably save on time by adding colormap to
            %centroidTable in the spotTableObject.
            centroidsInView = p.spotTable.centroidTableInRect(p.channelIdx, p.viewRect);
            set(p.viewObj.mainAxes, 'Xlim', [p.viewRect(2),  p.viewRect(2)+p.viewRect(4)])
            set(p.viewObj.mainAxes, 'Ylim', [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)])
            hold(p.viewObj.mainAxes, 'on')
            p.scatterH = scatter(centroidsInView.y, centroidsInView.x,...
                20, centroidsInView.GroupCount, 'filled',...
                'Parent', p.viewObj.mainAxes, 'HitTest','off');
            hold(p.viewObj.mainAxes, 'off')
            %colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
        end
        
        function figWindowDown(p, ~, ~)
            if d2utils.pointInSideRect(p.viewRect, get(p.viewObj.mainAxes, 'CurrentPoint'))
                switch(p.getSelectionType)
                    case 'normal'
                        p.zoomStart = get(p.viewObj.mainAxes, 'CurrentPoint');
                        set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopDragFcn});
                        set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.mainAxesZoom});
                    case 'open'
                        p.viewRect = [1, 1, size(p.scanObj.dapiStitch)];
                        p.updateImageInView();
                        p.updateMainAxes();
                end
            end
        end
        
        function mainAxesZoom(p, ~, ~)
            if ishandle(p.viewObj.zoomH)
                delete(p.viewObj.zoomH)
            end
            p.getSelectedRectangleCoords()
            if p.zoomRect(3) > 4 && p.zoomRect(4) > 4
                p.viewObj.zoomH = rectangle(...
                    'Position', p.zoomRect, ...
                    'EdgeColor', 'r', ...
                    'Parent', p.viewObj.mainAxes,...
                    'Hittest', 'off');
             end
        end
        
        function getSelectedRectangleCoords(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)),...
                abs(p.zoomStart(1,1:2) - currentPoint(1,1:2))]);
        end
        
        function stopDragFcn(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            if ~isempty(p.zoomRect) && all(p.zoomRect > 0)
                p.viewRect  = d2utils.coordToPixelRect(p.zoomRect); % Should notify event "viewChange"
                %Update imagesInView
                p.updateImageInView();
                %delete(p.scatterH)
                %Update main axes
                p.updateMainAxes();
                p.zoomRect = [];
            end
            delete(p.viewObj.zoomH)
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
       
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
            set(p.viewObj.mainAxes, 'Xlim', xlimits)
            set(p.viewObj.mainAxes, 'Ylim', ylimits)
            hold(p.viewObj.mainAxes, 'on')
            p.imageH = imshow(tmpRGB, 'XData', xlimits, 'YData', ylimits, 'Parent', p.viewObj.mainAxes);
            %axis fill
            %pbaspect auto
            set(p.viewObj.mainAxes, 'Visible', 'on')
            hold(p.viewObj.mainAxes, 'off')
            
            
        end
        
%         function resizeImageInView(p)
%            %Decide if overlay with DAPI?
%            %Resize if too large
%            if p.viewRect(3) * p.viewRect(4) < 2250001
%                p.resizedImg = p.imagesInView{p.channelIdx};
%                p.resizedDapi = p.dapiInView;
%            elseif p.viewRect(3) * p.viewRect(4) < 6250001
%                p.resizedImg = im2uint8(p.imagesInView{p.channelIdx});
%                p.resizedDapi = im2uint8(p.dapiInView);
%            elseif p.viewRect(3) * p.viewRect(4) < 1000000001
%                p.resizedImg = im2uint8(imresize(p.imagesInView{p.channelIdx}, 1/4));
%                p.resizedDapi = im2uint8(imresize(p.dapiInView, 1/4));
%            else
%                disp('View is too large to display image. Plotting scatter')
%                set(p.viewObj.scatterCheckBox, 'Value', 1)
%                p.updateMainAxes();
%            end
%         end
        
        function scatterCallback(p, ~, ~)
            if p.viewObj.scatterCheckBox.Value == 0
                p.updateImageInView();
            end
            p.updateMainAxes();
        end
        
        function zoomInPressed(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '');
            p.zoomROI = drawrectangle(p.viewObj.mainAxes); %Can add callback for clicking on ROI
        end
        
        function type = getSelectionType(p)
            type = get(p.viewObj.figHandle, 'SelectionType');
        end
        
        function keyPressFcns(p, ~, evt)
            keyPressed = evt.Key;
            switch(keyPressed)
                case 'return'
                    if isvalid(p.zoomROI)
                        p.viewRect = d2utils.coordToPixelRect(p.zoomROI.Position);
                        delete(p.zoomROI)
                        delete(p.scatterH)
                        set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown});
                        p.updateImageInView
                        %Update main axes
                        p.updateMainAxes();
                        %drawrectangle('Position', pos, 'Color', 'r', 'InteractionsAllowed', 'none');
                    end
                    
                
                    
            end
                
        end
        
        
        
        
        

        
        
    end
    
end