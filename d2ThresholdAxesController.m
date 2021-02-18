classdef d2ThresholdAxesController < handle
    
    properties
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        
        %View
        viewObj
        
        %Other controllers
        mainAxesCntrlr
        thumbCntrlr
        
        ThreshAxisMin
        ThreshAxisMax
        threshZoom = false
        histogramLineH
        thresholdLineH
        zoomThreshRectH
        includeMaskedSpots = true
        
    end
    
    methods
        function p = d2ThresholdAxesController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
        end
        
        function startup(p)
            p.spotTable.allIntensities();
            p.spotTable.allIntensitiesNoMasked();
            p.plotIntensityHistogram()
            p.plotIntensityThreshold()
        end
        
        function plotIntensityHistogram(p)
            if ~isempty(p.histogramLineH)
                delete(p.histogramLineH)
            end
            
            if p.includeMaskedSpots
                intensities = p.spotTable.spotsIntensitiesWithMasked{p.mainAxesCntrlr.channelIdx};
            else
                intensities = p.spotTable.spotsIntensitiesNoMasked{p.mainAxesCntrlr.channelIdx};
            end
            
            logRank = log(numel(intensities):-1:1);
            p.histogramLineH = line(p.viewObj.threshAxes, intensities, logRank, ...
                'HitTest', 'off', ...
                'Color', 'k');
            
            p.ThreshAxisMin = intensities(1);
            p.ThreshAxisMax = intensities(end) * 1.05;
            
            set(p.viewObj.threshAxes, 'XLim', [p.ThreshAxisMin p.ThreshAxisMax]);
            
            yaxismax = logRank(1)*1.1;
          
            set(p.viewObj.threshAxes, 'YLim', [0 yaxismax]);
        end
        
        function plotIntensityThreshold(p)
            if ~isempty(p.thresholdLineH) && ishandle(p.thresholdLineH)
                delete(p.thresholdLineH)
            end
            if isempty(p.spotTable.thresholds)
                p.spotTable.defaultThresholds();
            end
            threshold = p.spotTable.thresholds{p.mainAxesCntrlr.channelIdx};
            yaxis = get(p.viewObj.threshAxes, 'Ylim');
            p.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis,...
                'Color', 'b', 'HitTest', 'off');
            p.viewObj.threshValue.String = num2str(threshold);
        end
        
        function threshZoomButtonDown(p, ~, ~)
            p.threshZoom = true;
        end
        
        function thresholdButtonDown(p, ~, ~)
            switch(p.mainAxesCntrlr.getSelectionType)
                case 'normal'
                    currentPoint = get(p.viewObj.threshAxes, 'CurrentPoint');
                    if ~logical(p.threshZoom)
                        if abs(currentPoint(1,1) - p.thresholdLineH.XData(1)) < 150
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
            set(p.thresholdLineH, 'XData', [newThresh, newThresh])
            p.viewObj.threshValue.String = num2str(round(newThresh));
        end
        
        function stopThreshDrag(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            newThresh = str2double(p.viewObj.threshValue.String);
            p.spotTable.setThreshold(p.spotTable.spotChannels{p.mainAxesCntrlr.channelIdx}, newThresh);
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.mainAxesCntrlr.channelIdx}.GroupCount);
            p.mainAxesCntrlr.updateMainAxes();
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
            set(p.viewObj.threshAxes, 'XLim', [min(p.zoomThreshRectH.XData) max(p.zoomThreshRectH.XData)]);
            delete(p.zoomThreshRectH)
            p.threshZoom = false;
            %Should notify threshold change
            %newThresh = str2double(p.viewObj.threshValue.String);
            %p.spotTable.setThreshold(p.spotTable.spotChannels{p.mainAxesCntrlr.channelIdx}, newThresh);
            %p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.mainAxesCntrlr.channelIdx}.GroupCount);
            %p.mainAxesCntrlr.updateMainAxes();
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
      
        function threshValueChange(p, ~, ~)
            if ~isnan(str2double(p.viewObj.threshValue.String)) && isreal(str2double(p.viewObj.threshValue.String))
                newThresh = round(str2double(p.viewObj.threshValue.String));
                set(p.thresholdLineH, 'XData', [newThresh, newThresh])
                %Update threshold for spots
                p.spotTable.setThreshold(p.spotTable.spotChannels{p.mainAxesCntrlr.channelIdx}, newThresh);
                p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.mainAxesCntrlr.channelIdx}.GroupCount);
%                 p.updateColorMap(); %Specifying colors directly in scatter instead
                p.mainAxesCntrlr.updateMainAxes();
            else
                %If not number, return to previous threshold or default
                %threshold
                if isvalid(p.thresholdLineH)
                    oldThresh = get(p.thresholdLineH, 'XData');
                    set(p.threshValue, 'String', num2str(oldThresh(1)))
                else
                    threshold = p.spotTable.thresholds{p.mainAxesCntrlr.channelIdx};
                    yaxis = get(p.viewObj.threshAxes, 'Ylim');
                    p.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis,...
                        'Color', 'b', 'HitTest', 'off');
                    p.viewObj.threshValue.String = num2str(threshold);
                    p.spotTable.setThreshold(p.spotTable.spotChannels{p.mainAxesCntrlr.channelIdx}, num2str(threshold));
                    p.mainAxesCntrlr.updateMainAxes();
                end
            end
            
        end
        
        function filterMasksPushed(p, ~, ~)
            p.includeMaskedSpots = ~p.includeMaskedSpots;
            if p.includeMaskedSpots
                set(p.viewObj.filterMasksThresh, 'String', 'filter masked spots')
            else
                p.spotTable.allIntensitiesNoMasked(); %Update intensities excluding masked spots
                set(p.viewObj.filterMasksThresh, 'String', 'remove filter')
            end
            p.plotIntensityHistogram();
        end
%         function threshValueKey(p, ~, evt) %Could use this to toggle up and down with arrow keys
%             if strcmp(evt.Key, 'return')
%                 disp(p.viewObj.threshValue.String)
%             end
%         end

    end
end