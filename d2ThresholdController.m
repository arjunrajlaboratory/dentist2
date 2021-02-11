classdef d2ThresholdController < handle
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        
        viewObj
        threshZoom
        localRect %Image coords not axes coords
        zoomROI
        
    end
    
    methods
        function p = d2ThresholdController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
            
            p.localRect = [1, 1, max(p.spotTable.centroidLists{1}.x), max(p.spotTable.centroidLists{1}.y)];
            p.startup()
        end
        
        function startup(p)
            channelIdx = p.viewObj.channelPopup.Value;
            p.spotTable.makeIntensitiesToPlot();
            %p.plotScatterMain(channelIdx)
            p.plotIntensityHistogram(channelIdx)
            p.plotIntensityThreshold(channelIdx)
        end
        
        function p = changeChannel(p, src, evt)
            channelIdx = p.viewObj.channelPopup.Value;
            %Update centroid listbox
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{channelIdx}.GroupCount);
            updateMainAxes(p, channelIdx);
            
        end
        
        function p = centroidSelected(p, src, evt)
            if strcmp(src.Parent.SelectionType, 'open') %Respond to double mouse click
                channelIdx = p.viewObj.channelPopup.Value;
                cellIdx = p.viewObj.centroidList.Value;
                disp(p.spotTable.centroidLists{channelIdx}(cellIdx, :));
            end
        end
        
        function plotIntensityHistogram(p, channelIdx)
            %channelIdx = p.viewObj.channelPopup.Value;
            if ~isempty(p.viewObj.histogramLineH) && ishandle(p.viewObj.histogramLineH)
                delete(p.viewObj.histogramLineH)
            end
            intensities = p.spotTable.intensitiesToPlot{channelIdx};
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
        
        function plotIntensityThreshold(p, channelIdx)
            %channelIdx = p.viewObj.channelPopup.Value;
            if ~isempty(p.viewObj.thresholdLineH) && ishandle(p.viewObj.thresholdLineH)
                delete(p.viewObj.thresholdLineH)
            end
            if isempty(p.spotTable.thresholds)
                p.spotTable.defaultThresholds();
            end
            threshold = p.spotTable.thresholds{channelIdx};
            yaxis = get(p.viewObj.threshAxes, 'Ylim');
            p.viewObj.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis, ...
                'Color', 'b', 'HitTest', 'off');
            p.viewObj.threshValue.String = num2str(threshold);
        end
        
        function plotScatterMain(p, channelIdx)
            %channelIdx = p.viewObj.channelPopup.Value;
            centroidsInView = p.spotTable.centroidTableInRect(channelIdx, p.localRect);
            scatter(p.viewObj.mainAxes, centroidsInView.y, centroidsInView.x,...
                20, centroidsInView.GroupCount, 'filled',...
                'HitTest','off')
            colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
            set(p.viewObj.mainAxes, 'Xlim', [p.localRect(2)-1,  p.localRect(2)+p.localRect(4)+1])
            set(p.viewObj.mainAxes, 'Ylim', [p.localRect(1)-1,  p.localRect(1)+p.localRect(3)+1])
            set(p.viewObj.mainAxes, 'Ydir','reverse');
            set(p.viewObj.mainAxes, 'Interactions',[]);
            set(p.viewObj.mainAxes.Toolbar, 'Visible','off');

        end
        
        function plotScatterThumbnail(p)
            
        end
        
        function zoomInPressed(p, ~, ~)
            p.zoomROI = drawrectangle(p.viewObj.mainAxes);
            

        end
        
        function keyPressFcns(p, src, evt)
            keyPressed = evt.Key;
            switch(keyPressed)
                case 'return'
                    if isvalid(p.zoomROI)
                        pos = p.zoomROI.Position
                        delete(p.zoomROI)
                        drawrectangle('Position', pos, 'Color', 'r', 'InteractionsAllowed', 'none');
                    end
                    
                
                    
            end
                
        end
        
        
        
        
        

        
        
    end
    
end