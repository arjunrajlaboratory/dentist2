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
        zoomROI
        zoomStart
        zoomRect
        channelIdx
        scatterH
        
    end
    
    methods
        function p = d2ThresholdController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
            
            p.viewRect = [1, 1, max(p.spotTable.centroidLists{1}{:,{'x', 'y'}})];
            p.startup()
        end
        
        function startup(p)
            p.channelIdx = p.viewObj.channelPopup.Value;
            p.spotTable.makeIntensitiesToPlot();
            p.plotScatterMain(p.channelIdx)
            p.plotIntensityHistogram(p.channelIdx)
            p.plotIntensityThreshold(p.channelIdx)
        end
        
        function p = changeChannel(p, src, evt)
            p.channelIdx = p.viewObj.channelPopup.Value;
            %Update centroid listbox
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
            updateMainAxes(p, p.channelIdx);
            
        end
        
        function p = centroidSelected(p, src, evt)
            if strcmp(src.Parent.SelectionType, 'open') %Respond to double mouse click
                p.channelIdx = p.viewObj.channelPopup.Value;
                cellIdx = p.viewObj.centroidList.Value;
                disp(p.spotTable.centroidLists{p.channelIdx}(cellIdx, :));
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
            centroidsInView = p.spotTable.centroidTableInRect(channelIdx, p.viewRect);
            p.scatterH = scatter(centroidsInView.y, centroidsInView.x,...
                20, centroidsInView.GroupCount, 'filled',...
                'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'scatterH');
            set(p.viewObj.mainAxes, 'Xlim', [p.viewRect(2)-1,  p.viewRect(2)+p.viewRect(4)+1])
            set(p.viewObj.mainAxes, 'Ylim', [p.viewRect(1)-1,  p.viewRect(1)+p.viewRect(3)+1])
            set(p.viewObj.mainAxes, 'Ydir','reverse');
            colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
            set(p.viewObj.mainAxes, 'Interactions',[]);
            set(p.viewObj.mainAxes.Toolbar, 'Visible','off');

        end
        
        function figWindowDown(p, ~, ~)
            %disp(get(p.viewObj.mainAxes, 'CurrentPoint'));
            if ~(p.selectionOutSideMainView(get(p.viewObj.mainAxes, 'CurrentPoint')))
                disp('mainAxesClicked')
                %k = p.getSelectionType
                switch(p.getSelectionType)
                    case 'normal'
                        p.zoomStart = get(p.viewObj.mainAxes, 'CurrentPoint');
                        set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopDragFcn});
                        set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.mainAxesZoom});
                    case 'open'
                        p.viewRect = [1, 1, max(p.spotTable.centroidLists{1}{:,{'x', 'y'}})];
                        p.plotScatterMain(p.viewObj.channelPopup.Value)
                        %(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
                end
            end
        end
        
        function mainAxesZoom(p, ~, ~)
            if ishandle(p.viewObj.zoomH)
                delete(p.viewObj.zoomH)
                %delete(p.zoomRect)
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
        
        function pos = selectionOutSideMainView(p, point)
            pos = p.viewRect(1) > point(1,2) || p.viewRect(1)+ p.viewRect(3) < point(1,2)...
                || p.viewRect(2) > point(1,1) || p.viewRect(2)+ p.viewRect(4) < point(1,1)
        end
        
        function getSelectedRectangleCoords(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            %xLimits = xlim(p.axH);
            %yLimits = ylim(p.axH);
            
            %currentPoint(1,1) = max(min(currentPoint(1,1), xLimits(2)), xLimits(1));
            %currentPoint(1,2) = max(min(currentPoint(1,2), yLimits(2)), yLimits(1));
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)),...
                abs(p.zoomStart(1,1:2) - currentPoint(1,1:2))]);
        end
        
        function stopDragFcn(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            if ~isempty(p.zoomRect)
                p.viewRect  = p.zoomRect; % Should notify event "viewChange"
                p.plotScatterMain(p.viewObj.channelPopup.Value)
                p.zoomRect = [];
            end
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
        
        function figButtonDown(p, ~, ~)
            disp('fig button down')
        end
        
%         function figWindowDown(p, ~, ~)
%             disp('window button down')
%             disp(get(p.viewObj.mainAxes, 'CurrentPoint'));
%         end
        
        function plotScatterThumbnail(p)
            
        end
        
        function zoomInPressed(p, ~, ~)
            p.zoomROI = drawrectangle(p.viewObj.mainAxes);
            

        end
        
        function type = getSelectionType(p)
            type = get(p.viewObj.figHandle, 'SelectionType');
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