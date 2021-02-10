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
        
    end
    
    methods
        function p = d2ThresholdController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
            
            p.localRect = [1, 1, max(p.spotTable.centroidLists{1}.x), max(p.spotTable.centroidLists{1}.y)];
            plotScatterMain(p)
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
        
        function plotScatterMain(p)
            channelIdx = p.viewObj.channelPopup.Value;
            centroidsInView = p.spotTable.centroidTableInRect(channelIdx, p.localRect);
            scatter(p.viewObj.mainAxes,...
                centroidsInView.y, centroidsInView.x,...
                20, centroidsInView.GroupCount, 'filled')
            set(p.viewObj.mainAxes, 'Xlim', [p.localRect(2)-1,  p.localRect(2)+p.localRect(4)+1])
            set(p.viewObj.mainAxes, 'Ylim', [p.localRect(1)-1,  p.localRect(1)+p.localRect(3)+1])
            set(p.viewObj.mainAxes, 'Ydir', 'reverse')
            set(p.viewObj.mainAxes, 'Interactions',[]);
            p.viewObj.mainAxes.Toolbar.Visible = 'off';

        end
        
        function plotScatterThumbnail(p)
            
        end
        
        function startUpPlots(p)
            channelIdx = p.viewObj.channelPopup.Value;
            scatter(p.viewObj.mainAxes,...
                p.spotTable.centroidLists{channelIdx}.y, p.spotTable.centroidLists{channelIdx}.x,...
                20, p.spotTable.centroidLists{channelIdx}.GroupCount, 'filled')
            xlim(p.viewObj.mainAxes, [0 max(p.spotTable.centroidLists{channelIdx}.y)+5])
            ylim(p.viewObj.mainAxes, [0 max(p.spotTable.centroidLists{channelIdx}.x)+5])
            %p.viewObj.mainAxes.Xlim = [0 max(p.spotTable.centroidLists{channelIdx}.y)+5];
            %p.viewObj.mainAxes.Ylim = [0 max(p.spotTable.centroidLists{channelIdx}.x)+5];
        end
        
        function p = updateMainAxes(p, channelIdx)
%             channelIdx = p.viewObj.channelPopup.Value;
            scatter(p.viewObj.mainAxes,...
                p.spotTable.centroidLists{channelIdx}.y, p.spotTable.centroidLists{channelIdx}.x,...
                20, p.spotTable.centroidLists{channelIdx}.GroupCount, 'filled')
        end
        
        function p = updateThumbnail(p)
            
        end
        
        function p = plotIntensity(p)
            
        end
        
        function p = plotSpotCounts(p)
            
        end
        
        function updateThresholdPlot(p)
            
        end
        
        
    end
    
end