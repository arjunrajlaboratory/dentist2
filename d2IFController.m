classdef d2IFController < handle
    properties (Access = public)
        
        %Model objects
        scanObj
        maskObj
        IFboundaries
        IFtable
        
        %View
        viewObj
        
        thumbCntrlr
        
        viewRect %Image coords not axes coords
        panRect
        cellViewRadius = 500
        zoomROI
        zoomStart
        zoomRect
        quantMetricDict
        quantMetric
        fixedZoom = false
        panStart
        channelIdx
        zoomMode = true
        scatterH %Not sure if we need these handles. 
        imageH
        nucleiPlotH
        cellPlotH
        nucH
        cellH
        maskH
        imagesInView
        dapiInView
        
        numCells
    end
    
    methods
        function p = d2IFController(view, scanObj, maskObj, IFboundaries, IFtable)
            p.viewObj = view;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.IFboundaries = IFboundaries;
            p.IFtable = IFtable;
            p.viewRect = [1, 1, min(p.scanObj.stitchDim, [25000 25000])];
            p.startup()
        end
        
        function startup(p)
            p.channelIdx = p.viewObj.channelPopup.Value;
            p.plotScatterMain();
            p.quantMetricDict = containers.Map({'mean nucleus', 'mean cytoplasm', 'sum nucleus', 'sum cytoplasm'},...
                            {'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}); %Hard-coded for now. Can make this more functional in the future
            p.quantMetric = get(p.viewObj.quantMetric.SelectedObject, 'String');
            p.updateImageInView();
            p.numCells = height(p.IFtable.centroidLists{p.channelIdx});
        end
        
        function p = changeChannel(p, ~, ~)
            p.channelIdx = p.viewObj.channelPopup.Value;
            set(p.viewObj.IFCheckBox, 'String', p.IFtable.channels{p.channelIdx});
            %Update centroid listbox
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function changeColormap(p, ~, ~)
            p.IFtable.paletteIdx = p.viewObj.colormapPopup.Value;
%             p.IFtable.updateExpressionColors(); %Add var (e.g. meanNuc,
%             meanCyto etc.)
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
                p.overlayNuclei();
                p.overlayCells();
%                 p.overlayMasks();
            end
        end
        
        function updateCentroidListView(p)
            set(p.viewObj.centroidList, 'String', string(p.IFtable.centroidLists{p.channelIdx}{:,4})); %  %Add var (e.g. meanNuc, meanCyto etc.)
            p.numCells = height(p.IFtable.centroidLists{p.channelIdx}); %Although this value doesn't depend on the channel. Easy to put this here rather than everywhere where cell # can change. 
        end
        
        function p = centroidSelected(p, ~, ~)
            if strcmp(get(p.viewObj.figHandle, 'SelectionType'), 'open') %Respond to double mouse click
                cellIdx = get(p.viewObj.centroidList, 'Value');
                cellPos = p.IFtable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                set(p.viewObj.scatterCheckBox, 'Value', 0)
                p.updateImageInView();
                p.updateMainAxes();
                p.thumbCntrlr.overlayThumbnailRect();
            end
        end

        function quantMetricChanged(p, ~, evt)
            p.quantMetric = evt.NewValue.String;
            cellIdx = get(p.viewObj.centroidList, 'Value');
            cellID = p.IFtable.centroidLists{p.channelIdx}.cellID(cellIdx);
            p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            newCellIdx = find(p.IFtable.centroidLists{p.channelIdx}.cellID == cellID);
            set(p.viewObj.centroidList, 'Value', newCellIdx);
            p.updateMainAxes();
        end
        
        function plotScatterMain(p)
            centroidsInView = p.IFtable.centroidTableInRect(p.channelIdx, p.viewRect);
            centroidsInView = flipud(centroidsInView); %In order to plot cells with 0 expression in the back
            set(p.viewObj.mainAxes, 'XLim', [p.viewRect(2),  p.viewRect(2)+p.viewRect(4)])
            set(p.viewObj.mainAxes, 'YLim', [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)])
            hold(p.viewObj.mainAxes, 'on')
            %Note that it is slightly slower to specify RGB colors with expression_color
            %instead of changing the colormap on the axes however, specifying RGB colors
            %ensures that the colors don't change with the view 
            p.scatterH = scatter(centroidsInView.y, centroidsInView.x,...
                30, centroidsInView.expression_color, 'filled',...
                'Parent', p.viewObj.mainAxes, 'HitTest','off');
            hold(p.viewObj.mainAxes, 'off')
            %colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
        end
        
        function zoomInPressed(p, ~, ~)
            p.zoomMode = true;
            %iptPointerManager(p.viewObj.figHandle, 'disable');
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
            %Below are old functions that draw rectangle using matlab
            %ROI instead of mainAxesZoom fcn. Possibly more efficient but doesn't allow using
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
                        p.viewRect = [1, 1, min(p.scanObj.stitchDim, [25000 25000])];
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                    case 'alt'
                        p.viewRect = d2utils.expandView2x(p.viewRect, p.scanObj.stitchDim);
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
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
                p.thumbCntrlr.overlayThumbnailRect();
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
            if any(p.viewRect(3:4) < p.scanObj.stitchDim) %No need to pan if view is of entire scan
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
            p.viewRect = d2utils.updateViewPanning(p.viewRect, displacement, p.scanObj.stitchDim);
            p.updateImageInView;
            p.updateMainAxes;
            p.thumbCntrlr.overlayThumbnailRect();
        end
        
        function p = updateImageInView(p)
            p.imagesInView = cell(0, numel(p.IFtable.channels));
            if prod(p.viewRect(3:4)) < 4000001
                for i = 1:numel(p.IFtable.channels)
                    p.imagesInView{i} = p.scanObj.getImageRect(p.IFtable.channels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getDapiImage(p.viewRect);
            elseif prod(p.viewRect(3:4)) < 64000001
                for i = 1:numel(p.IFtable.channels)
                    p.imagesInView{i} = p.scanObj.getSmallImageRect(p.IFtable.channels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getSmallDapiImage(p.viewRect);
            else
                disp('View is too large to display image. Plotting scatter')
                set(p.viewObj.scatterCheckBox, 'Value', 1)
            end
        end
        
        function showImage(p)
            if get(p.viewObj.IFCheckBox, 'Value') 
                contrastIn = [get(p.viewObj.lowerContrastSlider, 'Value'), get(p.viewObj.upperContrastSlider, 'Value')];
                tmpIm = imadjust(p.imagesInView{p.channelIdx}, contrastIn, []);
%                 tmpRGB = cat(3, tmpIm, tmpIm, tmpIm+p.dapiInView);
            else
                tmpIm = zeros(size(p.imagesInView{p.channelIdx}), 'uint16');
            end
            
            if get(p.viewObj.dapiCheckBox, 'Value')
                tmpRGB = cat(3, tmpIm, tmpIm, tmpIm+p.dapiInView);
            else
                tmpRGB = repelem(tmpIm, 1, 1, 3);
            end
            
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
        
        function overlayCells(p, ~, ~)
            if logical(p.viewObj.cellBordersCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                cellsInView = p.IFboundaries.getCellBoundariesInRect(p.IFboundaries.channels{p.channelIdx}, p.viewRect);
                [cellXmat, cellYmat] = d2utils.polyvect2patchMats(cellsInView.cellBoundary); 
                hold(p.viewObj.mainAxes, 'on')
                p.cellPlotH = patch(cellYmat, cellXmat, reshape(cellsInView.colors, [], 1, 3),...
                    'FaceAlpha', 0.3, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'NucBoundary');
%                 for i = 1:height(cellsInView)
%                     patch(cellsInView.cellBoundary(i).Vertices(:,2), cellsInView.cellBoundary(i).Vertices(:,1), cellsInView.colors(i,:),...
%                     'FaceAlpha', 0.3, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'CellBoundary');
%                 end
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.cellPlotH)
%                 delete(findobj(p.viewObj.mainAxes, 'Tag', 'CellBoundary')) 
            end
        end
        
        function overlayNuclei(p, ~, ~)
            if logical(p.viewObj.nucleiBordersCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                nucleiInView = p.IFboundaries.getNucBoundariesInRect(p.IFboundaries.channels{p.channelIdx}, p.viewRect);
                [nucXmat, nucYmat] = d2utils.polyvect2patchMats(nucleiInView.nucBoundary);
                hold(p.viewObj.mainAxes, 'on')
                p.nucleiPlotH = patch(nucYmat, nucXmat, reshape(nucleiInView.colors, [], 1, 3),...
                    'FaceAlpha', 0.3, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'NucBoundary');
%                 for i = 1:height(nucleiInView)
%                     patch(nucleiInView.nucBoundary(i).Vertices(:,2), nucleiInView.nucBoundary(i).Vertices(:,1), nucleiInView.colors(i,:),...
%                     'FaceAlpha', 0.3, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'NucBoundary');
%                 end
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.nucleiPlotH)
%                 delete(findobj(p.viewObj.mainAxes, 'Tag', 'NucBoundary'))
            end
        end
%         
%         function overlayMasks(p, ~, ~)
%             masks = findobj(p.viewObj.mainAxes, 'type', 'images.roi.freehand');
%             delete(masks);
%             if logical(p.viewObj.masksCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
%                 maskTableTmp = p.maskObj.getChannelMasksInRect(p.viewRect, p.spotTable.spotChannels{p.channelIdx});
%                 maskIDs = unique(maskTableTmp.maskID);
%                 maskIDs(maskIDs == 0) = [];
%                 for i = 1:numel(maskIDs)
%                     drawfreehand(p.viewObj.mainAxes, 'Position', maskTableTmp{maskTableTmp.maskID == maskIDs(i), {'y', 'x'}},...
%                         'Color', 'red', 'InteractionsAllowed', 'none');
%                 end
%                 cellMasksTmp = p.maskObj.getChannelMasksInRect(p.viewRect, 'dapi');
%                 maskIDs = unique(cellMasksTmp.maskID);
%                 maskIDs(maskIDs == 0) = [];
%                 for i = 1:numel(maskIDs)
%                     drawfreehand(p.viewObj.mainAxes, 'Position', cellMasksTmp{cellMasksTmp.maskID == maskIDs(i), {'y', 'x'}},...
%                         'Color', 'blue', 'InteractionsAllowed', 'none');
%                 end
%             end
%         end
        
        function scatterCallback(p, ~, ~)
            if p.viewObj.scatterCheckBox.Value == 0
                p.updateImageInView();
            end
            p.updateMainAxes();
        end
        
        function shuffleColorsInView(p, ~, ~)
            if ~logical(p.viewObj.scatterCheckBox.Value)
                outTableTmp = p.spotTable.centroidTableInRect(p.channelIdx, p.viewRect);
                randomColors  = single(d2utils.distinguishable_colors(50));
                outTableTmp.colors = randomColors(randi(50, height(outTableTmp), 1),:);
                
                %Update centroid list with new colors
                idx = ismember(p.spotTable.centroidLists{p.channelIdx}.nucID, outTableTmp.nucID);
                p.spotTable.centroidLists{p.channelIdx}.colors(idx, :) = outTableTmp.colors;
                %Update nuclei table wtih new colors
                [idxA, idxB] = ismember(p.nucleiObj.nuclei.nucID, outTableTmp.nucID);
                idxB(idxB == 0) = [];
                p.nucleiObj.nuclei.colors(idxA, :) = outTableTmp.colors(idxB,:);
                %Update spots table with new colors
                [idxA, idxB] = ismember(p.spotTable.spots.nearestNucID, outTableTmp.nucID);
                idxB(idxB == 0) = [];
                p.spotTable.spots.colors(idxA, :) = outTableTmp.colors(idxB,:);
                %Update plots
                p.overlayNuclei;
                p.overlaySpots;
            end
        end
        
        function addNucAndCell(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.deleteCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            p.nucH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.nucH.Position) %Allows 'escape' from ROI
                p.IFtable.addNucleus(fliplr(round(p.nucH.Position)), [0 0], p.IFtable.channels{p.channelIdx}, true)
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            delete(p.nucH)
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.addCellButton, p.viewObj.deleteCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function addEmptyNuc(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.deleteCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            p.nucH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.nucH.Position) %Allows 'escape' from ROI
                p.IFtable.addNucleus(fliplr(round(p.nucH.Position)), [0 0], p.IFtable.channels{p.channelIdx}, false);
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            delete(p.nucH)
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.addCellButton, p.viewObj.deleteCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function addCell(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.deleteCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            p.cellH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.cellH.Position) %Allows 'escape' from ROI
                p.IFtable.addCell(fliplr(round(p.cellH.Position)), [0 0], p.IFtable.channels{p.channelIdx});
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            delete(p.cellH)
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.addCellButton, p.viewObj.deleteCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function deleteNuc(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.addCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
            if ~isempty(x)
                p.IFtable.deleteNuc([x, y], p.IFboundaries.channels{p.channelIdx}, p.viewRect);
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.addCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function deleteCell(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.addCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
            if ~isempty(x)
                p.IFtable.deleteCell([x, y], p.IFboundaries.channels{p.channelIdx}, p.viewRect);
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
%             set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.maskCellButton, p.viewObj.addCellButton,...
%                 p.viewObj.maskSpotsButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function type = getSelectionType(p)
            type = get(p.viewObj.figHandle, 'SelectionType');
        end
        
        function keyPressFunctions(p, ~, evt)
            keyPressed = evt.Key;
            modifierPressed = evt.Modifier;
            if isempty(modifierPressed)
                switch(keyPressed)
                    case 'z'
                        p.zoomInPressed();
                    case 'p'
                        p.panViewPressed();
%                     case 's'
%                         p.overlaySpots();
                    case 'n'
                        p.overlayNuclei();
                    case 's'
                        p.scatterCallback();
                    case 'm'
%                         p.addSpotMask();
                    case 'd'
                        p.deleteMask();
                    case 'f' 
                        p.fixedZoom = true;
                    case 'x'
                        set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
                    case 'uparrow'
                        cellIdx = max(1, get(p.viewObj.centroidList, 'Value')-1);
                        cellPos = p.spotTable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                        p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                        set(p.viewObj.scatterCheckBox, 'Value', 0)
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                    case 'downarrow'
                        cellIdx = min(get(p.viewObj.centroidList, 'Value')+1, p.numCells);
                        cellPos = p.spotTable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                        p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                        set(p.viewObj.scatterCheckBox, 'Value', 0)
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                end
            elseif strcmp(modifierPressed{1}, 'shift')
                switch(keyPressed)
                    case 'm'
%                         p.addCellMask();
                    case 's'
                        fprintf('Saving mask table, cell table, and spot tables.\nThis may take a minute\n')
                        p.nucleiObj.saveNucleiTable;
                        p.spotTable.updateScanSummary;
                        p.spotTable.saveSpotsTable;
                        p.maskObj.saveMasksTable;
                    case 'e'
                        p.spotTable.exportSpotsSummary;
                end
            end
        end
    end
end