classdef d2ThresholdView2 < handle
    
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        
        mainAxesCntrlr
        threshAxesCntrlr
        thumbAxesCntrlr
        
        figHandle
        mainAxes 
        thumbAxes
        threshAxes
        centroidList %List box 
        spotsCheckBox %Check box
        centroidsCheckBox %Check box
        scatterCheckBox %Check box
        masksCheckBox  %Check box
        
        channelPopup %Drop down
        colormapPopup %Drop down
        saveButton %Push button
        exportButton %Push button
        zoomAxes %Push button
        panAxes %Push button
        
        zoomThresh %Push button
        filterMasksThresh %Push button
        threshValue %Edit field
        
        addCellButton %Toggle button
        deleteCellButton %Toggle button
        maskCellButton %Push button
        maskSpotsButton %Push button
        deleteMaskButton %Push button
        
        upperContrastSlider %Slider
        lowerContrastSlider %Slider
        sliderLabel
        
        zoomH
        zoomRect
        histogramLineH
        thresholdLineH
        
        showSpots = true %Not sure if these handles are necessary. 
        showCentroids = true
        showScatter= true
        showMasks = true
        
    end
    
    % GUI startup and deletion
    methods (Access = public)

        % Construct app
        function p = d2ThresholdView2(scanObject, spotTable, maskObject, nucleiObject)
            p.scanObj = scanObject;
            p.spotTable = spotTable;
            p.maskObj = maskObject;
            p.nucleiObj = nucleiObject;

            % Create UIFigure and components
            createComponents(p)
            
            %registerApp(p, p.figHandle)
            startupFcn(p)
            %runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

%         % Code that executes before app deletion
%         function delete(p)
%             %save tables
% 
%             % Delete UIFigure when app is deleted
%             delete(app.UIFigure)
%         end
    end
    
    % GUI construction
    methods 
        
        function createComponents(p)

            p.figHandle = figure('Visible', 'off', 'Position', [100 100 1550 900]);
            %iptPointerManager(p.figHandle);
            p.mainAxes = axes('Parent', p.figHandle, 'Position', [0.025 0.025 0.60 0.95], 'Ydir','reverse', 'XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Interactions',[]);
            set(p.mainAxes.Toolbar, 'Visible','off');
            p.thumbAxes = axes('Parent', p.figHandle, 'XTickLabel', '', 'YTickLabel', '', 'Ydir', 'reverse', 'Position', [0.64 0.645 0.21 0.33], 'Interactions', []);
            set(p.thumbAxes.Toolbar, 'Visible', 'off')
            set(p.thumbAxes, 'Colormap', flipud(hot))
            p.threshAxes = axes('Parent', p.figHandle, 'Position', [0.64 0.025 0.34 0.34], 'Interactions', []);
            set(p.threshAxes.Toolbar, 'Visible', 'off')
            
            p.channelPopup = uicontrol('Style', 'popupmenu', 'String', p.spotTable.spotChannels, 'Units', 'normalized', 'Position', [0.64 0.5067 0.1111 0.0367]);
            p.colormapPopup = uicontrol('Style', 'popupmenu', 'String', p.spotTable.expressionColorPal, 'Units', 'normalized', 'Position', [0.755 0.5067 0.1111 0.0367]);
            p.centroidList = uicontrol('Style', 'listbox', 'String', string(p.spotTable.centroidLists{1}.GroupCount),'Units', 'normalized', 'Position', [0.86 0.645 0.12 0.33]);
            p.spotsCheckBox = uicontrol('Style', 'checkbox', 'String', 'spots (s)', 'Value', p.showSpots,'Units', 'normalized', 'Position', [0.65 0.585 0.0722 0.0333]);
            p.centroidsCheckBox = uicontrol('Style', 'checkbox', 'String', 'nuclei (n)', 'Value', p.showCentroids, 'Units', 'normalized', 'Position', [0.695 0.585 0.0722 0.0333]);
            p.masksCheckBox = uicontrol('Style', 'checkbox', 'String', 'masks', 'Value', p.showMasks, 'Units', 'normalized', 'Position', [0.74 0.585 0.0722 0.0333]);
            p.scatterCheckBox = uicontrol('Style', 'checkbox', 'String', 'scatter (c)', 'Value', p.showScatter, 'Units', 'normalized', 'Position', [0.785 0.585 0.0722 0.0333]);
            
            p.addCellButton = uicontrol('Style', 'pushbutton', 'String', 'add cells', 'Units', 'normalized', 'Position', [0.755 0.45 0.1111 0.0367]);
            p.deleteCellButton = uicontrol('Style', 'pushbutton', 'String', 'delete cells', 'Units', 'normalized', 'Position', [0.755 0.41 0.1111 0.0367]);

            p.maskSpotsButton = uicontrol('Style', 'pushbutton', 'String', 'mask spots (m)', 'Units', 'normalized', 'Position', [0.64 0.37 0.1111 0.0367]);
            p.maskCellButton = uicontrol('Style', 'pushbutton', 'String', 'mask cells (M)', 'Units', 'normalized', 'Position', [0.755 0.37 0.1111 0.0367]);
            p.deleteMaskButton = uicontrol('Style', 'pushbutton', 'String', 'delete mask (d)', 'Units', 'normalized', 'Position', [0.87 0.37 0.1111 0.0367]); 

            p.zoomAxes = uicontrol('Style', 'pushbutton', 'String', 'zoom (z)', 'Units', 'normalized', 'Position', [0.64 0.45 0.1111 0.0367]);
            p.panAxes = uicontrol('Style', 'pushbutton', 'String', 'pan view (p)', 'Units', 'normalized', 'Position', [0.64 0.41 0.1111 0.0367]);
            p.saveButton = uicontrol('Style', 'pushbutton', 'String', 'save (S)', 'Units', 'normalized', 'Position', [0.87 0.45 0.1111 0.0367]);
            p.exportButton = uicontrol('Style', 'pushbutton', 'String', 'export (E)', 'Units', 'normalized', 'Position', [0.87 0.41 0.1111 0.0367]);

            p.zoomThresh = uicontrol('Style', 'pushbutton', 'String', 'zoom', 'Units', 'normalized', 'Position', [0.91 0.32 0.0700 0.0333]);
            p.filterMasksThresh = uicontrol('Style', 'pushbutton', 'String', 'filter masked spots', 'Units', 'normalized', 'Position', [0.91 0.28 0.0700 0.0333]);
            p.threshValue = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.835 0.32 0.0700 0.0333]);

            p.upperContrastSlider = uicontrol('Style', 'slider', 'Value', 0.5, 'Units', 'normalized', 'Position', [0.8444 0.61 0.1233 0.0050]);
            p.lowerContrastSlider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.8444 0.58 0.1233 0.0050]);
            p.sliderLabel = uicontrol('Style', 'text', 'String', 'Adjust contrast', 'Units', 'normalized', 'Position', [0.8750 0.525 0.07 0.03]);

            p.figHandle.Visible = 'on';
        end
        
        %attach controller
        function startupFcn(p)
            p.mainAxesCntrlr = d2MainAxesController(p, p.scanObj, p.spotTable, p.maskObj, p.nucleiObj);  
            p.attatchMainAxesController(p.mainAxesCntrlr);        
            
%             p.threshAxesCntrlr = d2ThresholdAxesController(p, p.scanObj, p.spotTable, p.maskObj, p.nucleiObj);
%             p.attachThresholdController(p.threshAxesCntrlr)
%             
%             p.thumbAxesCntrlr = d2ThumbnailAxesController(p, p.scanObj, p.spotTable, p.maskObj, p.nucleiObj);  
%             p.attachThumbnailController(p.thumbAxesCntrlr)
            
            %Connect controllers;  
        end
        
        function p = attatchMainAxesController(p, controller)
            p.channelPopup.Callback = {@controller.changeChannel};
            p.colormapPopup.Callback = {@controller.changeColormap};
            p.centroidList.Callback = {@controller.centroidSelected};
            p.spotsCheckBox.Callback = {@controller.overlaySpots};
            p.centroidsCheckBox.Callback = {@controller.overlayNuclei};
            p.scatterCheckBox.Callback = {@controller.scatterCallback};
            p.masksCheckBox.Callback = {@controller.overlayMasks};
            p.addCellButton.Callback = {@controller.addCells};
            p.deleteCellButton.Callback = {@controller.deleteCells};
            p.maskSpotsButton.Callback = {@controller.addSpotMask};
            p.maskCellButton.Callback = {@controller.addCellMask};
            p.deleteMaskButton.Callback = {@controller.deleteMask};
            %p.zoomInAxes.Callback = {@contrsoller.};
            %p.zoomOutAxes.Callback = {@controller.};
            %p.saveButton.Callback = {@controller.};
            %p.exportButton.Callback = {@controller.};
            p.zoomAxes.Callback = {@controller.zoomInPressed};
            p.panAxes.Callback = {@controller.panViewPressed};
            p.upperContrastSlider.Callback = {@controller.updateMainAxes};
            p.lowerContrastSlider.Callback = {@controller.updateMainAxes};
            %p.saveButton.Callback = {@controller.};
            %p.exportButton.Callback = {@controller.};
            p.threshValue.Callback = {@controller.threshValueChange};
            %p.threshValue.KeyPressFcn = {@controller.threshValueKey};
            %p.mainAxes.ButtonDownFcn = {@controller.mainAxesButtonDown};
            p.threshAxes.ButtonDownFcn = {@controller.thresholdButtonDown};
            p.zoomThresh.Callback = {@controller.threshZoomButtonDown};
            p.thumbAxes.ButtonDownFcn = {@controller.thumbAxesButtonDown};
            p.figHandle.WindowButtonDownFcn = {@controller.figWindowDown};
            p.figHandle.KeyPressFcn = {@controller.keyPressFcns};
            p.figHandle.CloseRequestFcn = {@p.closeFigFcn};
        end
        
        function p = attachThumbnailController(p, controller)
        end
        
        function p = attachThresholdController(p, controller)
        end
        
        function closeFigFcn(p, ~, ~)
            delete(p.figHandle)
        end
        
    end

end
