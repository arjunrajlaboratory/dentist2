classdef d2stitchingView < handle
     
    properties (Access = public)
        %Model objects
        scanObj
        
        stichCntrlr
        
        figHandle
        imageAxesTopL
        imageAxesBottomR
        imageAxesBottomL
        
        startPos
        startPosButtons
        snake
        snakeButtons
        direction
        directionButtons
        
        newPositionButton
        rowValue
        colValue
        
        cpSelectRowButton
        cpSelectColButton
        doneButton
    end
     
    methods
        
        function p = d2stitchingView(scanFile, scanDim)
            p.scanObj = scanObject('scanFile', scanFile);
            p.scanObj.scanDim = scanDim;
            
            createComponents(p)
            
%             startupFcn(p)
        end
        
        function createComponents(p)
            p.figHandle = figure('Visible', 'off', 'Position', [100 100 800 800]);
            p.imageAxesBottomL = axes('Parent', p.figHandle, 'XTickLabel', '', 'YTickLabel', '', 'Position', [0.025 0.025 0.475 0.475], 'Interactions',[]);
            p.imageAxesTopL = axes('Parent', p.figHandle,'XTickLabel', '', 'YTickLabel', '', 'Position', [0.025 0.5 0.475 0.475], 'Interactions',[]);
            p.imageAxesBottomR = axes('Parent', p.figHandle,'XTickLabel', '', 'YTickLabel', '', 'Position', [0.5 0.025 0.475 0.475], 'Interactions',[]);
            
            p.startPos = uibuttongroup(p.figHandle, 'Units', 'normalized', 'Position', [0.525 0.85 0.15 0.1]);
            p.startPosButtons(1) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'top left', 'Units', 'normalized', 'Position', [0.01 0.8 0.9 0.2]);
            p.startPosButtons(2) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'top right' ,'Units', 'normalized', 'Position', [0.01 0.55 0.9 0.2]);
            p.startPosButtons(3) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'bottom left' ,'Units', 'normalized', 'Position', [0.01 0.3 0.9 0.2]);
            p.startPosButtons(4) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'bottom right' ,'Units', 'normalized', 'Position', [0.01 0.05 0.9 0.2]);
            
            p.snake = uibuttongroup(p.figHandle, 'Units', 'normalized', 'Position', [0.525 0.76 0.15 0.06]);
            p.snakeButtons(1) = uicontrol(p.snake, 'Style','radiobutton', 'String', 'snake', 'Units', 'normalized', 'Position', [0.01 0.5 0.9 0.4]);
            p.snakeButtons(2) = uicontrol(p.snake, 'Style','radiobutton', 'String', 'no snake' ,'Units', 'normalized', 'Position', [0.01 0.1 0.9 0.4]);
            
            p.direction = uibuttongroup(p.figHandle, 'Units', 'normalized', 'Position', [0.525 0.67 0.15 0.06]);
            p.directionButtons(1) = uicontrol(p.direction, 'Style','radiobutton', 'String', 'horizontal', 'Units', 'normalized', 'Position', [0.01 0.5 0.9 0.4]);
            p.directionButtons(2) = uicontrol(p.direction, 'Style','radiobutton', 'String', 'vertical' ,'Units', 'normalized', 'Position', [0.01 0.1 0.9 0.4]);
            
%             p.newPositionButton = uicontrol('Style', 'pushbutton', 'String', 'new positions', 'Units', 'normalized', 'Position', [0.755 0.45 0.1111 0.0367]);
%             p.rowValue = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.6 0.6 0.0700 0.0333]);
%             p.colValue = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.65 0.65 0.0700 0.0333]);
%             
%             p.cpSelectRowButton = uicontrol('Style', 'pushbutton', 'String', 'select row control points', 'Units', 'normalized', 'Position', [0.755 0.45 0.1111 0.0367]);
%             p.cpSelectColButton = uicontrol('Style', 'pushbutton', 'String', 'select col control points', 'Units', 'normalized', 'Position', [0.755 0.45 0.1111 0.0367]);
%             p.doneButton = uicontrol('Style', 'pushbutton', 'String', 'done', 'Units', 'normalized', 'Position', [0.755 0.45 0.1111 0.0367]);

            
%             function startupFcn(p)
%             end

            p.figHandle.Visible = 'on';
            
        end
    end
    
end