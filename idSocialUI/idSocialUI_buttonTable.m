function [panelH,propTH, panel2H] = idSocialUI_buttonTable(parentH,position,dataCell,callbackCell,tooltipCell)
fontsize = 12;

if nargin<1 || isempty(parentH)
    parentH = ...
        figure('Units','centimeters','Position',[16 2 19 17],'Name','Options', 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off','WindowStyle','normal');
end

panelH = [];
panel2H = [];
sliderH = [];
% parentH is a cell of handles of buttons. "Update mode"
parentStraight = parentH(:);
if numel(parentH)>1 && iscell(parentH) && ~all(cellfun(@(x) isempty(x),parentStraight)) && ...
        all(cellfun(@(x) ishghandle(x),parentStraight(cellfun(@(x) ~isempty(x),parentStraight))))
    panel2H = get(parentH{1},'Parent');
    panelH = get(get(parentH{1},'Parent'),'Parent');
    dummyH = parentH;
    parentH = get(get(get(parentH{1},'Parent'),'Parent'),'Parent');
    delete([dummyH{:}])
    sliderH = findobj(findobj(get(panelH,'Children')),'Style','slider');
end

if nargin<2 || isempty(position)
    if isempty(panelH)
        position = [.6 .1 .3 .3];
    else
        position = get(panelH,'Position');
    end
end


if nargin<3 || isempty(dataCell)
    
    dataCell = [strcat(repmat({'Option_Name_'},[9,1]),cellstr(num2str((1:9)'))), num2cell(1:9)'];
    dataCell{3,2} = [];
    no_cols = 2;
    
else
    no_cols = size(dataCell,2);
end

if nargin<4 || isempty(callbackCell)
    callbackCell = repmat({@helloworld},[size(dataCell,1),no_cols]);
end

if nargin<5 || isempty(tooltipCell)
    tooltipCell = [];
end


sliderWidthFixedCM = .5;
rowHeightCm =.7;



no_rows = size(dataCell,1);
colWidth = [.7 .3];
if ~sum(colWidth) == 1
    warning([mfilename ': Column width does not sum up to 1. Set to 0.5'])
    colWidth = [.5 .5];
end

last_figpos = -1;



if isempty(panelH)
panelH = uipanel(parentH,'Units','normalized',...
    'Position', position);
end

% Want fixed slider width, but norm. slider height:
set(panelH,'Units','centimeters');
figSizeCM = get(panelH,'Position');
set(panelH,'Units','normalized');

sliderWidthNorm = sliderWidthFixedCM/figSizeCM(3);

if isempty(panel2H)
    panel2H = uipanel(panelH,'Units','normalized',...
        'Position', [0 0 1-sliderWidthNorm 1]);
end

% set(parentH,'Resize','off');%'CloseRequestFcn',@closeOptionWindow);


if isempty(sliderH)
    sliderH = uicontrol('Style','Slider','Parent',panelH,...
        'Units','normalized','Position',[1-sliderWidthNorm 0 sliderWidthNorm 1],...
        'Value',1,'Callback',{@slider_callback1,panel2H});
end

if no_rows*rowHeightCm < figSizeCM(4)
    set(sliderH,'Enable','off')
end



% set(parentH,'WindowScrollWheelFcn', {@mouseScroll,panel2H,sliderH});

% gridy = 1-rowHeightCm:-rowHeightCm:0;
propTH = renderTable(panel2H,dataCell,rowHeightCm);


if verLessThan('matlab', '8.4')
    set(parentH,'ResizeFcn',{@changeSizFunc,panelH,panel2H,sliderH,propTH,dataCell});
else
    set(parentH,'SizeChangedFcn',{@changeSizFunc,panelH,panel2H,sliderH,propTH,dataCell} );
end

% Scale normalized to centimeter
set(panelH,'Units','centimeters');
panel1posCM = get(panelH,'Position');
set(panelH,'Units','normalized');

set(panel2H,'Units','centimeters');
panel2posCM = get(panel2H,'Position');
set(panel2H,'Position',[panel2posCM(1) -(rowHeightCm*no_rows-panel1posCM(4)) panel2posCM(3) rowHeightCm*no_rows]);
panel2posCM = get(panel2H,'Position');


    function slider_callback1(src,eventdata,arg1)
        set(src,'Units','normalized')
        val = get(src,'Value');
        set(arg1,'Units','centimeter')
        act_pos = get(arg1,'Position');
        
        
        set(arg1,'Position', [act_pos(1) -val*abs(panel2posCM(2)) act_pos(3) act_pos(4)])
%         set(src,'Value',abs(act_pos(2))/abs(panel2posCM(2)))
        set(arg1,'Units','normalized')
    end
    function mouseScroll(src,eventdata, handles,arg1)
        set(handles,'Units','centimeter')
        act_pos = get(handles,'Position');
            if eventdata.VerticalScrollCount > 0 && act_pos(2) <0
                val = rowHeightCm;
                    act_pos(2) = min(act_pos(2)+val,0);
            elseif eventdata.VerticalScrollCount < 0 && act_pos(2) > panel2posCM(2);%-(rowHeightCm*no_rows-panel2posCM(4))
                val = -rowHeightCm;
                act_pos(2) = max(act_pos(2)+val,panel2posCM(2));
            end
        
            
            set(handles,'Position',[act_pos(1)  act_pos(2) act_pos(3) act_pos(4)])

      
            set(arg1,'Value',abs(act_pos(2))/abs(panel2posCM(2)))
            set(handles,'Units','normalize')
    end
    
    function changeSizFunc(src,eventdata,handles1, handles, sliderH,propTHInt,rowNames)
        
        % Adjust panels (? WHHHYYYY is this necessary?)
        set(handles1,'Units','normalized');
        set(handles1,'Position', position)
        set(handles,'Units','normalized');
        set(handles,'Position', [0 0 1-sliderWidthNorm 1])
        
        
        % Adjust slider width:
        set(handles1,'Units','centimeters');
        figSizeCM = get(handles1,'Position');
        set(handles1,'Units','normalized');
        sliderWidthNorm = sliderWidthFixedCM/figSizeCM(3);
        set(sliderH,'Position',[1-sliderWidthNorm 0 sliderWidthNorm 1]);
        
        if no_rows*rowHeightCm <= figSizeCM(4)
            set(sliderH,'Enable','off')
        else
            set(sliderH,'Enable','on')
        end
        
        % Adjust table column width
        set(handles,'Units','centimeters');
        
        set(handles,'Units','centimeters');
        panel2posCMIntern = get(handles,'Position');
        set(handles,'Position',[panel2posCMIntern(1) -(rowHeightCm*no_rows-panel2posCMIntern(4)) panel2posCMIntern(3) rowHeightCm*no_rows]);
        panel2posCMIntern = get(handles,'Position');
        colWidthCmInt = panel2posCMIntern(3)*colWidth;
        gridyInt = panel2posCMIntern(4)-rowHeightCm:-rowHeightCm:0;%panel2posCMIntern(2);
        %         gridxInt = panel2posCMIntern(1):colWidthCmInt:panel2posCMIntern(1)+panel2posCMIntern(3);
%         if numel(gridyInt)<60
%             keyboard
%         end
        gridxInt = panel2posCMIntern(1)+[ 0 colWidthCmInt];
        
        if ~isempty(gridyInt) && ~isempty(gridxInt)
            for row = 1:no_rows
                for col = 1:no_cols
                    emptcols = find(cellfun(@(x) isempty(x),rowNames(row,col:end)),1,'first');
                    if ~isempty(emptcols)
                        col_effect = col + emptcols - 1;
                        colWidth_effective = sum(colWidthCmInt(col:col_effect));
                    else
                        colWidth_effective = colWidthCmInt(col);
                    end
                    set(propTHInt{row,col},'Units','centimeters');
                    set(propTHInt{row,col},'Position',[gridxInt(col) gridyInt(row)  colWidth_effective rowHeightCm])
                    set(propTHInt{row,col},'Units','normalized');
                end
            end
        end
        
        set(handles,'Units','normalized');
        
        % Adjust panel2 size
        set(handles1,'Units','centimeters');
        panel1posCM = get(handles1,'Position');
        set(handles1,'Units','normalized');
        set(handles,'Units','centimeters');
        panel2posCM = get(handles,'Position');
        set(handles,'Position',[panel2posCM(1) -(rowHeightCm*no_rows-panel1posCM(4)) panel2posCM(3) rowHeightCm*no_rows]);
        panel2posCM = get(handles,'Position');
        set(handles,'Units','normalized');

     
    end

    function propTHInt = renderTable(panel2HInt,rowNames,rowHeightCm)
        
        no_rowsIntern = size(rowNames,1);
        set(panel2HInt,'Units','centimeters');
        panel2posCMIntern = get(panel2HInt,'Position');
        set(panel2HInt,'Position',[panel2posCMIntern(1) -(rowHeightCm*no_rowsIntern-panel2posCMIntern(4)) panel2posCMIntern(3) rowHeightCm*no_rowsIntern]);
        panel2posCMIntern = get(panel2HInt,'Position');
        
        
        propTHInt = cell(no_rowsIntern,2);
        
        colWidthCmInt = panel2posCMIntern(3)*colWidth;
        gridyInt = max(0,panel2posCMIntern(4)-rowHeightCm):-rowHeightCm:0;
        
        gridxInt = panel2posCMIntern(1)+[ 0 colWidthCmInt];
        
        for row = 1:no_rowsIntern
            for col = 1:no_cols
                if ~(isempty(rowNames{row,col}) && isempty(callbackCell{row,col}))
                    
                    % check if columns to the right is
                    % empty. Stretch this column if this is
                    % the case.
                    emptcols = find(cellfun(@(x) isempty(x),rowNames(row,col:end)),1,'first');
                    if ~isempty(emptcols)
                        col_effect = col + emptcols - 1;
                        colWidth_effective = sum(colWidthCmInt(col:col_effect));
                    else
                        colWidth_effective = colWidthCmInt(col);
                    end
                    propTHInt{row,col} = uicontrol(panel2HInt,'String',rowNames(row,col), ...
                        'Style','pushbutton',...
                        'Units','centimeters', ...
                        'Callback',callbackCell{row,col}, ...
                        'Position',[gridxInt(col) gridyInt(row)  colWidth_effective rowHeightCm], ...
                        'Fontsize',fontsize,'Enable','on','HorizontalAlignment','left','BackgroundColor','w');
                    uistack(propTHInt{row,col},'bottom')
                    uistack(panel2HInt,'bottom')
                                        
                    if ~isempty(tooltipCell)
                        set(propTHInt{row,col},'TooltipString',tooltipCell{row,col})
                    end
                end
            end
        end
    end


    function helloworld(src,event)
        h=msgbox('Hello!');
        waitfor(h);
        
    end
end