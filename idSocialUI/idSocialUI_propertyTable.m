function varargout = idSocialUI_propertyTable(parentH,dataCell,descrFlag)
varargout{1} = [];%dataCell;
fontsize = 12;
propTH = [];
warnH = [];
pathChooseString = 'Choose Path';
allOk = true;
if nargin<1 || isempty(parentH)
    parentH = ...
        figure('Units','centimeters','Position',[16 2 19 17],'Name','Options', 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off','WindowStyle','normal');
end

if nargin<3 || isempty(descrFlag)
    descrFlag = false;
end


if nargin<2 || isempty(dataCell)
    if descrFlag
        dataCell = {'Option Name', 'Rendered', 'Description'};
        
        no_cols = 3;
    else
        dataCell = [strcat(repmat({'Option_Name_'},[60,1]),cellstr(num2str((1:60)'))), num2cell(1:60)'];
        dataCell{3,2} = [];
        no_cols = 2;
    end
else
    no_cols = 2;
end

tooltipCell = [];
if isstruct(dataCell)
    if size(dataCell,2)>1
        tooltipCell  = struct2cell(dataCell(2));
        dataCell = dataCell(1);
    end
    fieldnms = fieldnames(dataCell);
    fvals = struct2cell(dataCell);
    if descrFlag
        dataCell2 = cell(size(fieldnms,1),3);
        
    else
        dataCell2 = cell(size(fieldnms,1),2);
        
    end
    dataCell2(:,1) = fieldnms;
    dataCell2(:,2) = fvals;
    dataCell = dataCell2;
end
varType = check_VarType(dataCell(:,2));

no_rows = size(dataCell,1);
outval = cell(no_rows,1);



last_figpos = -1;



panelH = uipanel(parentH,'Units','normalized',...
    'Position', [0 .1 1 .9]);
panel2H = uipanel(panelH,'Units','normalized',...
    'Position', [0 0 0.96 1]);

% if verLessThan('matlab', '8.4')
%     set(parentH,'ResizeFcn',{@changeSizFunc,panelH,panel2H});
% else
%     set(parentH,'SizeChangedFcn',{@changeSizFunc,panelH,panel2H});
% end
set(parentH,'Resize','off','CloseRequestFcn',@closeOptionWindow);

sliderH = uicontrol('Style','Slider','Parent',panelH,...
    'Units','normalized','Position',[0.96 0 0.04 1],...
    'Value',1,'Callback',{@slider_callback1,panel2H});

set(gcf,'WindowScrollWheelFcn', {@mouseScroll,panel2H,sliderH});

% propTH = cell(size(dataCell,1),no_cols);

% gridy = 1-rowHeight:-rowHeight:1-rowHeight-no_rows*rowHeight;
rowHeightCm =.7;
% gridy = 1-rowHeightCm:-rowHeightCm:0;
propTH = renderTable(panel2H,dataCell(:,1),varType,rowHeightCm);


% Scale normalized to centimeter
set(panelH,'Units','centimeters');
panel1posCM = get(panelH,'Position');
set(panelH,'Units','normalized');

set(panel2H,'Units','centimeters');
panel2posCM = get(panel2H,'Position');
set(panel2H,'Position',[panel2posCM(1) -(rowHeightCm*no_rows-panel1posCM(4)) panel2posCM(3) rowHeightCm*no_rows]);
panel2posCM = get(panel2H,'Position');

% colWidthCm = panel2posCM(3)/2;
% gridy = panel2posCM(4)-rowHeightCm:-rowHeightCm:panel2posCM(2);
% gridx = panel2posCM(1):colWidthCm:panel2posCM(1)+panel2posCM(3);


% panelOkH = uipanel(parentH,'Units','normalized',...
%     'Position', [0 0 1 .1]);
fakeH = uicontrol(parentH,'Units','normalized','Style','pushbutton','String','', ...
    'Callback','','Fontsize',fontsize,'Enable','off'); %Use a fake push button instead of a panel because else in Matlab 2012 the scroll panel goes on top  
set(fakeH,'Position', [0 0 1 .1])
cancelH = uicontrol(parentH,'Units','normalized','Style','pushbutton','String','Cancel', ...
    'Callback','delete(gcf)','Fontsize',fontsize);
set(cancelH,'Position', [.2 0.01 .2 .08])
okH = uicontrol(parentH,'Units','normalized','Style','pushbutton','String','Ok', ...
    'Callback',@(src,evnt) okCallback(src,evnt,propTH,varType),'Enable','on','Fontsize',fontsize);
set(okH,'Position', [.6 0.01 .2 .08])


% for row = 1:no_rows
%     for col = 1:no_cols
%         if col ==1
%             propTH{row} = uicontrol(panel2H,'String',dataCell(row,col), ...
%                 'Style','edit',...
%                 'Units','centimeters', ...
%                 'Position',[gridx(col) gridy(row)  colWidthCm rowHeightCm], ...
%                 'Fontsize',fontsize,'Enable','off','HorizontalAlignment','left');
%         else
%             no_fields = varType(row).dims;
%            
%             if ~isempty(no_fields)
%                 propTH{row} = cell(no_fields);
%                 for frow = 1:no_fields(1)
%                     for fcol = 1:no_fields(2)
%                         try
%                            
%                         propTH{row}{frow,fcol} = uicontrol(panel2H, ...
%                             'Style',varType(row).renderer{frow,fcol},...
%                             'Units','centimeters', ...
%                             'Position',[gridx(col)+(fcol-1)*colWidthCm/no_fields(2) gridy(row)  colWidthCm/no_fields(2) rowHeightCm], ...
%                             'Fontsize',fontsize);
%                         if strcmpi(varType(row).renderer{frow,fcol},'edit') || strcmpi(varType(row).renderer{frow,fcol},'popupmenu')
%                             set(propTH{row}{frow,fcol},'String',varType(row).string{frow,fcol});
%                         end
%                         if strcmpi(varType(row).renderer{frow,fcol},'checkbox')
%                             extnd = get(propTH{row}{frow,fcol},'Extent');
%                             set(propTH{row}{frow,fcol},'Value',varType(row).value{frow,fcol},'HorizontalAlignment','center',...
%                                 'Position',[colWidthCm/no_fields(2)*.5+gridx(col)+(fcol-1)*colWidthCm/no_fields(2)-extnd(3) gridy(row)  colWidthCm/no_fields(2) rowHeightCm]);
%                         end
%                         catch
%                             keyboard
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

    function slider_callback1(src,eventdata,arg1)
        set(src,'Units','normalized')
        val = get(src,'Value');
        set(arg1,'Units','centimeter')
        act_pos = get(arg1,'Position');
        
        set(arg1,'Position', [act_pos(1) -val*abs(panel2posCM(2)) act_pos(3) act_pos(4)])
%         set(src,'Value',abs(act_pos(2))/abs(panel2posCM(2)))

    end
    function mouseScroll(src,eventdata, handles,arg1)
        set(handles,'Units','centimeter')
        act_pos = get(handles,'Position');
            val=0;
            if eventdata.VerticalScrollCount > 0 && act_pos(2) <0
                val = rowHeightCm;
                    act_pos(2) = min(act_pos(2)+val,0);
            elseif eventdata.VerticalScrollCount < 0 && act_pos(2) > panel2posCM(2);%-(rowHeightCm*no_rows-panel2posCM(4))
                val = -rowHeightCm;
                act_pos(2) = max(act_pos(2)+val,panel2posCM(2));
            end
        
            
            set(handles,'Position',[act_pos(1)  act_pos(2) act_pos(3) act_pos(4)])

      
            set(arg1,'Value',abs(act_pos(2))/abs(panel2posCM(2)))

    end
    
    function changeSizFunc(src,eventdata,handles1, handles)
%         set(parentH,'Resize','on');

        set(handles1,'Units','centimeters');
        panel1posCM = get(handles1,'Position');
        set(handles1,'Units','normalized');
        set(handles,'Units','centimeters');
        panel2posCM = get(handles,'Position');
        if panel1posCM(4)<=panel2posCM(4)+rowHeightCm
            set(handles,'Position',[panel2posCM(1) -(rowHeightCm*no_rows-panel2posCM(4)) panel2posCM(3) rowHeightCm*no_rows]);
            set(parentH,'Units','centimeters');
            last_figpos = get(parentH,'Position');
            set(parentH,'Units','normalized');
        else
            set(parentH,'Units','centimeters');
            figpos = get(parentH,'Position');
            set(parentH,'Position',[figpos(1) last_figpos(2) figpos(3) last_figpos(4)])
            set(parentH,'Units','normalized');
            %             set(parentH,'Resize','off');
        end
    end

    function varType = check_VarType(varCell)
        no_vars = size(varCell,1);
        varType = struct([]);
        
        for actVar = 1:no_vars
            if ~isempty(varCell{actVar})
                actVarVal = varCell{actVar};
                varSize = size(varCell{actVar});
                
                if isfloat(actVarVal) % Float/doubles
                    varType(actVar).dims = size(varCell{actVar});
                    
                    
                    for varRow = 1:varSize(1)
                        for varCol = 1:varSize(2)
                            varType(actVar).type{varRow,varCol} = 'float';
                            varType(actVar).renderer{varRow,varCol} = 'edit';
                            varType(actVar).string{varRow,varCol} = num2str(actVarVal(varRow,varCol));
                            varType(actVar).value{varRow,varCol} = actVarVal(varRow,varCol);

                        end
                    end
                    
                elseif isinteger(actVarVal) % Integer
                    varType(actVar).dims = size(varCell{actVar});
                    varType(actVar).type{1,1} = 'Integer';
                    varType(actVar).renderer{1,1} = 'edit';
                    varType(actVar).string{1,1} = num2str(actVarVal);
                    varType(actVar).value{1,1} = actVarVal;
                elseif islogical(actVarVal) % Logical
                    varType(actVar).dims = size(varCell{actVar});
                    varType(actVar).type{1,1} = 'logical';
                    varType(actVar).renderer{1,1} = 'checkbox';
                    varType(actVar).string{1,1} = num2str(actVarVal);
                    varType(actVar).value{1,1} = actVarVal;
                elseif ischar(actVarVal) && (strcmpi(actVarVal(1),'<') && strcmpi(actVarVal(end),'>')) || ...
                        (isletter(actVarVal(1)) &&  strcmpi(actVarVal(3),filesep)) % Pathname
                    varType(actVar).dims = [1 1];
                    varType(actVar).type{1,1} = 'path';
                    varType(actVar).renderer{1,1} = 'pushbutton';
                    varType(actVar).string{1,1} = pathChooseString;
                    varType(actVar).value{1,1} = '';
                elseif ischar(actVarVal) % Char
                    varType(actVar).dims = [1 1];
                    varType(actVar).type{1,1} = 'char';
                    varType(actVar).renderer{1,1} = 'edit';
                    varType(actVar).string{1,1} = actVarVal;
                    varType(actVar).value{1,1} = '';
                elseif iscell(actVarVal) % Cell
                    varType(actVar).dims = size(actVarVal);
                    % Combo list
                    if (iscellstr(actVarVal) || (iscellstr(actVarVal(1:end-1)) && iscell(actVarVal{end}))) && min(varSize)==1 && max(varSize)>1
                        varType(actVar).dims = [1 1];
                        varType(actVar).type{1,1} = 'cellstr';
                        varType(actVar).renderer{1,1} = 'popupmenu';
                        if iscell(actVarVal{end})
                            try
                            fidx = cellfun(@(x) ~isempty(x),strfind(actVarVal(1:end-1),actVarVal{end}{1}));
                            catch
                                keyboard
                            end
                            varType(actVar).value{1,1} = [actVarVal(fidx) actVarVal(~fidx)];
                            varType(actVar).string{1,1} = [actVarVal(fidx) actVarVal(~fidx)];
                        else
                            varType(actVar).value{1,1} = actVarVal;
                            varType(actVar).string{1,1} = actVarVal;
                        end
                    else % cell of numerical arrays (e.g., 2d - edges)
%                     
                            varType(actVar) = check_varTypeCell(actVarVal,varType(actVar));
                            
%                        
                        %                         for varRow = 1:varSize(1)
                        %                             for varCol = 1:varSize(2)
                        %                                 actVarVal2 = varCell{actVar}{varRow,varCol};
                        %     %                                 varSize = size(varCell{actVar});
                        %                                 if mixedType
                        %
                        %                                     for varRow2 = 1:size(actVarVal2,1)
                        %                                         for varCol2 = 1:size(actVarVal2,2)
                        %                                             varType(actVar).type{varRow,varCol}{varRow2,varCol2} = 'float';
                        %                                             varType(actVar).renderer{varRow,varCol}{varRow2,varCol2} = 'edit';
                        %                                             varType(actVar).string{varRow,varCol}{varRow2,varCol2} = num2str(actVarVal2(varRow2,varCol2));
                        %                                             varType(actVar).value{varRow,varCol}{varRow2,varCol2} = actVarVal2(varRow2,varCol2);
                        %
                        %                                         end
                        %                                     end
                        %                                 elseif isfloat(actVarVal2) % Float/doubles
                        %                                     varType(actVar).type{varRow,varCol} = 'float';
                        %                                     varType(actVar).renderer{varRow,varCol} = 'edit';
                        %                                     varType(actVar).string{varRow,varCol} = num2str(actVarVal2);
                        %                                     varType(actVar).value{varRow,varCol} = actVarVal2;
                    %
                    %
                    %
                    %                                 elseif isinteger(actVarVal2) % Integer
                    %                                     varType(actVar).type{varRow,varCol} = 'Integer';
                    %                                     varType(actVar).renderer{varRow,varCol} = 'edit';
                    %                                     varType(actVar).string{varRow,varCol} = num2str(actVarVal2);
                    %                                     varType(actVar).value{varRow,varCol} = actVarVal2;
                    %                                 elseif islogical(actVarVal2) % Logical
                    %                                     varType(actVar).type{varRow,varCol} = 'logical';
                    %                                     varType(actVar).renderer{varRow,varCol} = 'checkbox';
                    %                                     varType(actVar).value{varRow,varCol} = actVarVal2;
                    %                                     varType(actVar).string{varRow,varCol} = num2str(actVarVal2);
                    %                                 elseif ischar(actVarVal2) % Char
                    %                                     varType(actVar).type{varRow,varCol} = 'char';
                    %                                     varType(actVar).renderer{varRow,varCol} = 'edit';
                    %                                     varType(actVar).string{varRow,varCol} = actVarVal2;
                    %                                     varType(actVar).value{varRow,varCol} = '';
                    %                                 elseif (iscellstr(actVarVal2) || (iscellstr(actVarVal2(1:end-1)) && iscell(actVarVal2{end})))% && min(varSize)==1 && max(varSize)>1
                    % %
                    %                                     varType(actVar).type{varRow,varCol} = 'cellstr';
                    %                                     varType(actVar).renderer{varRow,varCol} = 'popupmenu';
                    %                                     if iscell(actVarVal2{end})
                    %                                         try
                    %                                             fidx = cellfun(@(x) ~isempty(x),strfind(actVarVal2(1:end-1),actVarVal2{end}{1}));
                    %                                         catch
                    %                                             keyboard
                    %                                         end
                    %                                         varType(actVar).value{varRow,varCol} = [actVarVal2(fidx) actVarVal2(~fidx)];
                    %                                         varType(actVar).string{varRow,varCol} = [actVarVal2(fidx) actVarVal2(~fidx)];
                    %                                     else
                    %                                         varType(actVar).value{varRow,varCol} = actVarVal2;
                    %                                         varType(actVar).string{varRow,varCol} = actVarVal2;
                    %                                     end
                    %                                 end
                    %
                    %                             end
                    %                         end
                end
                
                
                
                %                 elseif iscellstr(actVarVal) % Cellstr
                %
                %                 elseif isstruct(actVarVal) % Struct
            else
                varType(actVar) = [];
            end
            
            else
                varType(actVar).dims = [1 1];
                varType(actVar).type{1,1} = 'unknown';
                varType(actVar).renderer{1,1} = 'edit';
                varType(actVar).string{1,1} = '';
                varType(actVar).value{1,1} = NaN;
            end
        end
        
        
        function varType = check_varTypeCell(actVarVal,varType)
            if nargin<2 || isempty(varType)
                varType = struct;
            end
            varSizeIntern = size(actVarVal);
            mixedType = ~all(cellfun(@(x) isequal(class(x),class(actVarVal(1))),actVarVal));
            classVal = cellfun(@(x) class(x),actVarVal,'UniformOutput',false);
            isArray = cellfun(@(x) numel(x)>1,actVarVal);
            sizeArray = cellfun(@(x) numel(x),actVarVal);
            
            arrayIdx = find(isArray & ~strcmp(classVal,'cell'));
            actVarValNewSize = sum(sizeArray(isArray & ~strcmp(classVal,'cell'))) + sum(~(isArray & ~strcmp(classVal,'cell')));
            actVarValNew = cell(1,actVarValNewSize);
            actIdx = 1;
            for newIdx = 1:varSizeIntern(2)
                if any(newIdx==arrayIdx)
                    actVarValNew(actIdx:actIdx+sizeArray(newIdx)-1) = num2cell([actVarVal{newIdx}(:)]);
                    actIdx = actIdx+sizeArray(newIdx);
                else
                    actVarValNew(actIdx) = actVarVal(newIdx);
                    actIdx = actIdx+1;
                end
            end
            actVarVal = actVarValNew;
            varSizeIntern = size(actVarVal);
            
            varType.dims = [1 varSizeIntern(2)];
            
            for varRowIntern = 1:varSizeIntern(1)
                for varColIntern = 1:varSizeIntern(2)
                    actVarVal2Intern = actVarVal{varRowIntern,varColIntern};
                    %                                 varSize = size(varCell{actVar});
                    if isfloat(actVarVal2Intern) % Float/doubles
                        varType.type{varRowIntern,varColIntern} = 'float';
                        varType.renderer{varRowIntern,varColIntern} = 'edit';
                        varType.string{varRowIntern,varColIntern} = num2str(actVarVal2Intern);
                        varType.value{varRowIntern,varColIntern} = actVarVal2Intern;
                        
                        
                        
                    elseif isinteger(actVarVal2Intern) % Integer
                        varType.type{varRowIntern,varColIntern} = 'Integer';
                        varType.renderer{varRowIntern,varColIntern} = 'edit';
                        varType.string{varRowIntern,varColIntern} = num2str(actVarVal2Intern);
                        varType.value{varRowIntern,varColIntern} = actVarVal2Intern;
                    elseif islogical(actVarVal2Intern) % Logical
                        varType.type{varRowIntern,varColIntern} = 'logical';
                        varType.renderer{varRowIntern,varColIntern} = 'checkbox';
                        varType.value{varRowIntern,varColIntern} = actVarVal2Intern;
                        varType.string{varRowIntern,varColIntern} = num2str(actVarVal2Intern);
                    elseif ischar(actVarVal2Intern) % Char
                        varType.type{varRowIntern,varColIntern} = 'char';
                        varType.renderer{varRowIntern,varColIntern} = 'edit';
                        varType.string{varRowIntern,varColIntern} = actVarVal2Intern;
                        varType.value{varRowIntern,varColIntern} = '';
                    elseif (iscellstr(actVarVal2Intern) || (iscellstr(actVarVal2Intern(1:end-1)) && iscell(actVarVal2Intern{end})))% && min(varSize)==1 && max(varSize)>1
                        %
                        varType.type{varRowIntern,varColIntern} = 'cellstr';
                        varType.renderer{varRowIntern,varColIntern} = 'popupmenu';
                        if iscell(actVarVal2Intern{end})
                            try
                                fidxIntern = cellfun(@(x) ~isempty(x),strfind(actVarVal2Intern(1:end-1),actVarVal2Intern{end}{1}));
                            catch
                                keyboard
                            end
                            varType.value{varRowIntern,varColIntern} = [actVarVal2Intern(fidxIntern) actVarVal2Intern(~fidxIntern)];
                            varType.string{varRowIntern,varColIntern} = [actVarVal2Intern(fidxIntern) actVarVal2Intern(~fidxIntern)];
                        else
                            varType.value{varRowIntern,varColIntern} = actVarVal2Intern;
                            varType.string{varRowIntern,varColIntern} = actVarVal2Intern;
                        end
                    end
                    
                end
            end
            
        end
    end

    

    function propTHInt = renderTable(panel2HInt,rowNames,varTypeIntern,rowHeightCm)
        
        no_rowsIntern = size(rowNames,1);
        no_colsIntern = 2;
        set(panel2HInt,'Units','centimeters');
        panel2posCMIntern = get(panel2HInt,'Position');
        set(panel2HInt,'Position',[panel2posCMIntern(1) -(rowHeightCm*no_rowsIntern-panel2posCMIntern(4)) panel2posCMIntern(3) rowHeightCm*no_rowsIntern]);
        panel2posCMIntern = get(panel2HInt,'Position');
        
        propTHInt = cell(no_rowsIntern,2);
        
        colWidthCmInt = panel2posCMIntern(3)/2;
        gridyInt = panel2posCMIntern(4)-rowHeightCm:-rowHeightCm:0;%panel2posCMIntern(2);
        gridxInt = panel2posCMIntern(1):colWidthCmInt:panel2posCMIntern(1)+panel2posCMIntern(3);
        
        for row = 1:no_rowsIntern
            for col = 1:no_colsIntern
                if col ==1
                    propTHInt{row,1} = uicontrol(panel2HInt,'String',rowNames(row), ...
                        'Style','edit',...
                        'Units','centimeters', ...
                        'Position',[gridxInt(col) gridyInt(row)  colWidthCmInt rowHeightCm], ...
                        'Fontsize',fontsize,'Enable','off','HorizontalAlignment','left');
                    if ~isempty(tooltipCell)
                        set(propTHInt{row,1},'TooltipString',tooltipCell{row})
                    end
                else
                    no_fields = varTypeIntern(row).dims;
                    
                    if ~isempty(no_fields)
                        propTHInt{row,2} = cell(no_fields);
                        for frow = 1:no_fields(1)
                            for fcol = 1:no_fields(2)
                                try
%                                     if iscell(varTypeIntern(row).renderer{frow,fcol}) % Even one more nesting level
%                                         no_fields2 = size(varTypeIntern(row).renderer{frow,fcol});
%                                         for frow2 = 1:no_fields2(1)
%                                             for fcol2 = 1:no_fields2(2)
%                                                 propTHInt{row,2}{frow,fcol}{frow2,fcol2} = uicontrol(panel2HInt, ...
%                                                     'Style',varTypeIntern(row).renderer{frow,fcol}{frow2,fcol2},...
%                                                     'Units','centimeters', ...
%                                                     'Position',[gridxInt(col)+(fcol-1)*colWidthCmInt/no_fields(2)+(fcol2-1)*colWidthCmInt/no_fields(2)/no_fields2(2) ...
%                                                     gridyInt(row)  colWidthCmInt/no_fields(2)/no_fields2(2) rowHeightCm],'BackgroundColor','w',... ...
%                                                     'Fontsize',fontsize,'Callback',{@checkInput,varTypeIntern,row,frow,fcol});
%                                             end
%                                         end
%                                     else
                                        propTHInt{row,2}{frow,fcol} = uicontrol(panel2HInt, ...
                                            'Style',varTypeIntern(row).renderer{frow,fcol},...
                                            'Units','centimeters', ...
                                            'Position',[gridxInt(col)+(fcol-1)*colWidthCmInt/no_fields(2) gridyInt(row)  colWidthCmInt/no_fields(2) rowHeightCm],'BackgroundColor','w',... ...
                                            'Fontsize',fontsize,'Callback',{@checkInput,varTypeIntern,row,frow,fcol});
                                 
%                                     end
                                catch
                                    keyboard
                                end
                                if iscell(varTypeIntern(row).renderer{frow,fcol})
                                    no_fields2 = size(varTypeIntern(row).renderer{frow,fcol});
                                    for frow2 = 1:no_fields2(1)
                                        for fcol2 = 1:no_fields2(2)
                                            if strcmpi(varTypeIntern(row).renderer{frow,fcol}{frow2,fcol2},'edit') || strcmpi(varTypeIntern(row).renderer{frow,fcol}{frow2,fcol2},'popupmenu')
                                                set(propTHInt{row,2}{frow,fcol}{frow2,fcol2},'String',varTypeIntern(row).string{frow,fcol}{frow2,fcol2});
                                            end
                                            if strcmpi(varTypeIntern(row).renderer{frow,fcol}{frow2,fcol2},'pushbutton')
                                                set(propTHInt{row,2}{frow,fcol}{frow2,fcol2},'String',varTypeIntern(row).string{frow,fcol}{frow2,fcol2},'Callback',@getpath);
                                            end
                                            if strcmpi(varTypeIntern(row).renderer{frow,fcol}{frow2,fcol2},'checkbox')
                                                extnd = get(propTHInt{row,2}{frow,fcol}{frow2,fcol2},'Extent');
                                                fakebackground = uicontrol(panel2HInt, ...
                                                    'Style','edit',...
                                                    'Enable','off',...
                                                    'Units','centimeters', ...
                                                    'Position',[gridxInt(col)+(fcol-1)*colWidthCmInt/no_fields(2)+(fcol2-1)*colWidthCmInt/no_fields(2)/no_fields2(2) ...
                                                    gridyInt(row)  colWidthCmInt/no_fields(2)/no_fields2(2) rowHeightCm], ...
                                                    'Fontsize',fontsize);
                                                uistack(propTHInt{row,2}{frow,fcol}{frow2,fcol2},'top')
                                                set(propTHInt{row,2}{frow,fcol}{frow2,fcol2},'Value',varTypeIntern(row).value{frow,fcol}{frow2,fcol2},'HorizontalAlignment','center','BackgroundColor',get(fakebackground,'BackgroundColor'),...
                                                    'Position',[colWidthCmInt/no_fields(2)*.5+gridxInt(col)+(fcol-1)*colWidthCmInt/no_fields(2)-extnd(3) gridyInt(row)+.1  colWidthCmInt/no_fields(2)-.1 rowHeightCm-.2]);
                                            end
                                            if ~isempty(tooltipCell)
                                                set(propTHInt{row,2}{frow,fcol}{frow2,fcol2},'TooltipString',tooltipCell{row})
                                            end
                                        end
                                    end
                                else
                                    if strcmpi(varTypeIntern(row).renderer{frow,fcol},'edit') || strcmpi(varTypeIntern(row).renderer{frow,fcol},'popupmenu')
                                        set(propTHInt{row,2}{frow,fcol},'String',varTypeIntern(row).string{frow,fcol});
                                    end
                                    if strcmpi(varTypeIntern(row).renderer{frow,fcol},'pushbutton')
                                        set(propTHInt{row,2}{frow,fcol},'String',varTypeIntern(row).string{frow,fcol},'Callback',@getpath);
                                    end
                                    if strcmpi(varTypeIntern(row).renderer{frow,fcol},'checkbox')
                                        extnd = get(propTHInt{row,2}{frow,fcol},'Extent');
                                        fakebackground = uicontrol(panel2HInt, ...
                                            'Style','edit',...
                                            'Enable','off',...
                                            'Units','centimeters', ...
                                            'Position',[gridxInt(col)+(fcol-1)*colWidthCmInt/no_fields(2) gridyInt(row)  colWidthCmInt/no_fields(2) rowHeightCm],...
                                            'Fontsize',fontsize);
                                        uistack(propTHInt{row,2}{frow,fcol},'top')
                                        set(propTHInt{row,2}{frow,fcol},'Value',varTypeIntern(row).value{frow,fcol},'HorizontalAlignment','center','BackgroundColor',get(fakebackground,'BackgroundColor'),...
                                            'Position',[colWidthCmInt/no_fields(2)*.5+gridxInt(col)+(fcol-1)*colWidthCmInt/no_fields(2)-extnd(3) gridyInt(row)+.1  colWidthCmInt/no_fields(2)-.1 rowHeightCm-.2]);
                                    end
                                    if ~isempty(tooltipCell)
                                        set(propTHInt{row,2}{frow,fcol},'TooltipString',tooltipCell{row})
                                    end
                                end
                            end
                        end
                    end
                end
                
            end
            
        end
        
        function getpath(src,evnt)
            newpath = uigetdir;
            set(src,'String',newpath,'TooltipString',newpath);
        end
        
        function checkInput(src,evnt,varTypeIntern,actRow,actFrow,actFcol)
            inp = get(src,'String');
            actType = varTypeIntern(actRow).type{actFrow,actFcol};
            allOk = true;
            if strcmpi(actType,'float') || strcmpi(actType,'integer')
                if isnan(str2double(inp))
                    allOk = false;
                    waitfor(warndlg('Please enter a number.','Wrong Parameter Type','modal'));
                   
                    

                    if ishandle(src)
                        set(src,'String',varTypeIntern(actRow).string{actFrow,actFcol});
                        uicontrol(src)
                    end
                    
                    
                end
            elseif strcmpi(actType,'char')
            end
            
        end
    end

    function okCallback(src,evnt,propTHInt,varTypeIntern)
        no_rowsInt = size(varTypeIntern,2);
%         outval = cell(no_rowsInt,1);
        for rowInt =  1:no_rowsInt
            %             if strcmpi(varTypeIntern(rowInt).type,'float'
            no_fields = varTypeIntern(rowInt).dims;
            if ~isempty(no_fields) && all(no_fields) > 0
                for frow = 1:no_fields(1)
                    for fcol = 1:no_fields(2)
                        if iscell(varTypeIntern(rowInt).type{frow,fcol}) % Mixed type, e.g., float and cellstr=popup
                            
                        else
                            act_string = get(propTHInt{rowInt,2}{frow,fcol},'String');
                            act_val = get(propTHInt{rowInt,2}{frow,fcol},'Value');
                            
                            if strcmpi(varTypeIntern(rowInt).type{frow,fcol},'unknown')
                                if ~isempty(act_string) && isnan(str2double(act_string))
                                    outval{rowInt}{frow,fcol} = act_string;
                                elseif ~isempty(act_string)
                                    outval{rowInt}(frow,fcol)= str2double(act_string);
                                elseif isempty(act_string)
                                    outval{rowInt}{frow,fcol} = '';
                                end
                            elseif strcmpi(varTypeIntern(rowInt).type{frow,fcol},'char')
                                outval{rowInt}{frow,fcol} = act_string;
                            elseif strcmpi(varTypeIntern(rowInt).type{frow,fcol},'path')
                                %                             if strcmpi(act_string(1),'<') && strcmpi(act_string(end),'>')
                                %                                 act_string = '';
                                %                             end
                                dirOk = false;
                                if exist(act_string,'dir')~=7 && ~strcmpi(pathChooseString,act_string)
                                    
                                    try
                                        dirOk = mkdir(act_string);
                                    catch
                                    end
                                end
                                
                                if exist(act_string,'dir')~=7 || ~dirOk
                                    outval{rowInt}{frow,fcol} = '';
                                else
                                    outval{rowInt}(frow,fcol) = act_string;
                                end
                            elseif strcmpi(varTypeIntern(rowInt).type{frow,fcol},'float')
                                if all(strcmpi(varTypeIntern(rowInt).type,'float'))
                                    outval{rowInt}(frow,fcol) = str2double(act_string);
                                else
                                    outval{rowInt}{frow,fcol} = str2double(act_string);
                                end
                            elseif strcmpi(varTypeIntern(rowInt).type{frow,fcol},'cellstr')
                                if isequal(no_fields,[1 1])
                                    outval{rowInt} = {act_string{:} act_string(act_val)};
                                else
                                    outval{rowInt}{frow,fcol} = {act_string{:} act_string(act_val)};
                                end
                                %                             act_string(act_val)
                            elseif strcmpi(varTypeIntern(rowInt).type{frow,fcol},'logical')
                                outval{rowInt}(frow,fcol) = logical(get(propTHInt{rowInt,2}{frow,fcol},'Value'));
                            end
                        end
                    end
                end
            elseif no_fields == 0
                act_string = get(propTHInt{rowInt,2},'String');
                testval = str2double(act_string);
                if isnan(testval)
                    outval{rowInt} = act_string;
                else
                    outval{rowInt} = testval;
                end
            end
            
            
            
        end
        out1 = cell2struct(outval,cellfun(@(x) strrep(x,' ','_'),dataCell(:,1),'UniformOutput',false));
        if ~isempty(tooltipCell)
            out1 =  cat(2,out1,cell2struct(tooltipCell,dataCell(:,1)));
        end
        
        varargout{1} = out1;
        
        %         if ~isempty(warnH)
        %             uiwait(warnH)
        %         end
        %         if allOk
%         uiwait
        if allOk 
            close(parentH)
        else
            allOk = true;
        end
        %         end
    end

    function closeOptionWindow(src,~)
        
        if isempty(varargout{1})
            out1 = cell2struct(dataCell(:,2),dataCell(:,1));
            if ~isempty(tooltipCell)
                out1 =  cat(2,out1,cell2struct(tooltipCell,dataCell(:,1)));
            end
            varargout{1} = out1;
        end
        delete(parentH) 
    end

uiwait
    
end