function idSocialUI


gui.table_fontsize = 12;
gui.pathColumn=4;
gui.groupColumn=1;
gui.subsetColumn=2;
gui.trialColumn=3;
gui.preOptColumn=5;
gui.removeColumn=6;
gui.indivColumn=7;
gui.framerateColumn=8;
gui.bodylengthColumn=9;

gui.fontsize_button = 14;

% Add path
act_path=mfilename('fullpath');

% C = strsplit(act_path,'\');
C = textscan(act_path,'%s','Delimiter',[filesep filesep]);
C=C{1};
idS = strfind(C(1:end-1),'idSocial');
idSocIdx = find(cellfun(@(x) ~isempty(x),idS));
act_path = strcat(C(1:idSocIdx),filesep);
act_path = [act_path{:}];

% idSocIdx=strfind(act_path,'\idSocial');
% act_path=act_path(1:idSocIdx+length('idSocial')+1);
try
    cd(act_path)
catch le
    %     le=lasterror;
    if strcmp(le.identifier,'MATLAB:cd:NonExistentDirectory')
        error([mfilename ': Failed to find idSocial directory.' act_path ' nonexistent or not a directory.'])
    end
end
addpath(genpath(act_path))
%%%

gui.tooltips_on = true;

gui.figure_windows = [];
[gui.sectionNames, gui.funcCell]=idSocial_auxiliaries_guiSectionsAndFunctions;
gui.no_sections=size(gui.sectionNames,1);

fh = figure('Units','normalized','Position',[.1 .1 .85 .7],'Name','idSocial', 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off');

if ~isprop(fh,'SizeChangedFcn') && isprop(fh,'ResizeFcn')
    set(fh,'ResizeFcn',@(src,evnt) trajectoryTableSiz(src,evnt));
else
    set(fh,'SizeChangedFcn',@(src,evnt) trajectoryTableSiz(src,evnt));
end
mh = uimenu(fh,'Label','File');
savh = uimenu(mh,'Label','Save as...',...
    'Callback',@(src,evnt) saveData(src,evnt,fh,1));
savh = uimenu(mh,'Label','Save',...
    'Callback',@(src,evnt) saveData(src,evnt,fh,2));
guidata(fh,gui);
loadh = uimenu(mh,'Label','Load',...
    'Callback',@(src,evnt) loadData(src,evnt,fh));
gui=guidata(fh);

confh = uimenu(mh,'Label','Set idSocial path',...
    'Callback',@(src,evnt) setPath(src,evnt,fh));

fpos = get(fh,'Position');
% Column names and column format
columnname = {'<html><font size=+1>Group','<html><font size=+1>Subset','<html><font size=+1>Trial','<html><font size=+1>Trajectory File','<html><font size=+1>Options','<html><font size=+1>Remove Trajectory'};
columnformat = {'char','char','char','char','char','char'};

% Define the data
gui.DefaultData =    {'Add' 'Add' 'Add' '' 'Set All' '';...
    };


% gui.DefaultData =    {'<html><h1>Add</h1></html>' '<html><h1>Add</h1></html>' '<html><h1>Add</h1></html>' '' '<html><h1>Set All</h1></html>' '';...
%     };

% Create the uitable: Data preprocessing (top)


t = uitable('Data', gui.DefaultData,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [false false false false false false],...
    'ColumnWidth',{75 75 75 100 75 75},...
    'RowName',[],...
    'FontSize',gui.table_fontsize,...
    'CellSelectionCallback',@(src,evnt) cellSelect(src,evnt,fh));





% Set width and height
set(t,'Units','normalized');
pos = get(t,'Position');
set(t,'Position',[pos(1) .55 .8 .4]);%fpos(3)*.9;%t.Extent(3);
% t.Position(4) = .9;%fpos(4)*.9;%t.Extent(4);
gui.trajectoryTable = t;

set(t,'Units','pixels');
pos = get(t,'Position');
colwidth = {pos(3)/18 pos(3)/18 pos(3)/18 pos(3)*7/10 pos(3)/10 pos(3)/10};
set(gui.trajectoryTable,'ColumnWidth',{colwidth{[1:3]} pos(3)-sum([colwidth{[1:3 5:end]}]) colwidth{[5:end]}});
set(t,'Units','normalized');
pos = get(t,'Position');

% Create listbox to choose analysis section
gui.secH = uicontrol(fh,'Units','normalized','Position', [pos(1) .05 .08 .45],'Style','listbox','String',gui.sectionNames, ...
    'Callback',@(src,evnt) sectionSelect(src,evnt,fh),'FontSize',gui.table_fontsize);
set(gui.secH,'Enable','off')
% Create listbox to choose analysis section
gui.funcH = uicontrol(fh,'Units','normalized','Position', [pos(1)+.09 .05 .19 .45],'Style','listbox','String',gui.funcCell{1,2}, ...
    'Callback',@(src,evnt) functionSelect(src,evnt,fh),'FontSize',gui.table_fontsize);
set(gui.funcH,'Enable','off')
% functionContextH = uicontextmenu(gui.secH);
% functionContextHelpH = uimenu(functionContextH,'Label','Help'); %,'Callback',@setlinestyle


% Create button to choose function options
gui.funcOptsH = uicontrol(fh,'Units','normalized','Style','pushbutton','String','Options', ...
    'Callback',@(src,evnt) functionOptions(src,evnt,fh),'FontSize',gui.fontsize_button);
set(gui.funcOptsH,'Enable','off','Position', [pos(1)+.08+.2+.01 .32 .1 .07])
% Create button to choose function adv. options
gui.funcAdvOptsH = uicontrol(fh,'Units','normalized','Style','pushbutton','String','Filters', ...
    'Callback',@(src,evnt) functionAdvOptions(src,evnt,fh),'FontSize',gui.fontsize_button);
set(gui.funcAdvOptsH,'Enable','off','Position', [pos(1)+.08+.2+.01 .24 .1 .07])
% Create help button
gui.funcHelpH = uicontrol(fh,'Units','normalized','Style','pushbutton','String','Help', ...
    'Callback',@(src,evnt) funcHelp(src,evnt,fh),'FontSize',gui.fontsize_button);
set(gui.funcHelpH,'Enable','off','Position', [pos(1)+.08+.2+.01 .16 .1 .07])



gui.selectedSection = 1;
gui.selectedFunc = 1;


% Create the uitable: Choose data (plot_mode)
gui.DataSetTableColumnNames = {'Set 1'};
gui.DataSetTableDefaultData =    {'<html><body bgcolor="#808080">Add</body></html>';' '};
gui.DataSetTable = uitable('Data', gui.DataSetTableDefaultData,...
    'ColumnName', gui.DataSetTableColumnNames,...
    'ColumnEditable', [false false false false false false],...
    'ColumnWidth',{130},...
    'RowName',[],...
    'FontSize',gui.table_fontsize);

set(gui.DataSetTable,'CellSelectionCallback',@(src,evnt) cellSelectDataSet(src,evnt,gui.DataSetTable)...
    );

if gui.tooltips_on
    s = sprintf('"Add" to select and start analysis.\n"Remove" to delete existing analysis.\nSelect specific row (Results,Indiv. Results, etc.) and click on "Display" or "Plot" to show results.');
    set(gui.DataSetTable,'TooltipString',s);
end

set(gui.DataSetTable,'Units','normalized');
pos = get(gui.DataSetTable,'Position');
set(gui.DataSetTable,'Position',[pos(1)+.08+.2+.01+.1+.01 .05 .2 .45]);
set(gui.DataSetTable,'Enable','off')

% Center text
jscrollpane = findjobj(gui.DataSetTable);
jTable = jscrollpane.getViewport.getView;
cellStyle = jTable.getCellStyleAt(0,0);
cellStyle.setHorizontalAlignment(cellStyle.CENTER);
% Table must be redrawn for the change to take affect
jTable.repaint;

gui.expTree = 0;
gui.recentPath = '';


gui.FigureTableTooltips = [];
gui.FigureTableDataDefault = {'New Figure' ''};
gui.FigureTableCallbackCell = {{@newFigure,fh} ''};

[gui.FigureTableH, gui.FigureTableButtonsH] = ...
    idSocialUI_buttonTable(fh,[pos(1)+.08+.2+.01+.1+.01+.2+.01 .05 .2-.01 .45],gui.FigureTableDataDefault,gui.FigureTableCallbackCell);

set([gui.FigureTableButtonsH{:}],'Enable','off')

% "Preprocess" Button

gui.preprocessButtonH=uicontrol('Style','pushbutton' ,'String','Preprocess', ...
    'Units', 'normalized','Position', [.85 .8 .1 .07],'FontSize',gui.fontsize_button,...
    'Callback',@(src,evnt) preProcess(src,evnt,fh));
set(gui.preprocessButtonH,'Enable','off')

% "Copy to workspace" Button
gui.export2workspaceButtonH=uicontrol('Style','pushbutton' ,'String','Export', ...
    'Units', 'normalized','Position', [.85 .32 .1 .07],'FontSize',gui.fontsize_button,...
    'Callback',@(src,evnt) export2workspace(src,evnt,fh));
set(gui.export2workspaceButtonH,'Enable','off')

% "Display" Button
gui.displayResultsButtonH=uicontrol('Style','pushbutton' ,'String','Display', ...
    'Units', 'normalized','Position', [.85 .24 .1 .07],'FontSize',gui.fontsize_button,...
    'Callback',@(src,evnt) displayResults(src,evnt,fh));
set(gui.displayResultsButtonH,'Enable','off')

% "Plot" Button
gui.plotResultsButtonH=uicontrol('Style','pushbutton' ,'String','Plot', ...
    'Units', 'normalized','Position', [.85 .16 .1 .07],'FontSize',gui.fontsize_button,...
    'Callback',@(src,evnt) plotResults(src,evnt,fh));
set(gui.plotResultsButtonH,'Enable','off')

% "New figure" Button
% gui.newFigureButtonH=uicontrol('Style','pushbutton' ,'String','New figure', ...
%     'Units', 'normalized','Position', [.85 .08 .1 .07],'FontSize',gui.fontsize_button,...
%     'Callback',@(src,evnt) newFigure(src,evnt,fh));
% set(gui.newFigureButtonH,'Enable','off')


defPreOpts = idSocial_loadData;
defPreOpts = rmfield(defPreOpts,'act_method');
gui.defPreOpts = defPreOpts;
gui.first_preprocess = true;
gui.plotResultsData = 3;
gui.selectedSet = 1;

% filename = 'D:\~idSocialTemp\NewProject_20170330T160850\input_data.mat';
guidata(fh,gui);

% load_inputData([],[],fh,filename)

end

%%%%%%%%%%%%%%%%%%%%%% Callbacks  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_inputData(src,evnt,fh,filename)
gui=guidata(fh);
if exist(filename,'file')==2
    inpload = load(filename);
    fn = fieldnames(inpload);
    good_idsocial = false(1,numel(fn));
    for fl = 1:numel(fn)
        if isstruct(inpload.(fn{fl})) && isfield(inpload.(fn{fl}),'movementdata') && ...
                isfield(inpload.(fn{fl}),'info') && isfield(inpload.(fn{fl}),'options')
            good_idsocial(fl) = true;
        end
    end
    if ~any(good_idsocial) 
        warndlg([filename ' does not contain a valid idSocial file.'])
        return
    else
        gui.input_data = inpload.(fn{good_idsocial});
    end
    
    % Set trajectory table
    
    if isfield(gui.input_data(1,1,1).info,'trajectory_origLocation')
        trList = gui.input_data(1,1,1).info.trajectory_origLocation;
    else
        warndlg('Could not find the original trajectory locations. Show preprocessed files instead.')
        trPre = gui.input_data.movementdata;
        for gr = 1:size(trPre,1)
            for sb = 1:size(trPre,2)
                for tr = 1:size(trPre,3)
                    if ~isempty(trPre{gr,sb,tr})
                        trList{gr}{sb}{tr} = trPre{gr,sb,tr};
                    end
                end
            end
        end
        
        
    end
    no_groups = numel(trList);
    
    
    
    datDef = get(gui.trajectoryTable,'Data');
    tData = cell(sum(~isnan([gui.input_data.info.bodylength_in_pixels(:)]))+1,size(datDef,2));
    row = 1;
    for gr = 1:no_groups
        no_subset = numel(trList{gr});
        for sb = 1:no_subset
            no_trials = numel(trList{gr}{sb});
            for tr = 1:no_trials
                if ~isempty(trList{gr}{sb}{tr})
                    tData{row,1} = gr;
                    tData{row,2} = sb;
                    tData{row,3} = tr;

                    tData{row,4} = trList{gr}{sb}{tr};
                    
                    colName = get(gui.trajectoryTable,'Columnname');
                    colName{gui.framerateColumn} = '<html><font size=+1>Framerate';
                    tData{row,gui.framerateColumn}=gui.input_data.info.framerate(gr,sb,tr);
                    
                    colName{gui.indivColumn} = '<html><font size=+1># Indivs.';
                    tData{row,gui.indivColumn}=gui.input_data.info.no_focals(gr,sb,tr);
                    colName{gui.bodylengthColumn} = '<html><font size=+1>Bodylength [PXL]';
                    tData{row,gui.bodylengthColumn}=gui.input_data.info.bodylength_in_pixels(gr,sb,tr);
                    
                    tData{row,gui.preOptColumn} = 'Set';
                    tData{row,gui.removeColumn} = 'Remove';
                    
                    gui.trajectoryList{gr}{sb}{tr} = trList{gr}{sb}{tr};
                    gui.loadOptsNoCombos{gr}{sb}{tr} = rmfield(gui.input_data.options{gr}{sb}{tr},'act_method');
                    gui.loadOpts{gr}{sb}{tr} = RestoreComboOptions(gui.loadOptsNoCombos{gr}{sb}{tr},gui.defPreOpts);
                    
                    row = row + 1;
                end
            end
        end
    end
    tData(end,1:size(gui.DefaultData,2)) = gui.DefaultData;
                   
    
    oldUnits = get(gui.trajectoryTable,'Units');
    set(gui.trajectoryTable,'Data',tData,'Columnname',colName)
    set(gui.trajectoryTable,'Units','pixels');
    tabPos = get(gui.trajectoryTable,'Position');
    colwidth = {tabPos(3)/18 tabPos(3)/18 tabPos(3)/18 tabPos(3)*7/10 tabPos(3)/10 tabPos(3)/10 tabPos(3)/10 tabPos(3)/10 tabPos(3)/10};
    set(gui.trajectoryTable,'ColumnWidth',{colwidth{[1:3]} tabPos(3)-sum([colwidth{[1:3 5:end]}]) colwidth{[5:end]}});
    set(gui.trajectoryTable,'Units',oldUnits);

    guidata(fh,gui);

    set(gui.trajectoryTable,'CellSelectionCallback',@(src,evnt) cellSelect(src,evnt,fh)...
    );
    set(gui.preprocessButtonH,'Enable','on')
    
    guidata(fh,gui);

    %% Set Section and Function (=1)
    sectionSelect([],[],fh)
    
    functionSelect([],[],fh)
    gui=guidata(src);

    set(gui.secH,'Enable','on')
    set(gui.funcH,'Enable','on')
    
    guidata(fh,gui);

    updateDataSetTable(fh)
    gui=guidata(src);
%     updateFigTable(fh)
    
%     set(gui.export2workspaceButtonH,'Enable','on')
%     set(gui.displayResultsButtonH,'Enable','on')
%     set(gui.plotResultsButtonH,'Enable','on')
end
guidata(fh,gui);

end

function trajectoryTableSiz(src,evnt)
gui=guidata(src);

oldUnits = get(src,'Units');
set(src,'Units','pixels');
figPos = get(src,'Position');
set(src,'Units',oldUnits);


if isfield(gui,'trajectoryTable')
    dat = get(gui.trajectoryTable,'Data');
    oldUnits = get(gui.trajectoryTable,'Units');
    set(gui.trajectoryTable,'Units','pixels');
    tabPos = get(gui.trajectoryTable,'Position');
    if size(dat,2)<=6
        colwidth = {tabPos(3)/18 tabPos(3)/18 tabPos(3)/18 tabPos(3)*7/10 tabPos(3)/10 tabPos(3)/10};
        
    else
        colwidth = {tabPos(3)/18 tabPos(3)/18 tabPos(3)/18 tabPos(3)*7/10 tabPos(3)/10 tabPos(3)/10 tabPos(3)/10 tabPos(3)/10 tabPos(3)/10};
    end
    set(gui.trajectoryTable,'ColumnWidth',{colwidth{[1:3]} tabPos(3)-sum([colwidth{[1:3 5:end]}]) colwidth{[5:end]}});
    set(gui.trajectoryTable,'Units',oldUnits);
end
end


function cellSelect(src,evnt,fh)
gui=guidata(fh);
data = get(src,'Data');
rows = evnt.Indices(:,1);
cols = evnt.Indices(:,2);
set(src,'UserData',evnt.Indices)
no_cols = size(get(src,'Data'),2);
no_rows = size(get(src,'Data'),1);

if numel(rows)==1 && numel(cols)==1
    
    if (cols==gui.groupColumn || cols==gui.subsetColumn || cols==gui.trialColumn) && rows==no_rows
        pnl = openFolderFileChooser(src,evnt,fh);
        waitfor(pnl)%(dummyFig)
        if ishandle(fh) % In case window is closed with menu open
            gui=guidata(fh);
            
            PathName = [];
            IndexList = [];
            if isfield(gui,'openFolderFileChooserList') && ~isempty(gui.openFolderFileChooserList) && ...
                    ~(numel(gui.openFolderFileChooserList)==1 && strcmpi(gui.openFolderFileChooserList,'canceled'))
                
                
                PathName = gui.openFolderFileChooserList;
                if ~(ischar(PathName) && strcmp(PathName,'canceled')) && size(PathName,2)>1 % If coming from "folder tree reconstruction"
                    IndexList = cell2mat(PathName(:,1:3));
                    PathName = PathName(:,4);
                end
                gui.openFolderFileChooserList = [];
            elseif isfield(gui,'openFolderFileChooserList')  && numel(gui.openFolderFileChooserList)==1 && strcmpi(gui.openFolderFileChooserList,'canceled')
                if isfield(gui,'recentPath') && ~isempty(gui.recentPath)
                    PathName = gui.recentPath;
                end
                gui.openFolderFileChooserList = [];
            elseif isfield(gui,'openFolderFileChooserList') && ~isempty(gui.openFolderFileChooserList)
                warndlg('Sorry, could not open the location.')
            end
            
            data_dummy = [];
            set(src,'Data',data_dummy);
            set(src,'Data',data);

            
            % Check input
            if ~(ischar(PathName) && strcmp(PathName,'canceled'))
                good_file = false(size(PathName,1),1);
                for act_file = 1:size(PathName,1)
                    good_file(act_file) = strcmp(PathName{act_file}(end-2:end),'mat');
                    good_file(act_file) = good_file(act_file) && exist(PathName{act_file},'file')==2;
                    s = whos('-file',PathName{act_file});
                    % No datosegm:
                    good_file(act_file) = good_file(act_file) && ~any(cellfun(@(x) strcmpi(x,'datosegm'),{s.name}));
                    % Exactly one variable with ndims==3 in the struct
                    good_file(act_file) = good_file(act_file) && sum(cellfun(@(x) size(x,2),{s.size})==3)==1;
                    if good_file(act_file)
                        % Check dimensions of variable with ndims==3
                        good_file(act_file) = good_file(act_file) && cellfun(@(x) x(1)>x(2) & x(1)>x(3) & (x(3)==2 | x(3)==3),{s(cellfun(@(x) size(x,2),{s.size})==3).size});
                    end
                    
                end
                if ~any(good_file)
                    warndlg('Sorry, no trajectories found in the selected file(s).')
                end
                PathName = PathName(good_file);
                no_files_act = size(PathName,1);
                if ~isempty(IndexList)
                    IndexList = IndexList(good_file,:);
                end
                
                % END check input
                
                %         [FileName,PathName] = uigetfile(gui.recentPath);
                if ~isempty(PathName) && isempty(IndexList)%%&& ~(isnumeric(PathName)&& PathName==0) && (strcmpi(FileName,'trajectories.mat') || strcmpi(FileName,'trayectorias.mat'))
                    %             gui.recentPath = [PathName FileName];
                    gui.recentPath = PathName;
                    
                    for act_file = 1:no_files_act
                        data{rows+act_file-1,gui.pathColumn} = PathName{act_file};
                    end
                    set(src,'Data',data)
                    colWidth = max(cellfun(@(x) numel(x),data(:,gui.pathColumn)));
                    oldColWidth = get(src,'ColumnWidth');
                    oldColWidth{gui.pathColumn} = 7*colWidth;
                    set(src,'ColumnWidth',oldColWidth);
                    data_new = cell(size(data)+[1 0]);
                    data_new(1:size(data,1),1:size(data,2)) = data;
                    data_new(end,1:size(gui.DefaultData,2)) = gui.DefaultData;
                    data = data_new;
                    
                    %             data = vertcat(data,gui.DefaultData);
                    %             set(src,'Data',vertcat(data,gui.DefaultData));
                    if no_rows>1
                        parameters = handleComboMenus(gui.defPreOpts);
                        if cols==gui.groupColumn % Add Group
                            act_group = data{rows-1,gui.groupColumn}+1;
                            gui.expTree = vertcat(gui.expTree,NaN(1,size(gui.expTree,2))); % Add new row to gui.expTree
                            gui.expTree(end,1)=1; % First trial
                            data{rows,cols} = size(gui.expTree,1);
                            data{rows,gui.subsetColumn} = 1;
                            data{rows,gui.trialColumn} = 1;
                            set(src,'Data',data);
                            for act_file = 1:no_files_act
                                gui.trajectoryList{act_group}{1}{act_file} = PathName{act_file};
                                data{rows+act_file-1,gui.groupColumn} = act_group;
                                data{rows+act_file-1,gui.subsetColumn} = 1;
                                data{rows+act_file-1,gui.trialColumn} = act_file;
                                gui.trajectoryList{act_group}{1}{act_file} = PathName{act_file};
                                gui.loadOptsNoCombos{act_group}{1}{act_file} = parameters;
                                gui.loadOpts{act_group}{1}{act_file} = gui.defPreOpts;
                                
                            end
                            set(src,'Data',data)
                        elseif cols==gui.subsetColumn % Add Subset
                            act_group = data{rows-1,gui.groupColumn};
                            newSubset = find(isnan(gui.expTree(act_group,:)),1,'first');
                            if isempty(newSubset) % Add new columns to gui.expTree
                                gui.expTree=cat(2,gui.expTree,NaN(size(gui.expTree,1),1));
                                newSubset = find(isnan(gui.expTree(act_group,:)),1,'first');
                            end
                            gui.expTree(act_group,newSubset) = 1;
                            
                            for act_file = 1:no_files_act
                                gui.trajectoryList{act_group}{newSubset}{act_file} = PathName{act_file};
                                data{rows+act_file-1,gui.groupColumn} = act_group;
                                data{rows+act_file-1,gui.subsetColumn} = newSubset;
                                data{rows+act_file-1,gui.trialColumn} = act_file;
                                gui.trajectoryList{act_group}{newSubset}{act_file} = PathName{act_file};
                                gui.loadOptsNoCombos{act_group}{newSubset}{act_file} = parameters;
                                gui.loadOpts{act_group}{newSubset}{act_file} = gui.defPreOpts;
                            end
                            set(src,'Data',data)
                        elseif cols==gui.trialColumn % Add trial
                            
                            act_group = data{rows-1,gui.groupColumn};
                            act_subset = data{rows-1,gui.subsetColumn};
                            act_trial = data{rows-1,gui.trialColumn}+1;%nansum([gui.expTree(act_group,act_subset),1]);
                            gui.expTree(act_group,act_subset) = gui.expTree(act_group,act_subset)+1;
                            %                     data{rows,gui.groupColumn} = act_group;
                            %                     data{rows,gui.subsetColumn} = act_subset;
                            %                     data{rows,gui.trialColumn} = act_trial;
                            for act_file = 1:no_files_act
                                data(rows+act_file-1,1:4) = {act_group act_subset act_trial+act_file-1 PathName{act_file}};
                                gui.trajectoryList{act_group}{act_subset}{act_trial+act_file-1} = PathName{act_file};
                                gui.loadOptsNoCombos{act_group}{act_subset}{act_trial+act_file-1} = parameters;
                                gui.loadOpts{act_group}{act_subset}{act_trial+act_file-1} = gui.defPreOpts;
                            end
                            set(src,'Data',data)
                        end
                        
                        
                        
                        
                        
                    else
                        gui.expTree = 1;
                        data = cell(no_files_act,size(gui.DefaultData,2));
                        parameters = handleComboMenus(gui.defPreOpts);
                        
                        for act_file = 1:no_files_act
                            data(act_file,1:4) = {1 1 act_file PathName{act_file}};
                            gui.trajectoryList{1}{1}{act_file} = PathName{act_file};
                            gui.loadOpts{1}{1}{act_file} = gui.defPreOpts;
                            gui.loadOptsNoCombos{1}{1}{act_file} = parameters;
                            
                        end
                        
                        data = vertcat(data,gui.DefaultData);
                        set(src,'Data',data)
                        
                        
                    end
                    %             act_group = data{rows,gui.groupColumn};
                    %             act_subset = data{rows,gui.subsetColumn};
                    %             act_trial = data{rows,gui.trialColumn};
                    for act_file = 1:no_files_act
                        data{rows+act_file-1,gui.preOptColumn} = 'Set';
                        data{rows+act_file-1,gui.removeColumn} = 'Remove';
                    end
                    set(src,'Data',data)
                    
                elseif ~isempty(PathName) && ~isempty(IndexList) %% From folder tree reconstruction
                    parameters = handleComboMenus(gui.defPreOpts);
                    
                    % Delete previous info:
                    data = cell(size(PathName,1),size(gui.DefaultData,2));
                    gui.trajectoryList = {};
                    gui.loadOptsNoCombos = [];
                    gui.loadOpts = [];
                    gui.expTree = NaN(max(IndexList(:,1)),max(IndexList(:,2)));
                    
                    for act_file = 1:size(PathName,1)
                        data(act_file,1:4) = {IndexList(act_file,1) IndexList(act_file,2) IndexList(act_file,3) PathName{act_file}};
                        gui.trajectoryList{IndexList(act_file,1)}{IndexList(act_file,2)}{IndexList(act_file,3)} = PathName{act_file};
                        gui.loadOptsNoCombos{IndexList(act_file,1)}{IndexList(act_file,2)}{IndexList(act_file,3)} = parameters;
                        gui.loadOpts{IndexList(act_file,1)}{IndexList(act_file,2)}{IndexList(act_file,3)} = gui.defPreOpts;
                        gui.expTree(IndexList(act_file,1),IndexList(act_file,2)) = gui.expTree(IndexList(act_file,1),IndexList(act_file,2))+1;
                        data{act_file,gui.preOptColumn} = 'Set';
                        data{act_file,gui.removeColumn} = 'Remove';
                    end
                    data = vertcat(data,gui.DefaultData);
                    set(src,'Data',data)
                end
                set(gui.preprocessButtonH,'Enable','on')
            end
        end
    elseif cols==gui.removeColumn && rows<no_rows
        act_group = data{rows,gui.groupColumn};
        act_subset = data{rows,gui.subsetColumn};
        act_trial = data{rows,gui.trialColumn};
        
        dat = get(src,'Data');
        act_idces = [dat{rows, [gui.groupColumn gui.subsetColumn, gui.trialColumn]}];
        
        dat(rows,:) = [];
        if ~isempty(gui.trajectoryList{act_group}{act_subset})
            gui.trajectoryList{act_group}{act_subset}(act_trial) =[];
            gui.loadOptsNoCombos{act_group}{act_subset}(act_trial) =[];
            gui.loadOpts{act_group}{act_subset}(act_trial) =[];
        end
        
        idces = cell2mat(dat(1:no_rows-2, [gui.groupColumn gui.subsetColumn, gui.trialColumn]));
        % Check if there are trials > then act trial, and if so, adjust count
        if ~isempty(idces)
            same_gr_and_sb = find(idces(:,gui.groupColumn) == act_idces(gui.groupColumn) & ...
                idces(:,gui.subsetColumn) == act_idces(gui.subsetColumn));
            for k=1:numel(same_gr_and_sb)
                if idces(same_gr_and_sb(k),gui.trialColumn)>act_idces(gui.trialColumn)
                    idces(same_gr_and_sb(k),gui.trialColumn) = idces(same_gr_and_sb(k),gui.trialColumn) - 1;
                end
                gui.expTree(act_group,act_subset) = gui.expTree(act_group,act_subset) - 1;
                
            end
            % Check if there are subsets > then act subset and act_subset will be empty after deleting, and if so, adjust count
            same_gr = find(idces(:,gui.groupColumn) == act_idces(gui.groupColumn));% & idces(:,gui.subsetColumn) == act_idces(gui.subsetColumn));
            same_sb = find(idces(:,gui.groupColumn) == act_idces(gui.groupColumn) & idces(:,gui.subsetColumn) == act_idces(gui.subsetColumn));
            
            for k=1:numel(same_gr)
                if idces(same_gr(k),gui.subsetColumn)>act_idces(gui.subsetColumn) && isempty(same_sb)
                    idces(same_gr(k),gui.subsetColumn) = idces(same_gr(k),gui.subsetColumn) - 1;
                    
                end
                
            end
            
            %         same_gr = find(idces(:,gui.groupColumn) == act_idces(gui.groupColumn) & idces(:,gui.subsetColumn) == act_idces(gui.subsetColumn));
            if isempty(same_sb)
                gui.expTree(act_group,act_subset) = NaN;
                for k2 = 1:size(gui.expTree(act_group,act_subset+1:end),2)
                    gui.expTree(act_group,act_subset+k2-1) = gui.expTree(act_group,act_subset+k2);
                end
                gui.expTree(act_group,end) = NaN;
            end
            
            % Check if there are groups > then act group and act_group will be empty after deleting, and if so, adjust count
            same_gr = find(idces(:,gui.groupColumn) == act_idces(gui.groupColumn));
            for k=1:size(idces,1)
                if idces(k,gui.groupColumn)>act_idces(gui.groupColumn) && isempty(same_gr)
                    idces(k,gui.groupColumn) = idces(k,gui.groupColumn) - 1;
                end
                
            end
            if isempty(same_gr)
                gui.expTree(act_idces(gui.groupColumn),:) = [];
            end
            
            
            
            dat(1:no_rows-2, [gui.groupColumn gui.subsetColumn, gui.trialColumn]) = num2cell(idces);
        end
        set(src,'Data',dat);
        data = dat;
    elseif cols==gui.preOptColumn && rows<no_rows % Set individual options
        act_group = data{rows,gui.groupColumn};
        act_subset = data{rows,gui.subsetColumn};
        act_trial = data{rows,gui.trialColumn};
        if isstruct(gui.loadOpts)
            %                 [hPropsPane,parameters] = propertiesGUI([],gui.loadOpts);
%             opts = RestoreComboOptions(gui.loadOpts,gui.defPreOpts);
            parameters = idSocialUI_propertyTable([],gui.loadOpts,true);
        elseif iscell(gui.loadOpts)
            %                 [hPropsPane,parameters] = propertiesGUI([],gui.loadOpts{act_group}{act_subset}{act_trial});
%             opts = RestoreComboOptions(gui.loadOpts{act_group}{act_subset}{act_trial},gui.defPreOpts);
            parameters = idSocialUI_propertyTable([],gui.loadOpts{act_group}{act_subset}{act_trial},true);
            
            
        end
        if isempty(parameters)
            parameters = gui.defPreOpts;
        end
        gui.loadOpts{act_group}{act_subset}{act_trial} = parameters;
        
        parameters = handleComboMenus(parameters);
        gui.loadOptsNoCombos{act_group}{act_subset}{act_trial} = parameters;
        
        
    elseif cols==gui.preOptColumn && rows==no_rows % Set overall options
        if isfield(gui,'loadOpts') && isstruct(gui.loadOpts)
            %                 [hPropsPane,parameters] = propertiesGUI([],gui.loadOpts);
            parameters = idSocialUI_propertyTable([],gui.loadOpts,true);
            
            
        else
            %                 [hPropsPane,parameters] = propertiesGUI([],gui.defPreOpts);
            parameters = idSocialUI_propertyTable([],gui.defPreOpts,true);
            
        end
        if ~isempty(parameters)
            %         if ~isfield(gui,'defPreOpts')
            gui.loadOpts = idSocial_recursiveOptions2OptionsCell(gui.loadOpts,parameters);
            % Handle combo menus
            parameters = handleComboMenus(parameters);
            
            gui.loadOptsNoCombos = idSocial_recursiveOptions2OptionsCell(gui.loadOpts,parameters);
        end
        %         end
        %         gui.loadOpts = parameters;
        %     elseif cols==gui.groupColumn && rows == no_rows+1
        %         loaddir(src,evnt,fh)
    end
    
    if ishandle(fh)
        %         set(src,'ColumnWidth',{75 75 75 300 75 75});
    end
    
end

gui.trajectoryData = data;

if ishandle(fh)
    guidata(fh,gui);
end
end


function [pnl, filelist] = openFolderFileChooser(object, eventdata,fh)
gui=guidata(fh);
% dummyFig = figure('Visible','off'); % Dummy figure to use uiwait
filelist = {};
%         idx2coord(object,eventdata);

if isfield(gui,'contextMenu') && ishandle(gui.contextMenu) && IsGraphicHandle(gui.contextMenu)
    delete(gui.contextMenu)
end
C = get(0,'PointerLocation');%get (gca, 'CurrentPoint');

% With help from
% http://stackoverflow.com/questions/2769430/matlab-convert-global-coordinates-to-figure-coordinates/2770174
% Get figure position
oldUnits = get(get(object,'Parent'),'Units');
set(get(object,'Parent'),'Units','Pixels');
figpos = get(get(object,'Parent'),'Position');
set(get(object,'Parent'),'Units',oldUnits);

oldUnits = get(object,'Units');
set(object,'Units','Pixels');
tablepos = get(object,'Position');
set(object,'Units',oldUnits);

%# Compute an offset and scaling for coords:
offset = figpos(1:2);%+tablepos(1:2);


%# Apply the offsets and scaling:
C = (C-offset);

pnl = uipanel(get(object,'Parent'),'Units','pixels',...
    'Position', [C(1,1)+30 C(1,2)-70 120 120]);
btn(1) = uicontrol(pnl,'Style', 'pushbutton', 'String', 'Open File',...
    'Units','normalized',...
    'FontSize',gui.table_fontsize,...
    'Position', [0 3/4 1 1/4],...
    'Callback', @(src,evnt,x) chooseFile(src,evnt,pnl));
btn(2) = uicontrol(pnl,'Style', 'pushbutton', 'String', 'Open Folder',...
    'Units','normalized',...
    'FontSize',gui.table_fontsize,...
    'Position', [0 2/4 1 1/4],...
    'Callback', @(src,evnt,x) chooseFolder(src,evnt,pnl));
btn(3) = uicontrol(pnl,'Style', 'pushbutton', 'String', 'Open Tree',...
    'Units','normalized',...
    'FontSize',gui.table_fontsize,...
    'Position', [0 1/4 1 1/4],...
    'Callback', @(src,evnt,x) chooseFolderTree(src,evnt,fh));
btn(4) = uicontrol(pnl,'Style', 'pushbutton', 'String', 'Cancel',...
    'Units','normalized',...
    'FontSize',gui.table_fontsize,...
    'Position', [0 0 1 1/4],...
    'Callback', @(src,evnt,x) closeChooseFolder(src,evnt,pnl));
gui.contextMenu = pnl;
guidata(fh,gui);
disp(filelist)

if ~isempty(filelist)
    data = get(object,'Data');
    data{1,1} = filelist{1};
    set(object,'Data',data);
end

    function chooseFolderTree(src,evnt,fh)
        delete(pnl)
        gui=guidata(fh);
        if isfield(gui,'trajectoryTable') && ~isempty(get(gui.trajectoryTable,'Data')) && ...
                ~isequal(get(gui.trajectoryTable,'Data'),gui.DefaultData)
            wh=warndlg('Previous trajectories will be removed from the list!');
            uiwait(wh);
        end
        idSocialUI_chooseFolderTree(fh);
        uiwait(gcf)
    end
    function chooseFolder(src,evnt,pnl)
        gui=guidata(fh);
        delete(pnl)
        if isfield(gui,'recentPath')&& ~isempty(gui.recentPath)
            PathName = uigetdir(gui.recentPath{1});
        else
            PathName = uigetdir;
        end
        
        if (isnumeric(PathName) && PathName==0)
            gui.openFolderFileChooserList = 'canceled';
        else
            filelist = idSocial_recursiveDir(PathName,'*.mat');
            gui.openFolderFileChooserList = filelist;
        end
        guidata(fh,gui);
    end

    function chooseFile(src,evnt,pnl)
        gui=guidata(fh);
        delete(pnl)
        if isfield(gui,'recentPath')&& ~isempty(gui.recentPath)
            [FileName,PathName] = uigetfile(gui.recentPath{1});
        else
            [FileName,PathName] = uigetfile;
        end
        
        if (isnumeric(PathName) && PathName==0) || (isnumeric(PathName) && FileName==0)
            gui.openFolderFileChooserList = 'canceled';
        else
            filelist = {[PathName FileName]};
            gui.openFolderFileChooserList = filelist;
        end
        guidata(fh,gui);
    end
    function closeChooseFolder(src,evnt,pnl)
        gui=guidata(fh);
        
        delete(pnl)
        gui.contextMenu = NaN;
        gui.openFolderFileChooserList = 'canceled';
        guidata(fh,gui);
        
    end
end
% function loaddir(src,evnt,fh)
% % Called when user activates popup menu
%
% gui=guidata(fh); %
% dirname=uigetdir;%(gui.path);
% fname=rdir([dirname '\**\trajectories.mat']);
% fname=vertcat({fname.name})';
% % if ~isempty(fname)
% %     if strcmp(get(lbh,'String'),'Please select file/directory')
% %         set(lbh,'String',fname);
% %     else
% %         set(lbh,'String',vertcat(get(lbh,'String'),fname));
% %     end
% %     gui.alldata.trajectory{get(lbh,'UserData')}=get(lbh,'String');
% % end
% guidata(fh,gui)
%
% end

function preProcess(src,evnt,fh)
gui=guidata(fh);
% preH = figure('Units','normalized','Position',[.3 .3 .4 .5]);
% prePrDefOpt = idSocial_loadData;
% data(:,1) = fieldnames(prePrDefOpt);
% data(:,2) = cellfun(@(x) num2str(x),struct2cell(prePrDefOpt),'UniformOutput',false);
% propertiesGUI(prePrDefOpt)
% t = uitable('Data', data,...
%             'ColumnName', {'Options' 'Value'},...
%             'ColumnFormat', {'char' 'char'},...
%             'ColumnEditable', [false true],...
%             'ColumnWidth',{475 275},...
%             'RowName',[],...
%             'FontSize',14 );
%         % Set width and height
% t.Units='normalized';
% t.Position(3) = .95;%fpos(3)*.9;%t.Extent(3);
% t.Position(4) = .7;
% t.Position(2) = .2;

% startPreProButtonH=uicontrol('Style','pushbutton' ,'String','Start', ...
%     'Units', 'normalized','Position', [.4 .05 .2 .1],'FontSize',gui.fontsize_button,...
%     'Callback',@(src,evnt) preProcess(src,evnt,fh));
% gui.prePrDefOpt = prePrDefOpt;
cont = true;

if ~gui.first_preprocess
    choice = questdlg('Previous, not-safed analyses will be lost. Continue?', ...
        'Continue preprocessing?', ...
        'Ok','Cancel','Cancel');
    switch choice
        case 'Ok'
            cont=true;
        case 'Cancel'
            cont=false;
    end
end



if cont==true
    set(gui.DataSetTable,'Columnname',{'New Set'});
    gui.DataSetTableColumnNames = {'New Set'};
    if isfield(gui,'DataSetTable')
        
        set(gui.DataSetTable,'Data',gui.DataSetTableDefaultData)
        
    end
    gui.input_data = [];
    fn = fieldnames(gui);
    for k=1:size(fn,1)
        if isfield(gui.(fn{k}),'DataSetTableData')
            
            gui = rmfield(gui,fn{k});
            if isfield(gui,'plot_mode') && isfield(gui.plot_mode,fn{k})
                gui.plot_mode = rmfield(gui.plot_mode,fn{k});
            end
        end
    end
    
    
    set(gui.funcH,'Enable','off')
    set(gui.secH,'Enable','off')
    set(gui.funcOptsH,'Enable','off')
    set(gui.funcAdvOptsH,'Enable','off')
    set(gui.funcHelpH,'Enable','off')
    set(gui.DataSetTable,'Enable','off')
    set(gui.preprocessButtonH,'Enable','off')
    set(gui.export2workspaceButtonH,'Enable','off')
    set(gui.displayResultsButtonH,'Enable','off')
%     set(gui.newFigureButtonH,'Enable','off')
    set(gui.trajectoryTable,'Enable','off')
    set(gui.plotResultsButtonH,'Enable','off')
    pause(.1)
    
    try
        gui.input_data=idSocial_loadData(gui.trajectoryList, [], gui.loadOptsNoCombos);
        
        set(gui.funcH,'Enable','on')
        set(gui.secH,'Enable','on')
        set(gui.funcOptsH,'Enable','on')
        set(gui.funcAdvOptsH,'Enable','on')
        set(gui.funcHelpH,'Enable','on')
        set(gui.DataSetTable,'Enable','off')
        set([gui.FigureTableButtonsH{:}],'Enable','off')
        set(gui.export2workspaceButtonH,'Enable','off')
        set(gui.displayResultsButtonH,'Enable','off')
        set(gui.plotResultsButtonH,'Enable','off')
    catch ex
        msg = ex.message;
        if strfind(ex.message,'idSocial_loadData: Could not create')
            msg = vertcat(msg,{'Please try a different folder location.'});
        end
        errordlg(msg,'idSocialUI error')
    end
    set(gui.trajectoryTable,'Enable','on')
    set(gui.preprocessButtonH,'Enable','on')
    
    
    tData = gui.trajectoryData;
    colName = get(gui.trajectoryTable,'Columnname');
    colName{gui.framerateColumn} = '<html><font size=+1>Framerate';
    for row = 1:size(tData,1)-1
        tData{row,gui.framerateColumn}=gui.input_data.info.framerate(tData{row,gui.groupColumn},tData{row,gui.subsetColumn},tData{row,gui.trialColumn});
    end
    colName{gui.indivColumn} = '<html><font size=+1># Indivs.';
    for row = 1:size(tData,1)-1
        tData{row,gui.indivColumn}=gui.input_data.info.no_focals(tData{row,gui.groupColumn},tData{row,gui.subsetColumn},tData{row,gui.trialColumn});
    end
    colName{gui.bodylengthColumn} = '<html><font size=+1>Bodylength [PXL]';
    for row = 1:size(tData,1)-1
        tData{row,gui.bodylengthColumn}=gui.input_data.info.bodylength_in_pixels(tData{row,gui.groupColumn},tData{row,gui.subsetColumn},tData{row,gui.trialColumn});
    end
    set(gui.trajectoryTable,'Data',tData,'Columnname',colName)
    set(gui.trajectoryTable,'ColumnWidth',{75 75 75 300 75 75 75 85 95});
    
    
    
    
    gui.first_preprocess = false;
    guidata(fh,gui);
end
end

function sectionSelect(src,evnt,fh)
gui=guidata(fh);
% if isfield(gui,'selectedSection')
%     gui.selectedSectionPrev = gui.selectedSection;
% end
% if isfield(gui,'selectedFunc')
%     gui.selectedFuncPrev = gui.selectedFunc;
% end
if isempty(src)
    gui.selectedSection = 1;
else
    gui.selectedSection = get(src,'Value');
end


set(gui.funcH,'Value',1);
set(gui.funcH,'String',gui.funcCell{gui.selectedSection,2});
set(gui.funcAdvOptsH,'Enable','off');
set(gui.funcHelpH,'Enable','off')
set(gui.funcOptsH,'Enable','off');
set(gui.DataSetTable,'Enable','off')
guidata(fh,gui);
end

function functionSelect(src,evnt,fh)
gui=guidata(fh);

if isempty(src)
    gui.selectedFunc = 1;
else
    gui.selectedFunc = get(src,'Value');
end

if isfield(gui,gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)) && isfield(gui.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]),'DataSetTable')
    set(gui.DataSetTable,'Data',gui.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).DataSetTable)
end

% Get default options and split function specific from general ones
defFuncOps = eval([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2)]);

if isfield(gui.input_data,gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2))
    fnames = fieldnames(gui.input_data.(gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)));
    setnames = cellfun(@(x) ~isempty(x), strfind(fnames,'Set'));
    snames = cellfun(@(x) ['Set ' x(4:end)],fnames(setnames),'UniformOutput',false);
    
    if any(setnames)
        set(gui.DataSetTable,'Columnname',snames);
    else
        snames = {'New Set'};
        set(gui.DataSetTable,'Columnname',snames);
    end
    gui.DataSetTableColumnNames = snames;
else
    set(gui.DataSetTable,'Columnname',{'New Set'});
    gui.DataSetTableColumnNames = {'New Set'};
end
% if size(defFuncOps,2)>1
%     defFuncOps = defFuncOps(1);
% end
% Check if random controls are available:
[~ , rdata]= idSocial_recursiveGetOptionsFromOptionsCell(gui.input_data.options,'random_data');

defFuncOpsNames = fieldnames(defFuncOps);
[def_optionsAll, def_optionsTypes] = idSocial_auxiliaries_createDefOptions;
for k=1:size(defFuncOpsNames,1)
    if isfield(def_optionsAll,defFuncOpsNames{k}) && ~strcmp(defFuncOpsNames{k},'timeintervals_in_min') && ...
            (~strcmp(defFuncOpsNames{k},'random_data') || ~rdata)
        defFuncOps = rmfield(defFuncOps,defFuncOpsNames{k});
    end
end
gui.act_method = defFuncOps(1).act_method;
defFuncOps = rmfield(defFuncOps,'act_method');
gui.funcOpts.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=defFuncOps;
parameters = handleComboMenus(defFuncOps);
gui.funcOptsNoCombos.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=parameters;
gui.funcOptsTypes.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=def_optionsTypes;

% Get default advanced options and split function specific from general ones
defFuncOps = eval([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2)]);
% if size(defFuncOps,2)>1
%     defFuncOps = defFuncOps(1);
% end
defFuncOpsNames = fieldnames(defFuncOps);
[def_optionsAll, def_optionsTypes] = idSocial_auxiliaries_createDefOptions;
for k=1:size(defFuncOpsNames,1)
    if ~isfield(def_optionsAll,defFuncOpsNames{k}) || strcmp(defFuncOpsNames{k},'timeintervals_in_min') || ...
            (strcmp(defFuncOpsNames{k},'random_data') && rdata)
        defFuncOps = rmfield(defFuncOps,defFuncOpsNames{k});
    end
end
gui.funcAdvOpts.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=defFuncOps;
parameters = handleComboMenus(defFuncOps);
gui.funcAdvOptsNoCombos.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=parameters;


% if isfield(gui,gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)) && ...
%         isfield(gui.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]),'DataSetTableData')
%     set(gui.DataSetTable,'Data',gui.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).DataSetTableData);
%     
% else
%     set(gui.DataSetTable,'Data',gui.DataSetTableDefaultData);
%     
% end
guidata(fh,gui);
updateDataSetTable(fh);
gui=guidata(fh);


set(gui.funcAdvOptsH,'Enable','on');
set(gui.funcHelpH,'Enable','on')
set(gui.funcOptsH,'Enable','on');
set(gui.DataSetTable,'Enable','on')
guidata(fh,gui);
% enStatus = get([gui.FigureTableButtonsH{:}],'Enable');

guidata(fh,gui);
updateFigTable(fh);
gui=guidata(fh);

dt=get(gui.DataSetTable,'Data');
parentStraight = [gui.FigureTableButtonsH{:}];
if ~isempty(parentStraight) && ...
        all(IsGraphicHandle(parentStraight))
    if size(dt,1)>2 && size(dt,2)>1
        set([gui.FigureTableButtonsH{:}],'Enable','on');
    else
        set([gui.FigureTableButtonsH{:}],'Enable','off');
    end
else
    set([gui.FigureTableButtonsH{:}],'Enable','off');
end
guidata(fh,gui);

end



function functionOptions(src,evnt,fh)
gui=guidata(fh);

try
    defFuncOps = gui.funcOpts.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);
catch ex
    warndlg(['Something went wrong: ' ex.message])
    
end
def_optionsTypes = ...
    gui.funcOptsTypes.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);

% [hPropsPane,parameters] = propertiesGUI([],defFuncOps);
parameters = idSocialUI_propertyTable([],defFuncOps,true);

if ~isempty(parameters)
    defFuncOpsNames = fieldnames(defFuncOps);
    for k=1:size(defFuncOpsNames,1)
        if ~isempty(parameters(1).(defFuncOpsNames{k})) && isfield(def_optionsTypes,defFuncOpsNames{k})
            switch def_optionsTypes.(defFuncOpsNames{k}){1}
                case 'numeric'
                    if ischar(parameters(1).(defFuncOpsNames{k}))
                        parameters(1).(defFuncOpsNames{k}) = str2num(parameters(1).(defFuncOpsNames{k}));
                    end
            end
        end
    end
    
    gui.funcOpts.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=parameters;
    parameters = handleComboMenus(parameters);
    gui.funcOptsNoCombos.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=parameters;
end
guidata(fh,gui);
end

function functionAdvOptions(src,evnt,fh)
gui=guidata(fh);


defFuncOps = gui.funcAdvOpts.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);
def_optionsTypes = ...
    gui.funcOptsTypes.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);

% [hPropsPane,parameters] = propertiesGUI([],defFuncOps);
parameters = idSocialUI_propertyTable([],defFuncOps,true);

if ~isempty(parameters)
    defFuncOpsNames = fieldnames(defFuncOps);
    for k=1:size(defFuncOpsNames,1)
        if ~isempty(parameters(1).(defFuncOpsNames{k})) && isfield(def_optionsTypes,defFuncOpsNames{k})
            switch def_optionsTypes.(defFuncOpsNames{k}){1}
                case 'numeric'
                    if ischar(parameters(1).(defFuncOpsNames{k}))
                        parameters(1).(defFuncOpsNames{k}) = str2num(parameters(1).(defFuncOpsNames{k}));
                    end
            end
        end
    end
    
    gui.funcAdvOpts.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=parameters;
    parameters = handleComboMenus(parameters);
    gui.funcAdvOptsNoCombos.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)])=parameters;
end
guidata(fh,gui);
end


function funcHelp(src,evnt,fh)
gui=guidata(fh);
funcstring = [gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}];
if ~isempty(funcstring)
    idSocial_auxiliaries_extractHelp(funcstring);
end

end

function setPlotMode(src,evnt,fh,actSet)
gui=guidata(fh);



panel_fontsize = 12;
gui.table_fontsize_html = 4;

% opts = gui.funcOpts.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);
% gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]) = ...
%     eval([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2) '(opts)']);


pmfh = figure('Units','normalized','Position',[.45 .3 .4 .5], 'ToolBar', 'none','MenuBar', 'none','Name','Set plot parameters','NumberTitle','off','WindowStyle','modal');


hsp0 = uipanel('Parent',pmfh,'Title','Select data','FontSize',panel_fontsize);


% 'CellSelectionCallback',@(src,evnt) cellSelectPlotMode(src,evnt,fh)
% Set width and height
set(hsp0,'Units','normalized');
pos = get(hsp0,'Position');
set(hsp0,'Position',[pos(1)+.01 pos(2)+.05 .9 .9]);

set(hsp0,'Units','centimeter');
pos = get(hsp0,'Position');
set(hsp0,'Position',[pos(1) pos(2) 8 pos(4)]);
set(hsp0,'Units','normalized');
pos = get(hsp0,'Position');

% Create the uitable: Choose data (plot_mode)
gui.plot_modeTable = uitable('Parent',hsp0,'Data', [],...
    'ColumnName', {'<html><font size="1">Group<font size="1"></html>','Subset','Trial','Time','Focal','Neighbor'},...
    'ColumnEditable', [true false false false false false false],...
    'ColumnWidth',{35 50 50 50 50 50 50 },...
    'ColumnFormat', {'logical' 'numeric', 'numeric', 'numeric',  'numeric', 'numeric','numeric'},...
    'RowName',[],...
    'FontSize',gui.table_fontsize,'Units','normalized','Position',[0.01 0.01 .98 .98]);
set(gui.plot_modeTable,'ColumnName', {['<html><font size="' num2str(gui.table_fontsize_html) '">Select<font size="1"></html>'], ...
    ['<html><font size="' num2str(gui.table_fontsize_html) '">Group<font size="1"></html>'], ...
    ['<html><font size="' num2str(gui.table_fontsize_html) '">Subset<font size="1"></html>'], ...
    ['<html><font size="' num2str(gui.table_fontsize_html) '">Trial<font size="1"></html>'], ...
    ['<html><font size="' num2str(gui.table_fontsize_html) '">Time<font size="1"></html>'], ...
    ['<html><font size="' num2str(gui.table_fontsize_html) '">Focal<font size="1"></html>'], ...
    ['<html><font size="' num2str(gui.table_fontsize_html) '">Neighbor<font size="1"></html>']});

if gui.tooltips_on
    s = sprintf('Select row(s) to include corresponding indices in analysis.\nPress Ctrl to select various rows.\nClick on table header to sort rows.');
    set(gui.plot_modeTable,'TooltipString',s);
end
% Sortable (http://undocumentedmatlab.com/blog/uitable-sorting)


hJScroll = findjobj(gui.plot_modeTable); % findjobj is from file exchange
hJTable = hJScroll.getViewport.getView; % get the table component within the scroll object
% hJTable.setNonContiguousCellSelection(false);
% hJTable.setColumnSelectionAllowed(false);
% hJTable.setRowSelectionAllowed(true);

hJTable = handle(hJTable, 'CallbackProperties');

jscrollpane = findjobj(gui.plot_modeTable);
jtable = jscrollpane.getViewport.getView;
jheader= jtable.getTableHeader(); % Here, you got JTableHeader object.
h=handle(jheader, 'callbackproperties');


set(gui.plot_modeTable,'ColumnWidth',{35 70 70 70 70 70 70});


ti = gui.funcOptsNoCombos.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).timeintervals_in_min;
if ~isempty(ti)
    fr = gui.input_data.info.framerate;
    no_frames=   gui.input_data.info.no_frames;
    no_frames_part_array=floor(ti.*fr*60)-1;
    no_parts=2*floor(no_frames./no_frames_part_array)+1;
    tsteps = 1:max(no_parts(:));
else
    tsteps = 2;
end

% Def. plot_mode:
opts = gui.funcOptsNoCombos.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);
gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]) = ...
    eval([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2) '(opts)']);

focaldim = true;
neighbordim=true;
if isfield(gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]),'extraDims') && ...
        ~isempty(gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).extraDims)
    extrDim = gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).extraDims;
    if  isempty(extrDim{1})
        focaldim=false; % Set 1, if indices are the same in all combinations they will not appear in table.
    end
    if size(extrDim,2)>1 && isempty(extrDim{2})
        neighbordim=false;
    end
end

idxcombs = NaN(1,6);
for gr=1:gui.input_data.info.no_groups
    for sb=1:gui.input_data.info.no_subsets(gr)
        for tr=1:gui.input_data.info.no_trials(gr,sb)
            for tm=tsteps
                if focaldim
                    no_focals=gui.input_data.info.no_focals(gr,sb,tr);
                else
                    no_focals = 1;
                end
                for fc=1:no_focals
                    if neighbordim
                        no_neighbors=gui.input_data.info.no_neighbors(gr,sb,tr);
                    else
                        no_neighbors=1;
                    end
                    for nb=1:no_neighbors
                        if (focaldim && neighbordim && (fc~=nb || no_focals==1)) || ...
                                ~focaldim || ~neighbordim
                            idxcombs = vertcat(idxcombs,[gr, sb, tr, tm, fc, nb]);
                        end
                    end
                end
            end
        end
    end
end
idxcombs(1,:)=[];
gr = idxcombs(:,1);
sb = idxcombs(:,2);
tr = idxcombs(:,3);
tm = idxcombs(:,4);
fc = idxcombs(:,5);
nb = idxcombs(:,6);


% idstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
% idstring = idstring(10:end);
% load(gui.input_data.(idstring).output)
% idx=find(cellfun(@(x) ~isempty(x),output_new));
% [gr, sb, tr, tm, fc, nb]=ind2sub(size(output_new),idx);

no_idx = size(gr,1);
colname = {'Select' 'Group', 'Subset','Trial', 'Time','Focal','Neighbor'};
t = [false true true true true true true];
if size(gr,1)>1
    if all(gr==gr(1))
        gr=[];
        t(2) = false;
    end
    if all(sb==sb(1))
        sb=[];
        t(3)=false;
    end
    if all(tr==tr(1))
        tr=[];
        t(4)=false;
    end
    if all(fc==fc(1))
        fc=[];
        t(6)=false;
    end
    if all(nb==nb(1))
        nb=[];
        t(7)=false;
    end
    if all(tm==tm(1))
        tm=[];
        t(5) = false;
    end
end
colname = colname(t);
gui.colname = colname;




% Show idx selection
gui.tmstring2idx = [];
if all(tm==2) || isempty(tm)
    %     set(gui.plot_modeTable,'Enable','on','Data',[true(size(idxcombs,1),1) gr, sb, tr, fc, nb], ...
    %         'ColumnName', ['Select' colname],'ColumnFormat', {'logical', 'numeric', 'numeric', 'numeric',  'numeric', 'numeric','numeric'});
    dat = [repmat({true},[size(idxcombs,1),1]) num2cell(gr), num2cell(sb), num2cell(tr), num2cell(fc), num2cell(nb)];
    set(gui.plot_modeTable,'Enable','on','Data',dat, ...
        'ColumnName', ['Select' colname],'ColumnFormat', {'logical', 'numeric', 'numeric', 'numeric',  'numeric', 'numeric','numeric'});
    
else
    ti = gui.funcOptsNoCombos.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).timeintervals_in_min;
    fr = gui.input_data.info.framerate;
    no_frames=   gui.input_data.info.no_frames;
    
    startint = (tm-2)*ti/2;
    endint = tm*ti/2;
    startint(startint<0)=0;
    no_frames_part_array=floor(ti.*fr*60)-1;
    no_parts=2*floor(no_frames./no_frames_part_array)+1;
    max_no_parts=max(max(max(no_parts)));
    endint(endint==max_no_parts) = ...
        (max_no_parts-1)*ti/2;
    tmstring = strcat(num2str(startint),  ' - ', num2str(endint));
    utstr =tmstring;
    utm = tm;
    gui.tmstring2idx = cell(size(utm,1),2);
    gui.tmstring2idx(:,1) = cellstr(utstr);
    gui.tmstring2idx(:,2) = num2cell(utm);
    %         dat = cell(no_idx,sum(t));
    %         for k=1:6
    %             dat(:,cnt)=
    %         end
    dat = [repmat({true},[size(idxcombs,1),1]) num2cell(gr), num2cell(sb), num2cell(tr), cellstr(tmstring), num2cell(fc), num2cell(nb)];
    %         dat(:,1)=num2cell(gr);
    %         dat(:,2)=num2cell(sb);
    %         dat(:,3)=num2cell(tr);
    %         dat(:,4)=cellstr(tmstring);
    % %         dat(:,5)=num2cell(endint);%cellstr(tmstring);
    %         dat(:,5)=num2cell(fc);
    %         dat(:,6)=num2cell(nb);
    
    set(gui.plot_modeTable,'Enable','on','Data',dat, ...
        'ColumnFormat', {'logical' 'numeric', 'numeric', 'numeric',  'numeric', 'numeric','numeric'},...
        'ColumnName', ['Select' colname]);
end
gui.colname = colname;
% Set a matlab function as MouseClickedCallback
set(h, 'MouseClickedCallback', {@onHeaderClick, jtable, fh, gui.plot_modeTable,gui.tmstring2idx});
set(hJTable, 'MouseReleasedCallback', {@(src,evnt) cellSelectPlotMode(src,evnt,fh,gui,actSet)});
set(hJTable, 'MouseClickedCallback', {@(src,evnt) cellDeSelectPlotMode(src,evnt,fh,gui)});

% idces = get(gui.plot_modeTable,'Data');
% idces = idces([idces{:,1}]',2:end);
% idces = cell2mat(idces);
% funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
% colname = gui.colname;
% sbidx=strcmp(colname,'Subset');
% if any(sbidx)
%     colname{sbidx} = 'Day';
% end
%
%
% filterstring = {};
% for k=1:size(idces,1)
%     act_fstring=[];
%     for k2=1:size(colname,2)
%         if strcmp(colname{k2},'Time')
%             act_fstring = [act_fstring colname{k2} gui.tmstring2idx(find(strcmp(idces(k,k2),gui.tmstring2idx(:,1)),1,'first'),2)];
%         else
%             if iscell(idces(k,k2)) % With time present, idces is a cell,else array
%                 act_fstring = [act_fstring colname{k2} idces(k,k2)];
%             else
%                 act_fstring = [act_fstring colname{k2} {idces(k,k2)}];
%             end
%         end
%     end
%     filterstring = vertcat(filterstring,act_fstring);
% end
plot_mode_def = gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);
% gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]) = plot_mode_def;
% gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).data_filter = filterstring;
guidata(fh,gui);

% Select x-axis
hsp = uipanel('Parent',pmfh,'Title','Choose X-Axis','FontSize',panel_fontsize,...
    'Position',[pos(1)+.5 pos(2)+.6 .4 .2]);
set(hsp,'Position',[pos(1)+.5 pos(2)+.6 .19 .15]);

xaxString = plot_mode_def.xaxis;%gui.colname;
gui.chooseXaxisH = uicontrol(hsp,'Units','normalized','Style','popup','String',xaxString);
set(gui.chooseXaxisH,'Position', [.1 .05 .8 .8])
if numel(xaxString)==1; set(gui.chooseXaxisH,'Enable','off'); end

% Select x-axis
hsp1 = uipanel('Parent',pmfh,'Title','Choose Y-Axis','FontSize',panel_fontsize,...
    'Position',[pos(1)+.5 pos(2)+.6 .4 .2]);
set(hsp1,'Position',[pos(1)+.5+.21 pos(2)+.6 .19 .15]);

yaxString = plot_mode_def.yaxis;%gui.colname;
if isempty(yaxString); yaxString = ' '; end
gui.chooseYaxisH = uicontrol(hsp1,'Units','normalized','Style','popup','String',yaxString);
set(gui.chooseYaxisH,'Position', [.1 .05 .8 .8])
if numel(yaxString)==1; set(gui.chooseYaxisH,'Enable','off'); end


if isfield(plot_mode_def,'normalization') && ~isempty(plot_mode_def.normalization)&& ...
        (strcmpi(plot_mode_def.display_mode,'Map') || strcmpi(plot_mode_def.display_mode,'MapPolar') || strcmpi(plot_mode_def.display_mode,'hist'))% for histograms and normalization
    indiv_stat_pos_panel = [pos(1)+.5 pos(2)+.4+.05 .4 .15];
    indiv_stat_pos = [.1 .05 .8 .8];
    overall_stat_pos_panel = [pos(1)+.5 pos(2)+.2+.05+.05 .4 .15];
    overall_stat_pos = [.1 .05 .8 .8];
    norm_pos_panel = [pos(1)+.5 pos(2)+.05+.05+.05 .4 .15];
    norm_pos = [.1 .05 .8 .8];
else %if all(cellfun(@(x) isempty(strfind(x,'Edge')),xaxString))
    indiv_stat_pos_panel = [pos(1)+.5 pos(2)+.4 .4 .15];
    indiv_stat_pos = [.1 .05 .8 .8];
    overall_stat_pos_panel = [pos(1)+.5 pos(2)+.2 .4 .15];
    overall_stat_pos = [.1 .05 .8 .8];
end

% Select indiv statistics
hsp2 = uipanel('Parent',pmfh,'Title','Statistics (Indiv.)','FontSize',panel_fontsize,...
    'Position',[pos(1)+.5 pos(2)+.3 .4 .2]);
set(hsp2,'Position',indiv_stat_pos_panel);


statString = plot_mode_def.statistics;%{'Mean', 'Median', 'Pool'};
if isempty(statString); statString = ' '; end
gui.chooseStatH = uicontrol(hsp2,'Units','normalized','Style','popup','String',statString,'Callback',@chooseStatCallback);
set(gui.chooseStatH,'Position',indiv_stat_pos )
if numel(statString)==1; set(gui.chooseStatH,'Enable','off'); end

% Select overall statistics
hsp3 = uipanel('Parent',pmfh,'Title','Statistics (Whole Set)','FontSize',panel_fontsize,...
    'Position',[pos(1)+.5 pos(2)+.3 .4 .2]);
set(hsp3,'Position',overall_stat_pos_panel);


statString2 = plot_mode_def.data;%{'Mean', 'Median'};
if isempty(statString2); statString2 = ' '; end
gui.chooseStatH2 = uicontrol(hsp3,'Units','normalized','Style','popup','String',statString2);
set(gui.chooseStatH2,'Position',overall_stat_pos)
if numel(statString2)==1; set(gui.chooseStatH2,'Enable','off'); end


if isfield(plot_mode_def,'normalization') && ~isempty(plot_mode_def.normalization)&& ...
        (strcmpi(plot_mode_def.display_mode,'Map') || strcmpi(plot_mode_def.display_mode,'MapPolar') || strcmpi(plot_mode_def.display_mode,'hist'))% for histograms and normalization
    % Select normalization (for histograms)
    hsp4 = uipanel('Parent',pmfh,'Title','Normalization','FontSize',panel_fontsize,...
        'Position',[pos(1)+.5 pos(2)+.3 .4 .2]);
    set(hsp4,'Position',norm_pos_panel);
    
    statString3 = plot_mode_def.normalization;%{'Mean', 'Median'};
    if isempty(statString3); statString3 = ' '; end
    gui.chooseStatH3 = uicontrol(hsp4,'Units','normalized','Style','popup','String',statString3);
    set(gui.chooseStatH3,'Position',norm_pos)
    if numel(statString3)==1; set(gui.chooseStatH3,'Enable','off'); end
end
% if ~isempty(strtrim(xaxString(get(gui.chooseXaxisH,'Value')))); plot_mode.xaxis = xaxString(get(gui.chooseXaxisH,'Value')); end;
% if ~isempty(strtrim(statString(get(gui.chooseStatH,'Value')))); plot_mode.statistics = statString(get(gui.chooseStatH,'Value')); end;
% if ~isempty(strtrim(statString2(get(gui.chooseStatH2,'Value')))); plot_mode.data = statString2(get(gui.chooseStatH2,'Value')); end;

guidata(fh,gui);

% Cancel and Ok-Button
gui.plot_modeTableHCancel = uicontrol(pmfh,'Units','normalized','Style','pushbutton','String','Cancel');
set(gui.plot_modeTableHCancel,'Position', [pos(1)+.52 pos(2) .15 .1])

gui.plot_modeTableHOk = uicontrol(pmfh,'Units','normalized','Style','pushbutton','String','Ok', ...
    'Enable','off');
set(gui.plot_modeTableHOk,'Position', [pos(1)+.73 pos(2) .15 .1])



    function chooseStatCallback(src,evnt,fh)
        if strcmpi(src.String(src.Value),'pool')
            statString2 = plot_mode_def.statistics;
            poolidx = cellfun(@(x) ~isempty(x),strfind(lower(statString2),'pool'));
            statString2(poolidx) = [];
        else
            statString2 = plot_mode_def.data;
            
        end
        set(gui.chooseStatH2,'String',statString2);
    end


set(gui.plot_modeTableHCancel,'Callback',@(src,evnt) plot_modeTableCancel(src,evnt,fh));
set(gui.plot_modeTableHOk,'Callback',@(src,evnt) plot_modeTableOk(src,evnt,gui,fh,actSet));
set(gui.plot_modeTableHOk,'Enable','on')

guidata(fh,gui);
end

% function clickDataSet(src,evnt,jtable, fh, pmTable)
%     cellSelectDataSet(pmTable,[],fh)
% end

function cellSelectDataSet(src,evnt,fh,act_set)

gui=guidata(fh);
act_indices = [];
% try % if evnt is java
act_indices = evnt.Indices;
% catch
% end
% if ~isempty(evnt) && isfield( evnt,'Indices')
%         act_indices = get(evnt,'Indices');
% end

%
% if isempty(act_indices) && ...
%         isfield(gui,'DataSetTableIndices') && ~isempty(gui.DataSetTableIndices)
%     act_indices = gui.DataSetTableIndices;
% elseif ~isempty(act_indices)
%     gui.DataSetTableIndices = act_indices;
% else
%     act_indices = [1 1];
%     gui.DataSetTableIndices = act_indices;
% end

% if isempty(src)
%     src=fh;
% end


if ~isempty(act_indices)
    data = get(src,'Data');
    rows = act_indices(1,1);
    cols = act_indices(1,2);
    set(src,'UserData',act_indices(1,2))
    no_cols = size(get(src,'Data'),2);
    no_rows = size(get(src,'Data'),1);
    
    % Try to find appropriate axis labels
    labels=[];
    funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
    
    last_row = sum(cellfun(@(x) ~isempty(x),data(:,cols)));
    
    act_data = [];
    if isfield(gui.input_data,funcstring(10:end)) && isfield(gui.input_data.(funcstring(10:end)),['Set' num2str(cols)] ) && ...
            isfield(gui.input_data.(funcstring(10:end)),'results')
        act_data = gui.input_data.(funcstring(10:end)).(['Set' num2str(cols)]);
        act_results = gui.input_data.(funcstring(10:end)).results;
        %         disp(['Set' num2str(cols)])
        %         disp(gui.input_data.(funcstring(10:end)))
        gui.selectedSet = cols;
    end
    
    
    try
        % Add new set
        if cols<=no_cols
            gui.selectedSet = cols;
        end
        if numel(cols)==1 && numel(rows)==1
            if no_cols==cols && rows==1
                %                 gui.selectedSet = cols;
                setPlotMode(src,evnt,fh,gui.selectedSet)
            end
            % Indices
            if no_cols>cols && rows==1
                funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                %                 dat = gui.plot_mode.(funcstring(10:end)).(['Set' num2str(cols)]).data_filter;
                dat = act_data.plot_mode.data_filter;
                %                 Table_DisplayData(src,evnt,fh,dat,cols)
                gui.displayResultsData = dat;
                gui.plotResultsData = 2;
            end
            % Results
            if no_cols>cols && rows==2
                %                 funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                %                 dat = gui.(funcstring(10:end)).results;
                dat = act_results;
                Table_DisplayData(src,evnt,fh,dat{3,cols+1},cols)
                gui.displayResultsData = dat{3,cols+1};
                gui.plotResultsData = 3;%
                gui.displayResultsIndivData = 4;%dat{4,cols+1};
            end
            % Indiv. Results
            if no_cols>cols && rows==3
                %                 funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                %                 dat = gui.(funcstring(10:end)).results{4,cols+1};
                dat = act_results{4,cols+1};
                if isnumeric(dat) && ndims(dat)<3
                    dat = dat(~all(isnan(dat),2),:);
                    %                     Table_DisplayData(src,evnt,fh,dat,cols)
                    gui.displayResultsData = dat;
                    gui.plotResultsData = 4;
                elseif ndims(dat)>=3
                    gui.displayResultsData = {'Cannot display data'};
                    gui.plotResultsData = -1;
                end
            end
            % Results (Distribution)
            if no_cols>cols && rows==4
                %                 funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                %                 dat = gui.(funcstring(10:end)).results;
                dat = act_results;
                %                 Table_DisplayData(src,evnt,fh,dat{5,cols+1},cols)
                gui.displayResultsData = dat{5,cols+1};
                gui.plotResultsData = 5;
            end
            % Indiv. Results (Distribution)
            if no_cols>cols && rows==5
                funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                %                 dat = gui.(funcstring(10:end)).results{6,cols+1};
                dat = act_results{6,cols+1};
                if isnumeric(dat)
                    dat = dat(~all(isnan(dat),2),:);
                    %                     Table_DisplayData(src,evnt,fh,dat,cols)
                    gui.displayResultsData = dat;
                    gui.plotResultsData = 6;
                end
            end
            if no_cols>cols && rows==6
                %                 dat = gui.(funcstring(10:end)).results{7,cols+1};
                dat = act_results{7,cols+1};
                gui.displayResultsData = dat;
                gui.plotResultsData = 7;
            end
            if no_cols>cols && rows==7 % Ticks
                funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                %                 dat = gui.(funcstring(10:end)).results{8,cols+1};
                dat = act_results{8,cols+1};
                %                 if isnumeric(dat)
                %                     dat = dat(~all(isnan(dat),2),:);
                %                     Table_DisplayData(src,evnt,fh,dat,cols)
                gui.displayResultsData = dat;
                gui.plotResultsData = 8;
                %                 end
            end
            if no_cols>cols && rows==8 % TickLabels
                %                 dat = gui.(funcstring(10:end)).results{9,cols+1};
                dat = act_results{9,cols+1};
                gui.displayResultsData = dat;
                gui.plotResultsData = 9;
            end
            
            % RAND %%%%
            funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
            %             dat = gui.(funcstring(10:end)).results;
            if  last_row>9
                % Results RAND
                if no_cols>cols && rows==9
                    funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                    %                     dat = gui.(funcstring(10:end)).results;
                    dat = act_results;
                    %                     Table_DisplayData(src,evnt,fh,dat{7,cols+1},cols)
                    gui.displayResultsData = dat{10,cols+1};
                    gui.plotResultsData = 10;
                end
                % Indiv. Results RAND
                if no_cols>cols && rows==10
                    funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                    %                     dat = gui.(funcstring(10:end)).results{8,cols+1};
                    dat = act_results{11,cols+1};
                    if isnumeric(dat)
                        if isnumeric(dat) && ndims(dat)<3
                            dat = dat(~all(isnan(dat),2),:);
                            %                         Table_DisplayData(src,evnt,fh,dat,cols)
                        elseif ndims(dat)>=3
                            gui.displayResultsData = {'Cannot display data'};
                            gui.plotResultsData = -1;
                        end
                        gui.displayResultsData = dat;
                        gui.plotResultsData = 11;
                    end
                end
                % Results (Distribution) RAND
                if no_cols>cols && rows==11
                    funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                    %                     Table_DisplayData(src,evnt,fh,dat{9,cols+1},cols)
                    dat = act_results;
                    gui.displayResultsData = dat{12,cols+1};
                    gui.plotResultsData = 12;
                end
                % Indiv. Results (Distribution) RAND
                if no_cols>cols && rows==12
                    funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
                    %                     dat = gui.(funcstring(10:end)).results{10,cols+1};
                    dat = act_results{13,cols+1};
                    if isnumeric(dat)
                        dat = dat(~all(isnan(dat),2),:);
                        %                         Table_DisplayData(src,evnt,fh,dat,cols)
                        gui.displayResultsData = dat;
                        gui.plotResultsData = 13;
                    end
                end
                if no_cols>cols && rows==13
                end
            end
            if no_cols>cols && rows==last_row % Remove set
                %                 while gui.blockDataSet
                %                     pause(.1)
                %                     gui=guidata(fh);
                %                 end
                %                 gui.blockDataSet = true;
                guidata(fh,gui);
                choice = questdlg(['Remove set' '?'], 'Remove set', ...
                    'Ok', ...
                    'Cancel','Cancel');
%                 set(gca,'FontSize',gui.fontsize_button)
%                 set(gcf,'WindowStyle','modal')
                if strcmpi(choice,'Ok')
                    data = get(src,'Data');
                    no_cols = size(get(src,'Data'),2);
                    colname = get(src,'Columnname');
                    data(:,cols) = [];
                    gui.input_data.(funcstring(10:end)).results(:,cols+1) = [];
                    for k=cols:no_cols-1
                        %                     colname{k}=['Set ' num2str(k)];
                        %                     if k>1
                        %                         actSetFromString = str2double(colname{k-1}(5:end))+1;
                        %                     else
                        %                         actSetFromString = 1;
                        %                     end
                        actSetFromString = k;
                        %                     colname{k}=colname{k+1};%['Set ' num2str(actSetFromString)];
                        
                        
                        if k<no_cols-1
                            
                            %                         colname{k}=colname{k+1};
                            if isfield(gui,'plot_mode') && isfield(gui.plot_mode.(funcstring(10:end)),['Set' num2str(actSetFromString+1)])
                                gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSetFromString)]) = ...
                                    gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSetFromString+1)]);
                            end
                            if isfield(gui.input_data.(funcstring(10:end)),['Set' num2str(actSetFromString+1)])
                                gui.input_data.(funcstring(10:end)).(['Set' num2str(actSetFromString)]) = ...
                                    gui.input_data.(funcstring(10:end)).(['Set' num2str(actSetFromString+1)]);
                                
                            end
                        end
                    end
                    colname(no_cols-1:end) = [];%repmat({''},1,numel(colname)-no_cols+2);
                    if isempty(colname)
                        colname{1} = 'New Set';
                    end
                    if isfield(gui,'plot_mode') && isfield(gui.plot_mode.(funcstring(10:end)),(['Set' num2str(no_cols-1)]))
                        gui.plot_mode.(funcstring(10:end)) = rmfield(gui.plot_mode.(funcstring(10:end)),(['Set' num2str(no_cols-1)]));
                    end
                    gui.input_data.(funcstring(10:end)) = rmfield(gui.input_data.(funcstring(10:end)),(['Set' num2str(no_cols-1)]));
                    gui.DataSetTableColumnNames = colname;
                    if size(data,2)==1
                        remIdx = find(cellfun(@(x) isempty(x), data(:,1)),1,'first');
                        data = data(1:remIdx,1);
                    end
                    set(src,'Data',data,'Columnname',colname);
                    
                    
                    
%                     pause(.9); % This seems to help against "multiple removals".
                    gui.(funcstring(10:end)).DataSetTableData = get(src,'Data');
                    %                 gui.(funcstring(10:end)).DataSetTableData
                    %                 gui.input_data.(funcstring(10:end)).results
                    %                 gui.input_data.(funcstring(10:end))
                    %                 waitfor(src,'Data');
                    no_cols = no_cols - 1; 
                end
            end
            if isfield(gui,'displayResultsData') && cols<no_cols
%                 gui.plotResultsData = gui.displayResultsData;
%                 gui.plotResultsIndivData = gui.displayResultsIndivData;
                
                if size(gui.input_data.(funcstring(10:end)).results,1)>9
                    tickidx = 11;
                else
                    tickidx = 7;
                end
                if iscell(act_results{tickidx,cols+1}) && ...
                        numel(gui.input_data.(funcstring(10:end)).results{tickidx,cols+1})>1
                    gui.plotResultsDataXTicks = gui.input_data.(funcstring(10:end)).results{tickidx,cols+1}{1};
                    gui.plotResultsDataYTicks = gui.input_data.(funcstring(10:end)).results{tickidx,cols+1}{2};
                    
                else
                    gui.plotResultsDataXTicks = gui.input_data.(funcstring(10:end)).results{tickidx,cols+1};
                end
                gui.plotResultsDataLabels = labels;
                set(gui.displayResultsButtonH,'Enable','on')
                set(gui.plotResultsButtonH,'Enable','on')
%                 set(gui.newFigureButtonH,'Enable','on')
                
            end
            
        elseif numel(cols)>1 && numel(rows)==1
            % %             gui.displayResultsData = [];
            %             set(gui.displayResultsButtonH,'Enable','off')
            %
            %             if all(no_cols>cols) && rows==10
            %                 funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
            %                 dat = gui.(funcstring(10:end)).results(3,cols+1);
            %                 figure; plot()
            %             end
        end
    catch
        keyboard
    end
end

actCellDat = get(src,'Data');
htmlStart = '<html><table border=0 width=400 bgcolor=#8abfe2 text-align: right;><TR><TD>';
htmlEnd = '</TD></TR> </table></right></html>';
if isfield(gui,'DataSetTableIndices') && ~isempty(gui.DataSetTableIndices) && ~isempty(act_indices) && act_indices(1)<last_row
    if ~isempty(actCellDat) && size(actCellDat,1)>=gui.DataSetTableIndices(1) && size(actCellDat,2)>=gui.DataSetTableIndices(2) && ...
            ~isempty(actCellDat(gui.DataSetTableIndices(1),gui.DataSetTableIndices(2))) && ...
            ~isempty(actCellDat{gui.DataSetTableIndices(1),gui.DataSetTableIndices(2)})
        %     actCellDat(gui.DataSetTableIndices(1),gui.DataSetTableIndices(2)) = ...
        %         {actCellDat{gui.DataSetTableIndices(1),gui.DataSetTableIndices(2)}(length(htmlStart)+1:end-length(htmlEnd))};
        try
            actCellDat(gui.DataSetTableIndices(1),gui.DataSetTableIndices(2)) = ...
                strrep(actCellDat(gui.DataSetTableIndices(1),gui.DataSetTableIndices(2)),htmlStart,'');
            actCellDat(gui.DataSetTableIndices(1),gui.DataSetTableIndices(2)) = ...
                strrep(actCellDat(gui.DataSetTableIndices(1),gui.DataSetTableIndices(2)),htmlEnd,'');
        catch
            keyboard
        end
    end
end
if ~isempty(act_indices) && act_indices(1)<last_row && ~(act_indices(2)==size(actCellDat,2) && act_indices(1))
    
    actCellDat{act_indices(1),act_indices(2)} = ...
        [htmlStart actCellDat{act_indices(1),act_indices(2)} htmlEnd];
    gui.DataSetTableIndices = act_indices;
    
end

% Deselect
set(src,'Data',cell(size(actCellDat)));
set(src,'Data', actCellDat );

guidata(fh,gui);
end

function cellDeSelectPlotMode(src,evnt,fh,gui)
idces = get(gui.plot_modeTable,'Data');
set(gui.plot_modeTable,'Data',cell(size(idces)));
set(gui.plot_modeTable,'Data',idces);
end

function cellSelectPlotMode(src,evnt,fh,gui,actSet)
% gui=guidata(fh);
row = get(src,'SelectedRows')+1;
col = get(src,'SelectedColumns');

idces = get(gui.plot_modeTable,'Data');


% single selection
if min(col)>0 && numel(col)==1 && numel(row)==1
    
    colname = get(gui.plot_modeTable,'ColumnName');
    timeIdx = find(strcmpi(colname,'time'));
    
    
    
    if ~isempty(timeIdx)
        notTimeIdx = 2:size(idces,2);
        
        idcesMat = NaN(size(idces,1),size(idces,2)-1);
        for k=1:size(idces,1)
            idcesMat(k,timeIdx-1)=gui.tmstring2idx{find(strcmp(idces(k,timeIdx),gui.tmstring2idx(:,1)),1,'first'),2};
            
        end
        if ~isempty(notTimeIdx)
            notTimeIdx(notTimeIdx==timeIdx) = [];
%             idces = idces([idces{:,1}]',2:end);
            idcesMat(:,notTimeIdx-1)= cell2mat(idces(:,notTimeIdx)); 
        end
    else
         idcesMat = cell2mat(idces(:,2:end));
    end
    
    
   
    % Act idx (selected)
    sameIdx = idcesMat(:,col)==idcesMat(row,col) & [idces{:,1}]';
    
    
    idces(~sameIdx,1) = repmat({false},[sum(~sameIdx),1]);
    idces(sameIdx,1) = repmat({true},[sum(sameIdx),1]);
end
% If various cells are selcted
if min(col)>0 && (numel(col)>1 || numel(row)>1)
    idces(:,1) = repmat({false},[size(idces,1),1]);
    
    idces(row,1) = repmat({true},[numel(row),1]);
end


% set(gui.plot_modeTable,'Data',cell(size(idces)));
set(gui.plot_modeTable,'Data',idces);

% single selection
if min(col)>0 && numel(col)==1 && numel(row)==1
    idces = idcesMat([idces{:,1}]',:);
end
% idces = idces([idces{:,1}]',2:end);
% idces = cell2mat(idces);
funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
colname = gui.colname;
sbidx=strcmp(colname,'Subset');
% if any(sbidx)
%     colname{sbidx} = 'Day';
% end

filterstring = {};
for k=1:size(idces,1)
    act_fstring=[];
    for k2=1:size(colname,2)
        if strcmp(colname{k2},'Time')
            act_fstring = [act_fstring colname{k2} gui.tmstring2idx(find(strcmp(idces(k,k2),gui.tmstring2idx(:,1)),1,'first'),2)];
        else
            if iscell(idces(k,k2)) % With time present, idces is a cell,else array
                act_fstring = [act_fstring colname{k2} idces(k,k2)];
            else
                act_fstring = [act_fstring colname{k2} {idces(k,k2)}];
            end
        end
    end
    filterstring = vertcat(filterstring,act_fstring);
end
% options = gui.funcOpts.(funcstring(10:end));
% if ~isempty(filterstring)
%     set(gui.plot_modeTableHOk,'Enable','on')
% end
plot_mode.data_filter = filterstring;
% actSet = gui.selectedSet;
gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]) = plot_mode;
% gui.input_data = eval([funcstring '(gui.input_data,options,plot_mode)']);
guidata(fh,gui);

end

function plot_modeTableOk(src,evnt,gui,fh,actSet)
% gui=guidata(fh);
% actSet = gui.selectedSet;
set(gui.funcH,'Enable','off')
set(gui.secH,'Enable','off')
set(gui.funcOptsH,'Enable','off')
set(gui.funcAdvOptsH,'Enable','off')
set(gui.funcHelpH,'Enable','off')
set(gui.DataSetTable,'Enable','off')
set(gui.preprocessButtonH,'Enable','off')
set(gui.export2workspaceButtonH,'Enable','off')
set(gui.displayResultsButtonH,'Enable','off')
% set(gui.newFigureButtonH,'Enable','off')
set(gui.trajectoryTable,'Enable','off')
set(gui.plotResultsButtonH,'Enable','off')
pause(.1)


funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);

% plot_mode_def = gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]);%gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);


idces = get(gui.plot_modeTable,'Data');
idces = idces([idces{:,1}]',2:end);
% ATTENTION: MIGHT BE NECESSARY,but for nowtry without idces = cell2mat(idces);
funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
colname = gui.colname;
sbidx=strcmp(colname,'Subset');
% if any(sbidx)
%     colname{sbidx} = 'Day';
% end


filterstring = {};
for k=1:size(idces,1)
    act_fstring=[];
    for k2=1:size(colname,2)
        if strcmp(colname{k2},'Time')
            act_fstring = [act_fstring colname{k2} gui.tmstring2idx(find(strcmp(idces(k,k2),gui.tmstring2idx(:,1)),1,'first'),2)];
        else
            if iscell(idces(k,k2)) % With time present, idces is a cell,else array
                act_fstring = [act_fstring colname{k2} idces(k,k2)];
            else
                act_fstring = [act_fstring colname{k2} {idces(k,k2)}];
            end
        end
    end
    filterstring = vertcat(filterstring,act_fstring);
end
plot_mode_def = gui.plot_modeDef.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);
gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]) = plot_mode_def;
gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).data_filter = filterstring;


plot_mode = gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]);
% if isfield(plot_mode_def,'xlabel') && isfield(plot_mode_def,'ylabel') && ...
%         ~isempty(plot_mode_def.xlabel)&& ~isempty(plot_mode_def.ylabel)
%     xaxString = {plot_mode_def.xlabel};
%     yaxString = {plot_mode_def.ylabel};
%
% else
xaxString = plot_mode_def.xaxis;
yaxString = plot_mode_def.yaxis;
% end
statString = plot_mode_def.statistics;
statString2 = get(gui.chooseStatH2,'String');%plot_mode_def.data;



if ~isempty(strtrim(xaxString(get(gui.chooseXaxisH,'Value')))); plot_mode.xaxis = xaxString(get(gui.chooseXaxisH,'Value')); end;
if ~isempty(yaxString) && ~isempty(strtrim(yaxString(get(gui.chooseYaxisH,'Value')))); plot_mode.yaxis = yaxString(get(gui.chooseYaxisH,'Value')); end;
if ~isempty(strtrim(statString(get(gui.chooseStatH,'Value')))); plot_mode.statistics = statString(get(gui.chooseStatH,'Value')); end;
if ~isempty(strtrim(statString2(get(gui.chooseStatH2,'Value')))); plot_mode.data = statString2(get(gui.chooseStatH2,'Value')); end;

if isfield(plot_mode_def,'normalization') && ~isempty(plot_mode_def.normalization) && ...
        (strcmpi(plot_mode_def.display_mode,'hist') || ...
        strcmpi(plot_mode_def.display_mode,'Map')|| ...
        strcmpi(plot_mode_def.display_mode,'MapPolar'))
    
    statString3 = plot_mode_def.normalization;
    if ~isempty(strtrim(statString3(get(gui.chooseStatH3,'Value')))); plot_mode.normalization = statString3(get(gui.chooseStatH3,'Value')); end;
end

gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).xaxis = plot_mode.xaxis;
if isfield(plot_mode,'yaxis')
    gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).xaxis = plot_mode.xaxis;
end
if ~isempty(yaxString); gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).yaxis = plot_mode.yaxis;end;
gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).data = plot_mode.data;
gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).statistics = plot_mode.statistics;

plot_mode.data_filter = gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).data_filter;
gui.plot_mode.(funcstring(10:end)).(['Set' num2str(actSet)]).display_mode = plot_mode_def.display_mode;
plot_mode.display_mode = plot_mode_def.display_mode;


options = gui.funcOptsNoCombos.(funcstring(10:end))(1);
fn = fieldnames(gui.funcAdvOptsNoCombos.(funcstring(10:end)));
for k = 1:size(fn,1)
    options.(fn{k}) = gui.funcAdvOptsNoCombos.(funcstring(10:end))(1).(fn{k});
end
options = idSocial_auxiliaries_generateDefOptions(options);

% if isfield(gui,funcstring(10:end)) && isfield(gui.(funcstring(10:end)),'results')
%     results_save = gui.(funcstring(10:end)).results;
% else
%     results_save = [];
% end
close(gcf)

try
    close(findall(0,'Tag','idSocialMessageBoxFigure'));
    
    progressBox = idSocial_auxiliaries_message; % idSocial_function_wrapper_FlatInput will look for the 'idSocialMessageBox'-Tag
    input_data = ...
        eval([funcstring '(gui.input_data,options,plot_mode)']);
    close(progressBox)
    gui.input_data = input_data;
    results_new = input_data.(funcstring(10:end)).results;
    if ~isfield(gui.input_data,funcstring(10:end)) || ~isfield(gui.input_data.(funcstring(10:end)),['Set' num2str(actSet)])
        gui.input_data.(funcstring(10:end)).(['Set' num2str(actSet)]) = ...
            input_data.(funcstring(10:end)).(['Set' num2str(actSet)]);
    end
    % results_new = input_data.(gui.act_method).results;
    
    
    
    results_new{1,2}=['Set ' num2str(actSet)];
    % if ~isempty(results_save)
    %     results_prev = results_save;
    %     results=results_prev;
    %     results(1:max(size(results,1),size(results_new,1)),actSet+1) = results_new(:,2);
    %     results(1:max(size(results,1),size(results_new,1)),1) = results_new(:,1);
    % else
    results = results_new;
    % end
    gui.(funcstring(10:end)).results = results;
    
    dat = get(gui.DataSetTable,'Data');
    dat(1:8,actSet) = {'Index Selection'; ...
        ['Results (' plot_mode.data{1} ')']; ...
        ['Indiv. Results(' plot_mode.statistics{1} ')' ]; ...
        ['Results (Distribution)']; ...
        ['Indiv. Results (Distribution)' ]; ...
        'Axis'; ...
        'Ticks'; ...
        'Tick Labels'; ...
        };
    dat(9,actSet) = {'<html><body bgcolor="#808080">Remove</body></html>'};
    if sum(cellfun(@(x) ~isempty(x),gui.(funcstring(10:end)).results(:,actSet+1)))>9
        dat(9:12,actSet) = {['Results (RAND)']; ...
            ['Indiv. Results(RAND)' ]; ...
            ['Results (Distribution,RAND)']; ...
            ['Indiv. Results (Distribution,RAND)' ]; ...
            };
        dat(13,actSet) = {'<html><body bgcolor="#808080">Remove</body></html>'};
        
    end
    
    dat(1,actSet+1) = {'<html><body bgcolor="#808080">Add</body></html>'};
    
    
    gui.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).DataSetTableData = dat;
    
    if numel(gui.DataSetTableColumnNames)==1 && strcmpi(gui.DataSetTableColumnNames{1},'New Set')
        gui.DataSetTableColumnNames = {'Set 1'};
    else
        if size(gui.DataSetTableColumnNames,1)> size(gui.DataSetTableColumnNames,2)
            gui.DataSetTableColumnNames = gui.DataSetTableColumnNames';
        end
        maxColName =  max(cellfun(@(x) str2double(x(5:end)),gui.DataSetTableColumnNames));
        gui.DataSetTableColumnNames = [gui.DataSetTableColumnNames, {['Set ' num2str(maxColName+1)]}];
    end
    set(gui.DataSetTable,'Data',dat,'ColumnName',gui.DataSetTableColumnNames);
    set(gui.export2workspaceButtonH,'Enable','on')
    set(gui.displayResultsButtonH,'Enable','on')
    set(gui.plotResultsButtonH,'Enable','on')
%     set(gui.newFigureButtonH,'Enable','on')
    set(gui.plotResultsButtonH,'Enable','on')
    set([gui.FigureTableButtonsH{:}],'Enable','on')

catch
    warndlg(['Sorry, could not execute ' ])
end
set(gui.funcH,'Enable','on')
set(gui.secH,'Enable','on')
set(gui.funcOptsH,'Enable','on')
set(gui.funcAdvOptsH,'Enable','on')
set(gui.funcHelpH,'Enable','on')
set(gui.DataSetTable,'Enable','on')
set(gui.preprocessButtonH,'Enable','on')
set(gui.export2workspaceButtonH,'Enable','on')
set(gui.displayResultsButtonH,'Enable','on')
% set(gui.newFigureButtonH,'Enable','on')
set(gui.trajectoryTable,'Enable','on')
% set(gui.plotResultsButtonH,'Enable','on')

guidata(fh,gui);
end

function updateDataSetTable(fh)
gui=guidata(fh);

funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);

if isfield(gui.input_data,funcstring(10:end)) && ~isempty(gui.input_data.(funcstring(10:end))) && ...
        isfield(gui.input_data.(funcstring(10:end)),'Set1') 
    
    fn = fieldnames(gui.input_data.(funcstring(10:end)));
    maxSet = -1;
    for k=1:size(fn,1)
        if ~isempty(strfind(fn{k},'Set'))
            maxSet_Temp = str2double(strrep(fn{k},'Set',''));
            if isfield(gui.input_data.(funcstring(10:end)).(fn{k}),'output_plot') && size(gui.input_data.(funcstring(10:end)).results,2)-1>=maxSet_Temp
                maxSet = max(maxSet,maxSet_Temp);
            end
        end
    end
    min_rows = inf;
    max_rows = -1;
    for actSet = 1:maxSet
        if isfield(gui.input_data.(funcstring(10:end)),'results')
            no_rows = sum(cellfun(@(x) ~isempty(x),gui.input_data.(funcstring(10:end)).results(:,actSet+1)));
            min_rows = min(min_rows,no_rows);
            max_rows = max(max_rows,no_rows);
        else
            maxSet = -1;
            min_rows = inf;
            max_rows = -1;
            
        end
            
        
    end
    no_rows = max(max_rows,min_rows);
    
    if maxSet>0
        dat = cell(no_rows,maxSet+1);
        colnames = cell(1,maxSet+1);
        
        for actSet = 1:maxSet
            
            plot_mode = gui.input_data.(funcstring(10:end)).(['Set' num2str(actSet)]).plot_mode;
            colnames{actSet} = ['Set ' num2str(actSet)];
            %         dat = get(gui.DataSetTable,'Data');
            dat(1:8,actSet) = {'Index Selection'; ...
                ['Results (' plot_mode.data{1} ')']; ...
                ['Indiv. Results(' plot_mode.statistics{1} ')' ]; ...
                ['Results (Distribution)']; ...
                ['Indiv. Results (Distribution)' ]; ...
                'Axis'; ...
                'Ticks'; ...
                'Tick Labels'; ...
                };
            dat(9,actSet) = {'<html><body bgcolor="#808080">Remove</body></html>'};
            if no_rows>9
                dat(9:12,actSet) = {['Results (RAND)']; ...
                    ['Indiv. Results(RAND)' ]; ...
                    ['Results (Distribution,RAND)']; ...
                    ['Indiv. Results (Distribution,RAND)' ]; ...
                    };
                dat(13,actSet) = {'<html><body bgcolor="#808080">Remove</body></html>'};
                
            end
            
            
            
        end
        colnames{maxSet+1} = 'New Set';
        dat(1,maxSet+1) = {'<html><body bgcolor="#808080">Add</body></html>'};
        gui.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).DataSetTableData = dat;
        gui.DataSetTableColumnNames = colnames;
        
        %     if numel(gui.DataSetTableColumnNames)==1 && strcmpi(gui.DataSetTableColumnNames{1},'New Set')
        %         gui.DataSetTableColumnNames = {'Set 1'};
        %     else
        %         if size(gui.DataSetTableColumnNames,1)> size(gui.DataSetTableColumnNames,2)
        %             gui.DataSetTableColumnNames = gui.DataSetTableColumnNames';
        %         end
        %         maxColName =  max(cellfun(@(x) str2double(x(5:end)),gui.DataSetTableColumnNames));
        %         gui.DataSetTableColumnNames = [gui.DataSetTableColumnNames, {['Set ' num2str(maxColName+1)]}];
        %     end
        set(gui.DataSetTable,'Data',dat,'ColumnName',gui.DataSetTableColumnNames);
        
    end
else
    dat(1,1) = {'<html><body bgcolor="#808080">Add</body></html>'};
    gui.DataSetTableColumnNames = {'New Set'};
    set(gui.DataSetTable,'Data',dat,'ColumnName',gui.DataSetTableColumnNames);

end
guidata(fh,gui);

end

function plot_modeTableCancel(src,evnt,fh)
gui=guidata(fh);
close(gcf)

guidata(fh,gui);
end

function displayResults(src,evnt,fh)
gui=guidata(fh);
if isfield(gui,'displayResultsData')
    dat = gui.displayResultsData;
    
    if iscell(dat) && numel(dat)>1 && all(cellfun(@(x) iscell(x)||isempty(x),dat(:)))
        if all(cellfun(@(x) numel(x)==1,dat)) && all(cellfun(@(x) isnumeric(x{1}),dat))
            dat2 = cell(numel(dat),max(cellfun(@(x) numel(x{1}),dat)));
            for k=1:numel(dat)
                if iscell(dat{k})
                    
                    dat2(k,1:numel(dat{k}{1}))= num2cell(dat{k}{1});
                    
                    
                end
            end
        else
            dat2 = cell(numel(dat),max(cellfun(@(x) numel(x),dat)));
            for k=1:numel(dat)
                if iscell(dat{k})
                    
                    dat2(k,1:numel(dat{k}))= dat{k};
                    
                elseif ~isempty(dat{k})
                    dat2(k,1:numel(dat{k}))= num2cell(dat{k});
                end
            end
        end
    elseif iscell(dat) && numel(dat)>1 && all(cellfun(@(x) iscell(x)||isempty(x),dat(:)))
        dat2 = cell(numel(dat),max(cellfun(@(x) numel(x),dat)));
        for k=1:numel(dat)
            if ~isempty(dat{k})
                dat2(k,1:numel(dat{k}))= num2cell(dat{k,:});
            end
        end
    elseif iscell(dat) && numel(dat)>1 && all(cellfun(@(x) isnumeric(x)||isempty(x),dat(:)))
        dat2 = cell(numel(dat),max(cellfun(@(x) numel(x),dat)));
        for k=1:numel(dat)
            if ~isempty(dat{k})
                dat2(k,1:numel(dat{k}))= num2cell(dat{k,:});
            end
        end
    else
        dat2=dat;
    end
    
    figure('Units','normalized','Position',[.45 .3 .3 .3],'Toolbar','none','Menubar','none','Name','Display data (Copy with ctrl-c)','NumberTitle','off');
    try
        uh = uitable('Data', dat2,'ColumnName','','RowName','','Units','normalized','Position',[0 0 1 1],'FontSize',12);
        if size(dat2,2)==1
            set(uh,'Columnwidth',{400})
        end
    catch
        warndlg('Sorry, cannot display the data');
    end
end
end

function plotResults(src,evnt,fh)
gui=guidata(fh);
rowIdx = gui.plotResultsData;

switch rowIdx
    case 3
        control = false;
    case 10
        control = true;
    otherwise
        control = false;
end

allFigs = get(0,'Children');
plotFigs = findobj(allFigs,'-regexp','Name','Plot Results');

% figure(plotFigs(1))
if ~isempty(plotFigs) && IsGraphicHandle(plotFigs(1))
    [gui.input_data, plotFigs] = idSocial_plotUI(gui.input_data,[gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)],gui.selectedSet,control,plotFigs(1));
else
    [gui.input_data, plotFigs] = idSocial_plotUI(gui.input_data,[gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)],gui.selectedSet,control);
end
   
guidata(fh,gui);

% % Update the figure table
updateFigTable(fh);


%     function closeFig(src,evnt)
%         ud = get(src,'UserData');
%         act_fig = ud.act_fig;
%         set(src,'CloseRequestFcn','closereq')
%         if verLessThan('matlab','8.2')
%             saveas(src, gui.input_data.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).(['Fig' num2str(act_fig)]).path,'fig');
%         else
%             savefig(src, gui.input_data.([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]).(['Fig' num2str(act_fig)]).path);
%         end
%         close
%     end

end

function updateFigTable(fh)

gui=guidata(fh);
% if nargin>2 && ~isempty(plotfigH)
%     figData = get(plotfigH,'UserData');
%     if ~isempty(figData) && isstruct(figData) && isfield(figData,'act_fig')
%         
%     end
% end

funcstring = ([gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)]);

    
    if isfield(gui.input_data,funcstring)
    fn = fieldnames(gui.input_data.(funcstring));
    fig_array = false(1,size(fn,1));
%     fig_array_no = NaN(1,size(fn,1));
    for k=1:size(fn,1)
        if ~isempty(strfind(fn{k},'Fig')) && ~isempty(gui.input_data.(funcstring).(fn{k})) && ...
                isstruct(gui.input_data.(funcstring).(fn{k}))
            fig_array(k) = true;
%             fig_array_no(k) = str2double(strrep(fn{k},'Fig',''));
        end
    end
    fn_figs = fn(fig_array);
    new_figDataCell = cell(sum(fig_array)+1,2);
    new_figCallbackCell = cell(sum(fig_array)+1,2);
 
    
    new_figDataCell{1,1} = 'New Figure';
    new_figDataCell{1,2} = '';
    new_figCallbackCell{1,1} = '';
    count = 2;
    for k = sum(fig_array):-1:1
        %         ud = get(gcf,'UserData');
        if isfield(gui.input_data.(funcstring).(fn_figs{k}),'act_set')
            strCell = strcat(cellstr(num2str(gui.input_data.(funcstring).(fn_figs{k}).act_set'))',',');
            if numel(gui.input_data.(funcstring).(fn_figs{k}).act_set) == 1
                strCell = ['Set ' strCell{:}];
            else
                strCell = ['Sets ' strCell{:}];
            end
            
            new_figDataCell{count,1} = ['Fig ' num2str(k) ': ' strCell(1:end-1) ];%fn_figs{k};
        else
            new_figDataCell{count,1} = ['Fig ' num2str(k)];%fn_figs{k};
        end
        
        %
        
        new_figDataCell{count,2} = 'X';
        new_figCallbackCell{count,1} = '';
        new_figCallbackCell{count,2} = '';
        count = count + 1;
    end
    else
        new_figDataCell = cell(1,2);
        new_figCallbackCell = cell(1,2);
        new_figDataCell{1,1} = 'New Figure';
        new_figDataCell{1,2} = '';
        new_figCallbackCell{1,1} = '';
        
    end
    
    [gui.FigureTableH, gui.FigureTableButtonsH] = ...
        idSocialUI_buttonTable(gui.FigureTableButtonsH,[],new_figDataCell,new_figCallbackCell,gui.FigureTableTooltips);
    
    
    % Have to repeat the deleteFigure Callbacks with the new
    % handles!
    set(gui.FigureTableButtonsH{1,1},'Callback',{@newFigure,fh});
    
    guidata(fh,gui);
    count = 2;
    for k = size(gui.FigureTableButtonsH,1)-1:-1:1
        propString = idSocial_auxiliaries_optstruct2string(gui.input_data,funcstring,k);
        set(gui.FigureTableButtonsH{count,1},'Callback',{@selectFigure,fh,k},'Tooltips',propString);
        set(gui.FigureTableButtonsH{count,2},'Callback',{@deleteFigure,fh,k});
        count = count + 1;

    end
    
%     gui=guidata(fh);
    
%     gui.FigureTableButtonsH = buttonHandles;
%     gui.input_data = input_data;
    guidata(fh,gui); % Did it work without this??!! 2017-03-30

end

function selectFigure(src,evt,fh,fig_no)
gui=guidata(fh);
funcstring=gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2);

open_figs = findall(0,'Type','figure');
fig_array_open = NaN(1,size(open_figs,1));
for k=1:size(open_figs,1)
    ud = get(open_figs(k),'UserData');
    if isfield(ud,'act_fig')
        fig_array_open(k) = ud.act_fig;
    end
end
open_figs = open_figs(~isnan(fig_array_open));
fig_array_open = fig_array_open(~isnan(fig_array_open));

if isfield(gui.input_data,funcstring) && isfield(gui.input_data.(funcstring),['Fig' num2str(fig_no)])
    
    
    
    % Figure is open, just bring it to front
%     if isfield(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]),'handle') && ...
%             IsGraphicHandle(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle)
%         figure(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle) 
    if any(fig_array_open==fig_no)
        
%         [gui.input_data, gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle] = ...
            idSocial_plotUI(gui.input_data,[gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)], ...
            [],[],open_figs(fig_array_open==fig_no));
        
    % Figure is not present, load from disk and adjust fig no    
    elseif isfield(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]),'path') && ...
            exist(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).path,'file')==2
        
%         [gui.input_data, gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle] = ...
            idSocial_plotUI(gui.input_data,[gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)], ...
            [],[],fig_no);
%         
%         gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle = ...
%             openfig(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).path);
%         ud = get(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle,'UserData');
%         ud.act_fig = fig_no;
%         set(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle,'UserData',ud)
%         set(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle,'CloseRequestFcn',{@closeFig2})
        
    end
%     set(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle,'Name',['Plot Results ' num2str(fig_no)]);
    
    
end

guidata(fh,gui);

%     function closeFig2(src,evnt)
%         
%         uds2 = get(src,'UserData');
%         save_path = uds2.save_path;
%         if verLessThan('matlab','8.2')
%             set(src,'CloseRequestFcn','')
%             saveas(src, save_path,'fig');
%             set(src,'CloseRequestFcn',{@closeFig2})
%         else
%             set(src,'CloseRequestFcn','')
%             uds2.act_fig=[];
%             set(src,'UserData',uds2);
%             get(src)
%             savefig(src, save_path);
%             set(src,'CloseRequestFcn',{@closeFig2})
%         end
%         delete(src)
%         
%     end


end

function deleteFigure(src,evt,fh,fig_no)
gui=guidata(fh);

funcstring=gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2);

if isfield(gui.input_data,funcstring) && isfield(gui.input_data.(funcstring),['Fig' num2str(fig_no)])
    if isfield(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]),'handle') && ...
            IsGraphicHandle(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle)
        delete(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).handle)
    end
    if isfield(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]),'path') && ...
            exist(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).path,'file')==2
        delete(gui.input_data.(funcstring).(['Fig' num2str(fig_no)]).path)
    end
    gui.input_data.(funcstring) = rmfield(gui.input_data.(funcstring),['Fig' num2str(fig_no)]);

end

% Check how many figures there are now
fn = fieldnames(gui.input_data.(funcstring));
fig_array = false(1,size(fn,1));
for k=1:size(fn,1)
    if ~isempty(strfind(fn{k},'Fig')) && ~isempty(gui.input_data.(funcstring).(fn{k})) && ...
            isstruct(gui.input_data.(funcstring).(fn{k}))
        fig_array(k) = true;
    end
end


open_figs = findall(0,'Type','figure');
fig_array_open = NaN(1,size(open_figs,1));
for k=1:size(open_figs,1)
    ud = get(open_figs(k),'UserData');
    if isfield(ud,'act_fig')
        fig_array_open(k) = ud.act_fig;
    end
end
open_figs = open_figs(~isnan(fig_array_open));
fig_array_open = fig_array_open(~isnan(fig_array_open));

count=1;
% Sort fig fields so they run from 1 to max without gaps and
% move files
for k = 1:numel(fig_array)
    if isfield(gui.input_data.(funcstring),['Fig' num2str(k)])
        if count ~= k
            gui.input_data.(funcstring).(['Fig' num2str(count)]) = gui.input_data.(funcstring).(['Fig' num2str(k)]);
            gui.input_data.(funcstring) = rmfield(gui.input_data.(funcstring),['Fig' num2str(k)]);
            
%             if isfield(gui.input_data.(funcstring).(['Fig' num2str(count)]),'handle') && ...
%                     IsGraphicHandle(gui.input_data.(funcstring).(['Fig' num2str(count)]).handle)
%                 ud = get(gui.input_data.(funcstring).(['Fig' num2str(count)]).handle,'UserData');
%                 ud.act_fig = count;
%                 set(gui.input_data.(funcstring).(['Fig' num2str(count)]).handle,'UserData',ud)
%             end
            if any(k == fig_array_open)
                ud = get(open_figs(k == fig_array_open),'UserData');
                ud.act_fig = count;
                set(open_figs(k == fig_array_open),'UserData',ud)
            end
            
            if isfield(gui.input_data.(funcstring).(['Fig' num2str(count)]),'path') && ...
                    exist(gui.input_data.(funcstring).(['Fig' num2str(count)]).path,'file')==2
                [~, temp_savepath]= idSocial_recursiveGetOptionsFromOptionsCell(gui.input_data.options,'temp_savepath');

                movefile(gui.input_data.(funcstring).(['Fig' num2str(count)]).path,...
                    [temp_savepath funcstring '_figure' num2str(count) '.fig']);
                gui.input_data.(funcstring).(['Fig' num2str(count)]).path = ...
                    [temp_savepath funcstring '_figure' num2str(count) '.fig'];
%                 ud = get(gui.input_data.(funcstring).(['Fig' num2str(count)]).handle,'UserData');
%                 ud.save_path =  [temp_savepath funcstring '_figure' num2str(count) '.fig'];
%                 set(gui.input_data.(funcstring).(['Fig' num2str(count)]).handle,'UserData',ud)

            end
            
            
            
        end
        
        count = count + 1;
        
    end
end



guidata(fh,gui);
updateFigTable(fh);
end

function newFigure(src,evnt,fh)
gui=guidata(fh);

rowIdx = gui.plotResultsData;

switch rowIdx
    case 3
        control = false;
    case 10
        control = true;
    otherwise
        control = false;
end

[gui.input_data, plotFigs] = ...
    idSocial_plotUI(gui.input_data,[gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2)],gui.selectedSet,control);
% plotfh = figure('Units','normalized','Position',[.45 .3 .5 .5],'NumberTitle','off','Color','w');
% 
% funcstring=gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(10:end-2);
% fn = fieldnames(gui.input_data.(funcstring));
% maxFig = 0;
% for k=1:size(fn,1)
%     if ~isempty(strfind(fn{k},'Fig'))
%         maxFig = max(maxFig,str2double(strrep(fn{k},'Fig','')));
%     end
% end
% act_fig = maxFig+1;
% [~, temp_savepath]= idSocial_recursiveGetOptionsFromOptionsCell(gui.input_data.options,'temp_savepath');
% 
% set(plotfh,'Name', ['Plot Results ' num2str(act_fig)])
% ud.act_fig = act_fig;
% ud.save_path =  [temp_savepath funcstring '_figure' num2str(act_fig) '.fig'];
% set(plotfh,'UserData',ud)
% % gui.figure_windows = [gui.figure_windows plotfh];
% 
% % set(plotfh,'CloseRequestFcn',{@closeFig})
% gui.input_data.(funcstring).(['Fig' num2str(act_fig)]).handle = plotfh;
% gui.input_data.(funcstring).(['Fig' num2str(act_fig)]).path = ud.save_path;
% 
% if verLessThan('matlab','8.2')
% %     set(plotfh,'CloseRequestFcn','')
%     saveas(plotfh, gui.input_data.(funcstring).(['Fig' num2str(act_fig)]).path,'fig');
% %     set(plotfh,'CloseRequestFcn',{@closeFig})
% else
%     set(plotfh,'CloseRequestFcn','')
%     %         pause(.1)
%     ud.act_fig=[];
%     set(plotfh,'UserData',ud);
%     savefig(plotfh, gui.input_data.(funcstring).(['Fig' num2str(act_fig)]).path);
%     set(plotfh,'CloseRequestFcn',{@closeFig})
% end


guidata(fh,gui);

updateFigTable(fh);
%     function closeFig(src,evnt)
%         
%         uds2 = get(src,'UserData');
%         save_path = uds2.save_path;
%         if verLessThan('matlab','8.2')
%             set(src,'CloseRequestFcn','')
%             saveas(src, save_path,'fig');
%             set(src,'CloseRequestFcn',{@closeFig})
%         else
%             set(src,'CloseRequestFcn','')
%             uds2.act_fig=[];
%             set(src,'UserData',uds2);
%             get(src)
%             savefig(src, save_path);
%             set(src,'CloseRequestFcn',{@closeFig})
%         end
%         delete(src)
%         
%     end


end



% function closeFig(src,evnt)
% 
% uds = get(src,'UserData');
% act_fig = uds.act_fig;
% set(src,'CloseRequestFcn','')
% 
% if verLessThan('matlab','8.2')
%     set(src,'CloseRequestFcn','')
%     saveas(src, gui.input_data.(funcstring).(['Fig' num2str(act_fig)]).path,'fig');
%     set(src,'CloseRequestFcn',{@closeFig})
% else
%     set(src,'CloseRequestFcn','')
%     %         pause(.1)
%     uds.act_fig=[];
%     set(src,'UserData',uds);
%     savefig(src, gui.input_data.(funcstring).(['Fig' num2str(act_fig)]).path);
%     set(src,'CloseRequestFcn',{@closeFig})
% end
% % close
% delete(src)
% 
% end

function Table_DisplayData(src,evnt,fh,dat,cols)
gui=guidata(fh);

gui.displayResultsData = dat;
guidata(fh,gui);
end

function export2workspace(src,evnt,fh)
gui=guidata(fh);
if isfield(gui,'out_results')
    out_results = gui.out_results;
end
funcstring = gui.funcCell{gui.selectedSection,1}{gui.selectedFunc}(1:end-2);
out_results.(funcstring(10:end)) = gui.(funcstring(10:end)).results;
disp(['idSocial_results exported to workspace. New results can be found in idSocial_results.' funcstring(10:end) ':'])
assignin('base','idSocial_results',out_results)
out_results.(funcstring(10:end))
guidata(fh,gui);
end

function onHeaderClick(src, evt, jtable, fh, pmTable,tmstring2idx)
gui=guidata(fh);
% if(get(evt, 'ClickCount') > 1)
col = jtable.convertColumnIndexToModel(src.columnAtPoint(evt.getPoint())) + 1;
dat = get(pmTable,'Data');
colw = get(pmTable,'ColumnWidth');
if iscell(dat)
    tcol = find(cellfun(@(x) ischar(x),dat(1,:)));
else
    tcol=[];
end
%     tmcol =ischar(dat{1,col});
if col>1
    if ~isempty(tcol)
        tmstr = unique(tmstring2idx(:,1));
        str2tm = cell(numel(tmstr),1);
        for k=1:numel(tmstr)
            str2tm(k,1) = tmstr(k);
            idx = strfind(tmstring2idx(:,1),tmstr{k});
            str2tm(k,2) = tmstring2idx(find(cellfun(@(x) ~isempty(x),idx),1,'first'),2);
        end
        for k = 1:size(dat,1)
            dat{k,tcol} = str2tm{strcmp(dat{k,tcol},str2tm(:,1)),2};
        end
        dat=cell2mat(dat);
    end
    [sortDat sortIdx] = sortrows(dat,col);
    if isequal(sortDat(:,col),dat(:,col)) % Check if already sorted
        %         if all(diff(sortDat(:,col))>=0) && numel(unique(sortDat(:,col)))>2 % Check if  ascending or descending
        %         if  numel(unique(sortDat(:,col)))>2 % Check if  ascending or descending
        sortDat = sortDat(end:-1:1,:);
        sortIdx = sortIdx(end:-1:1,:);
        %         elseif all(diff(sortDat(:,col))>=0) && numel(unique(sortDat(:,col)))<=2
        %         elseif numel(unique(sortDat(:,col)))<=2
        %             if isequal(dat(:,col),sortDat(:,col))
        %                 sortDat = sortDat(end:-1:1,:);
        %             end
        %         end
    end
    if ~isempty(tcol)
        sortDat= num2cell(sortDat);
        
        for k = 1:size(sortDat,1)
            sortDat{k,tcol} = str2tm{sortDat{k,tcol}==[str2tm{:,2}],1};
        end
    end
    set(pmTable,'Data',sortDat)
    set(pmTable,'ColumnWidth',colw);
elseif col==1
    if all(~[dat{:,1}])
        dat(:,1) = repmat({true},[size(dat,1),1]);
    else
        dat(:,1) = repmat({false},[size(dat,1),1]);
    end
    set(pmTable,'Data',dat);
end
% end
guidata(fh,gui);

end

function plot_modeTableContext(src, evt,fh)
gui=guidata(fh);


idxcombs = NaN(1,6);
for gr=1:gui.input_data.info.no_groups
    for sb=1:gui.input_data.info.no_subsets(gr)
        for tr=1:gui.input_data.info.no_trials(gr,sb)
            for tm=tsteps
                for fc=1:gui.input_data.info.no_focals(gr,sb,tr)
                    for nb=1:gui.input_data.info.no_neighbors(gr,sb,tr)
                        idxcombs = vertcat(idxcombs,[gr, sb, tr, tm, fc, nb]);
                    end
                end
            end
        end
    end
end
idxcombs(1,:)=[];
gr = idxcombs(:,1);
sb = idxcombs(:,2);
tr = idxcombs(:,3);
tm = idxcombs(:,4);
fc = idxcombs(:,5);
nb = idxcombs(:,6);


guidata(fh,gui);
end

function saveData(src, evt,fh,prev)
gui=guidata(fh);
%([gui.recentPath '*.mat']);
if prev<2 || ~isfield(gui,'savepath') || isempty(gui.savepath)
    if isfield(gui,'savepath') && ~isempty(gui.savepath)
        if ~strcmp(gui.savepath(end),filesep);
            savep = [gui.savepath filesep];
        end
%         [FileName,PathName] = uiputfile([savep '*.fig', ' MATLAB GUI']);
          [FileName,PathName] = uiputfile([savep '*.mat']);
        
    else
%         [FileName,PathName] = uiputfile(['*.fig', ' MATLAB GUI']);
        [FileName,PathName] = uiputfile(['*.mat']);

    end
    if ischar(PathName) && ischar(FileName)
        if ~strcmp(FileName(end-3:end),'.mat')
            FileName = [FileName '.mat'];
        end
        gui.savepath = [PathName FileName(1:end-4) '.mat'];
        
        
    end
end

if isfield(gui,'savepath') && ~isempty(gui.savepath)
    % save([PathName 'Data' FileName],'gui');
    guidata(fh,gui);
    gui=guidata(fh);
    input_data = gui.input_data;
    try
        save(gui.savepath,'input_data');
    catch
        warning(['Could not save ' gui.savepath '. Sorry!'])
    end
%     try
%         if verLessThan('matlab','8.2')
%             saveas(fh,savep)
%         else
%             savefig(fh,savep)
%         end
%         disp(['Saved ' savep])
%     catch
%         warning(['Could not save ' savep '. Matlab 2013b or higher is required. Sorry!'])
%     end
end
% fh2 = openfig([PathName FileName(1:end-4) '.fig'],'reuse');
end

function loadData(src, evt,fh)
gui=guidata(fh);


[FileName,PathName] = uigetfile({'*.mat';'*.fig';});

if ischar(PathName) && ischar(FileName) && exist([PathName FileName(1:end-4) '.fig'],'file') == 2
    pos = get(fh,'Position');
    close(fh)
    fh = openfig([PathName FileName(1:end-4) '.fig'],'reuse');
    gui=guidata(fh);
    % fh=fh2;
    set(fh,'Position',pos);
    % load([PathName FileName]);
elseif ischar(PathName) && ischar(FileName) && exist([PathName FileName(1:end-4) '.mat'],'file') == 2
    pos = get(fh,'Position');
    load_inputData(src,evt,fh,[PathName FileName(1:end-4) '.mat'])
    gui=guidata(fh);
    % fh=fh2;
    set(fh,'Position',pos);
else
    warning(['Could not load' PathName FileName(end-3:end) '.fig'])
end
guidata(fh,gui);
end

function setPath(src, evt,fh)
if exist([pwd filesep 'idSocial.conf'],'file')==2
    fileID = fopen([pwd filesep 'idSocial.conf'],'r');
    tempString = textscan(fileID,'%s','Delimiter','');
    fclose(fileID);
    temp_savepath = tempString{1}{1};
else
    temp_savepath = '';
end
temp_savepath = uigetdir(temp_savepath, ...
    sprintf('Please select directory for idSocial projects'));

% Open idSocial.conf and write folder information if the folder exists. If not
if ~(isempty(temp_savepath) || temp_savepath==0) && exist(temp_savepath,'dir')==7
    fileID = fopen([pwd filesep 'idSocial.conf'],'w');
    fprintf(fileID,'%s',temp_savepath);
    fclose(fileID);
end
if exist([pwd filesep 'idSocial.conf'],'file')~=2 || isempty(temp_savepath) || temp_savepath==0
    if exist(temp_savepath,'dir')~=7
        warning('Could not find folder for temporary files.')
        if exist([pwd filesep 'idSocial.conf'],'file')==2
            delete([pwd filesep 'idSocial.conf'])
        end
    end
    
    
end
end

% Function "outsourced"
% function parameters = handleComboMenus(parameters)
% % Handle combo menus
% % parameters = gui.defPreOpts;
% optNames = fieldnames(parameters(1));
% optVals = struct2cell(parameters(1));
% isCombo = cellfun(@(x) iscell(x)&size(x,2)>1,optVals);
% optCombos = optVals(isCombo);
% optComboName = optNames(isCombo);
% 
% for k=1:sum(isCombo)
% %     all(cellfun(@(x) ischar(x),optCombos{k}(1:end-1))) && ...
%     if all(cellfun(@(x) isequal(class(x),class(optCombos{k}{1})),optCombos{k}(1:end-1))) && ...
%             iscell(optCombos{k}{end}) && all(cellfun(@(x) isequal(class(x),class(optCombos{k}{end}{1})),optCombos{k}(1:end-1)))
%         
%         optCombos{k} = optCombos{k}{end};%optCombos{k}{optCombos{k}{end}{1}};
%         parameters(1).(optComboName{k}) = optCombos{k};
%     elseif any(cellfun(@(x) iscell(x),optCombos{k})) % combined "combo -non-combo" menu
%         
%         isNestedCombo = find(cellfun(@(x) iscell(x),optCombos{k}));
% %         optCombosNested = optCombos{k}(isNestedCombo);
%         for k2=isNestedCombo
%             if all(cellfun(@(x) isequal(class(x),class(optCombos{k}{k2}{1})),optCombos{k}{k2}(1:end-1))) && ...
%                 iscell(optCombos{k}{k2}{end}) && all(cellfun(@(x) isequal(class(x),class(optCombos{k}{k2}{end}{1})),optCombos{k}{k2}(1:end-1)))
%                 
% %                 iscell(optCombos{k}{k2}{end}) && ischar(optCombos{k}{k2}{end}{1})
%                 optCombos{k}{k2} = optCombos{k}{k2}{end}{1};%optCombos{k}{optCombos{k}{end}{1}};
%             else
%                 optCombos{k}{k2} = optCombos{k}{k2}{1};
%             end
%         end
%         parameters(1).(optComboName{k}) = optCombos{k};
%     else
%         parameters(1).(optComboName{k}) = optCombos{k}{1};
%     end
% end
% end
function opts = RestoreComboOptions(opts,defOpts)
% Handle combo menus
% parameters = gui.defPreOpts;
optNames = fieldnames(defOpts(1));
optVals = struct2cell(defOpts(1));
isCombo = cellfun(@(x) iscell(x)&size(x,2)>1,optVals);
optCombos = optVals(isCombo);
optComboName = optNames(isCombo);

for k=1:sum(isCombo)
    if iscell(optCombos{k}) && isfield(opts,optComboName{k})
        opts.(optComboName{k}) = [optCombos{k} {{opts.(optComboName{k})}}];
%         idx = strcmp(optCombos{k},opts.(optComboName{k}));
%         opts.(optComboName{k}) = [opts.(optComboName{k}) optCombos{k}(~idx)];

    end
end
end