% PathName = 'C:\Users\Robert\Syncplicity Folders\Onset3\data2\';
% ids = {'N2','day','WT','trajectories.mat'};
% sfilters = {};
% countable = [];%{false,false,true,false};
% countable = {false,true,false};
% cfilter = {[],[5 25],[]};
% 
% [good_files,idxOut] = idSocial_createFileTreeFromPath(PathName,ids,countable,cfilter,sfilters) 
%%
function folderTreeFig = idSocialUI_chooseFolderTree(parentFig)
if nargin<1 || isempty(parentFig)
    parentFig = [];
end

fontsize_button = 11;

no_parts=5;
margin = .05;
part_height = (1-6*margin)/5;
width = .9;

gridy = margin:(part_height+margin):1-(part_height+margin);


bheight = part_height*.9;
bwidth = .3;

fontsize = 12;


folderTreeFig=...
    figure('Units','normalized','Position',[.3 .3 .3 .4],'Name','Select trajectory folders', 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off','WindowStyle','modal');

folderTreeSetPathH = uicontrol(folderTreeFig,'Units','normalized','Style','pushbutton','String','Choose folder', ...
    'Callback',@(src,evnt) folderTreeSetPath(src,evnt),'FontSize',fontsize);
set(folderTreeSetPathH,'Position', [.05 gridy(5) bwidth bheight])

folderTreePathTextH = uicontrol(folderTreeFig,'Units','normalized','Style','text','String','...', ...
    'FontSize',fontsize);
set(folderTreePathTextH,'Position', [.05+bwidth gridy(5) width-bwidth-2*margin bheight])

% First identifier
folderTreeId1PanelH = uipanel('Title','First Identifier','FontSize',fontsize);
set(folderTreeId1PanelH,'Position', [.05 gridy(4) width part_height])

folderTreeId1H = uicontrol(folderTreeId1PanelH,'Units','normalized','Style','edit','String','', ...
    'FontSize',fontsize);
set(folderTreeId1H,'Position', [margin margin (1-4*margin)/3+.15 1-2*margin])

folderTreeSetId1NumH = uicontrol(folderTreeId1PanelH,'Units','normalized','Style','checkbox','String','Numbered', ...
    'Callback',@(src,evnt) enableLim(src,evnt,1),'FontSize',fontsize);
set(folderTreeSetId1NumH,'Position', [margin+margin+(1-4*margin)/3+.16 margin (1-4*margin)/3 1-2*margin])

folderTreeSetId1Lim1H = uicontrol(folderTreeId1PanelH,'Units','normalized','Style','edit','String','Lim1', ...
    'FontSize',fontsize,'Enable','off');
set(folderTreeSetId1Lim1H ,'Position', [3*margin+2*(1-4*margin)/3+.04 margin ((1-4*margin)/3-margin)/2 1-2*margin])
folderTreeSetId1Lim2H = uicontrol(folderTreeId1PanelH,'Units','normalized','Style','edit','String','Lim2', ...
    'FontSize',fontsize,'Enable','off');
set(folderTreeSetId1Lim2H,'Position', [4*margin+2*(1-4*margin)/3+((1-4*margin)/3-margin)/2 margin ((1-4*margin)/3-margin)/2 1-2*margin])

% Second identifier
folderTreeId2PanelH = uipanel('Title','Second Identifier','FontSize',fontsize);
set(folderTreeId2PanelH,'Position', [.05 gridy(3) width part_height])

folderTreeId2H = uicontrol(folderTreeId2PanelH,'Units','normalized','Style','edit','String','', ...
    'FontSize',fontsize);
set(folderTreeId2H,'Position', [margin margin (1-4*margin)/3+.15 1-2*margin])

folderTreeSetId2NumH = uicontrol(folderTreeId2PanelH,'Units','normalized','Style','checkbox','String','Numbered', ...
    'Callback',@(src,evnt) enableLim(src,evnt,2),'FontSize',fontsize);
set(folderTreeSetId2NumH,'Position', [margin+margin+(1-4*margin)/3+.16 margin (1-4*margin)/3 1-2*margin])

folderTreeSetId2Lim1H = uicontrol(folderTreeId2PanelH,'Units','normalized','Style','edit','String','Lim1', ...
    'FontSize',fontsize,'Enable','off');
set(folderTreeSetId2Lim1H ,'Position', [3*margin+2*(1-4*margin)/3+.04 margin ((1-4*margin)/3-margin)/2 1-2*margin])
folderTreeSetId2Lim2H = uicontrol(folderTreeId2PanelH,'Units','normalized','Style','edit','String','Lim2', ...
    'FontSize',fontsize,'Enable','off');
set(folderTreeSetId2Lim2H,'Position', [4*margin+2*(1-4*margin)/3+((1-4*margin)/3-margin)/2 margin ((1-4*margin)/3-margin)/2 1-2*margin])

% Filters
folderTreeId1PanelH = uipanel('Title','Additional Filters (use "," to separate)','FontSize',fontsize);
set(folderTreeId1PanelH,'Position', [.05 gridy(2) width part_height])

folderTreeFiltersH = uicontrol(folderTreeId1PanelH,'Units','normalized','Style','edit','String','', ...
    'FontSize',fontsize);
set(folderTreeFiltersH ,'Position', [margin margin (1-2*margin) 1-2*margin])


% Ok, cancel

folderTreeSetCancelH = uicontrol(folderTreeFig,'Units','normalized','Style','pushbutton','String','Cancel', ...
    'Callback','close(gcf)','FontSize',fontsize);
set(folderTreeSetCancelH,'Position', [width/2 - bwidth margin  bwidth bheight])

folderTreeSetOkH = uicontrol(folderTreeFig,'Units','normalized','Style','pushbutton','String','Ok', ...
    'Callback',@(src,evnt) folderTreeSetOk(src,evnt,folderTreeFig),'FontSize',fontsize);
set(folderTreeSetOkH,'Position', [width/2 + 2*margin  margin  bwidth bheight])

PathName = [];
% guidata(parentFig,gui);

%%%%%% Callbacks
    function folderTreeSetPath(src,evnt)
        PathName = uigetdir;
        if ~isempty(PathName)
            set(folderTreePathTextH,'String',PathName)
        end
    end
    function folderTreeSetOk(src,evnt,fh)
%         PathName = 'C:\Users\Robert\Syncplicity Folders\Onset3\data2\';
%         ids = {get(folderTreeId1H,'String'), get(folderTreeId2H,'String')};%{'N2','day','WT','trajectories.mat'};
        ids = cell(1,2);
        ids(1) = {get(folderTreeId1H,'String')};
        ids(2) = {get(folderTreeId2H,'String')};

        sfilters = get(folderTreeFiltersH,'String');
        if ~isempty(sfilters)
            sfilters = textscan(sfilters,'%s','Delimiter',',');
            sfilters =sfilters{1};
        end
        
        countable = {get(folderTreeSetId1NumH,'Value'),get(folderTreeSetId2NumH,'Value')};
        
        cf11 = get(folderTreeSetId1Lim1H,'String');
        cf12 = get(folderTreeSetId1Lim2H,'String');
        cf21 = get(folderTreeSetId2Lim1H,'String');
        cf22 = get(folderTreeSetId2Lim2H,'String');
        cfilter = {[],[]};
        if ~isempty(cf11) && ~isempty(cf12) 
            cf11 = str2double(cf11);
            cf12 = str2double(cf12);
            if ~isnan(cf11) && ~isnan(cf12) && isnumeric(cf11) && isnumeric(cf12)
                cfilter{1} = [cf11 cf12];
            end
        end
        
        if ~isempty(cf21) && ~isempty(cf22) 
            cf21 = str2double(cf21);
            cf22 = str2double(cf22);
            if ~isnan(cf21) && ~isnan(cf22) && isnumeric(cf21) && isnumeric(cf22)
                cfilter{2} = [cf21 cf22];
            end
        end
        
        [good_files,idxOut] = idSocial_createFileTreeFromPath(PathName,ids,countable,cfilter,sfilters);
        set(fh,'Visible','off');
        folderTreeDisplayList(src,evnt,parentFig,folderTreeFig,[num2cell(idxOut),good_files]);
    end

    function enableLim(src,evnt,id)
        if id==1
            if get(src,'Value')
                set(folderTreeSetId1Lim1H,'Enable','on','String','')
                set(folderTreeSetId1Lim2H,'Enable','on','String','')
            else
                set(folderTreeSetId1Lim1H,'Enable','off','String','Lim1')
                set(folderTreeSetId1Lim2H,'Enable','off','String','Lim2')
            end
            
        elseif id==2
            if get(src,'Value')
                set(folderTreeSetId2Lim1H,'Enable','on','String','')
                set(folderTreeSetId2Lim2H,'Enable','on','String','')
            else
                set(folderTreeSetId2Lim1H,'Enable','off','String','Lim1')
                set(folderTreeSetId2Lim2H,'Enable','off','String','Lim2')
            end
        end
    end



end

function folderTreeDisplayList(src,evnt,parentFig,folderTreeFig,act_data)
% gui=guidata(parentFig);
% act_data = gui.folderTreeDisplayList;
fontsize = 12;

% act_data(:,1:3) = cellfun(@(x) num2str(x),act_data(:,1:3),'UniformOutput',false); % Only for formatting reasons.

folderTreeDisplayFig=...
    figure('Units','normalized','Position',[.3 .3 .5 .4],'Name','Please check found files', 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off','WindowStyle','modal');

set(folderTreeDisplayFig,'Units','pixels');
figPos = get(folderTreeDisplayFig,'Position');
set(folderTreeDisplayFig,'Units','normalized');



columnname = {'Remove','Group','Subgroup','Trial','File'};
columnformat = {'logical','char','char','char','char'};
columnwidth = {figPos(3)*.3/4 figPos(3)*.3/4 figPos(3)*.3/4 figPos(3)*.3/4 figPos(3)*.7};
folderTreeSetPathH = uitable('Data',[repmat({true},[size(act_data,1) 1]) act_data], ...
    'Units','normalized', ...
    'Position',[0 .15 1 .85], ...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [true false false false false],...
    'ColumnWidth',columnwidth,...
    'RowName',[],...
    'Fontsize',fontsize);
    
if ~isprop(folderTreeDisplayFig,'SizeChangedFcn') && isprop(folderTreeDisplayFig,'ResizeFcn')
    set(folderTreeDisplayFig,'ResizeFcn',@(src,evnt) folderTreeDisplaySiz(src,evnt,folderTreeSetPathH));
else
    set(folderTreeDisplayFig,'SizeChangedFcn',@(src,evnt) folderTreeDisplaySiz(src,evnt,folderTreeSetPathH));
end
folderTreeDisplayCancelH = uicontrol(folderTreeDisplayFig,'Units','normalized','Style','pushbutton','String','Cancel', ...
    'Callback',@(src,evnt) folderTreeDisplayNewListCancel(src,evnt),'Fontsize',fontsize);
set(folderTreeDisplayCancelH,'Position', [.2 .01 .2 .15-.02])

folderTreeDisplayOkH = uicontrol(folderTreeDisplayFig,'Units','normalized','Style','pushbutton','String','Ok', ...
    'Callback',@(src,evnt) folderTreeDisplayNewList(src,evnt,folderTreeSetPathH),'Fontsize',fontsize);
set(folderTreeDisplayOkH,'Position', [1-.2-.2 .01 .2 .15-.02])

    function folderTreeDisplaySiz(src,evnt,folderTreeSetPathH)
            set(src,'Units','pixels'); 
            figPos = get(src,'Position');
            set(src,'Units','normalized');
            colwidth = get(folderTreeSetPathH,'ColumnWidth');
            set(folderTreeSetPathH,'ColumnWidth',{colwidth{1},colwidth{2} colwidth{3} colwidth{4} figPos(3)-colwidth{1}-colwidth{2}-colwidth{3}-colwidth{4}});
            
            
    end
    function folderTreeDisplayNewList(src,evnt,folderTreeSetPathH)
        gui = guidata(parentFig);
        act_dataNew = get(folderTreeSetPathH,'Data');
        checked = [act_dataNew{:,1}];
        folderTreeDisplayList = act_data(checked,:);
        
        folderTreeDisplayList = sortrows(folderTreeDisplayList,[1 2 3]);
        
        disp(folderTreeDisplayList)
        gui.openFolderFileChooserList=folderTreeDisplayList;
        guidata(parentFig,gui);
        close(folderTreeDisplayFig)
        close(folderTreeFig)
    end
    function folderTreeDisplayNewListCancel(src,evnt)
        gui = guidata(parentFig);
        gui.openFolderFileChooserList = 'canceled';
        guidata(parentFig,gui);
        close(get(src,'Parent'));
        close(folderTreeFig)
    end
%     guidata(parentFig,gui);

end



