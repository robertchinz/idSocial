function help_comment = idSocial_auxiliaries_extractHelp(file_act) 
if nargin <1
    file_act = 'C:\Users\Robert\Syncplicity Folders\Syncplicity\git_idSocial\core\idSocial_dynamicMapsPolar.m';
end
tempString= [];
if exist(file_act,'file')==2
    
    fileID = fopen(file_act,'r');
    tempString = textscan(fileID,'%s','Delimiter','');
    fclose(fileID);
    
end
if ~isempty(tempString)
    tempString = tempString{1};
    
%     comment_lines = cellfun(@(x) strcmp(x(1),'%') ,tempString);
%     comment_lines = cellfun(@(x) ~isempty(regexp(x,'^(%[^%%]*)')),tempString);
%      comment_lines = cellfun(@(x) ~isempty(x),regexp(tempString,'^(%[^%{2}])'));
     comment_lines = cellfun(@(x) ~isempty(x),regexp(tempString,'^%')) & ...
         cellfun(@(x) isempty(x),regexp(tempString,'^%{2}'));


    conn = bwconncomp(comment_lines);
    
    if ~isempty(conn) && conn.NumObjects>0
        comment_lines = conn.PixelIdxList{1};
        
        
        help_comment = tempString(comment_lines,:);
        help_comment = cellfun(@(x) x(2:end),help_comment,'UniformOutput',false);
    end
else
    help_comment = 'No help available.';
end

[~,file_str] = fileparts(file_act);
% file_act = strrep(file_act,'.m',''); fileparts
colwidth = 560;
fh = figure('Units','pixels','Position',[400 200 colwidth 600],'Name',[file_str ' Help'], 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off','Resize','off');

new_text = strcat(help_comment,'\n');
new_text = sprintf([new_text{:}]);
panelH = uipanel(fh,'Units','normalized',...
    'Position', [0.01 0.01 .98 .98],'BackgroundColor',[1 1 1]);
set(panelH,'Units','pixels');
pan1PosPxl =  get(panelH,'Position');
panel2H = uipanel(panelH,'Units','pixels',...
    'Position', [0 -100 pan1PosPxl(3) pan1PosPxl(4)+100],'BackgroundColor',[1 1 1]);



textBox = uicontrol(panel2H,'Units','normalized','Position',[0.02 0.02 0.96 0.96],'Style','text','String',help_comment,'HorizontalAlignment','left','FontSize',12,'BackgroundColor','w');
set(textBox,'Units','pixels');
margin = 5;
textExt = get(textBox,'Extent');
% textExt(3:4) = textExt(3:4)+margin;
set(panel2H,'Position', [0 pan1PosPxl(4)-textExt(4)-margin pan1PosPxl(3)+2*margin textExt(4)+2*margin])
% pan2PosPxl =  get(panel2H,'Position');
set(textBox,'Units','pixels','Position', [margin 0 textExt(3) textExt(4)])

set(textBox,'Units','normalized');
yrange = min((pan1PosPxl(4)-textExt(4)-margin),(-pan1PosPxl(4)+textExt(4)-margin));

if textExt(4)+margin > pan1PosPxl(4)
sliderH = uicontrol('Style','Slider','Parent',panelH,...
    'Units','normalized','Position',[0.96 0 0.04 1],...
    'Value',1,'Callback',{@slider_callback1,panel2H});

set(gcf,'WindowScrollWheelFcn', {@mouseScroll,panel2H,sliderH});
end

    function cellselect(src,evnt)
        %         old_data = get(src,'Data');
        %         set(src,'Data',{'_'})
        %         set(src,'Data',old_data)
        
        old_data = get(src,'String');
        set(src,'String',{'_'})
        set(src,'String',old_data)
        
    end


    function slider_callback1(src,eventdata,arg1)
        
        sliderMin=get(src,'Min');
        sliderMax=get(src,'Max');
        dummy=get(src,'Value');
        s=dummy/(sliderMax-sliderMin);
        dummy=get(arg1,'position');frameWidth=dummy(3);
        set(arg1,'position',[dummy(1) (s)*yrange dummy(3) dummy(4)])
        
    end
    function mouseScroll(src,eventdata, handles,arg1)
        set(handles,'Units','pixel')
        scrollHeight = 40;
        act_pos = get(handles,'Position');
            val=0;
            if eventdata.VerticalScrollCount > 0 && act_pos(2) <0
                val = scrollHeight;
                act_pos(2) = min(act_pos(2)+val,0);
            elseif eventdata.VerticalScrollCount < 0 && act_pos(2) > yrange;%-(rowHeightCm*no_rows-panel2posCM(4))
                val = -scrollHeight;
                act_pos(2) = max(act_pos(2)+val,yrange);
            end
        
            
            set(handles,'Position',[act_pos(1)  act_pos(2) act_pos(3) act_pos(4)])

      
            set(arg1,'Value',abs(act_pos(2))/-yrange)

    end

end