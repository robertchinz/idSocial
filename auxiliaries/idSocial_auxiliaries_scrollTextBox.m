function textBox = idSocial_auxiliaries_scrollTextBox(help_comment,textBox) 

if nargin<2
    textBox=[];
end

if nargin<1
    help_comment = '';
    colwidth = 560;
    fh = figure('Units','pixels','Position',[400 200 colwidth 100],'Name','Progress', 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off','Resize','off');
    set(fh,'Tag','idSocialMessageBox');
    % new_text = strcat(help_comment,'\n');
    % new_text = sprintf([new_text{:}]);
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
elseif ~isempty(textBox) && ishghandle(textBox)
    try
        actString = get(textBox,'String');
        if iscell(actString)
            set(textBox,'String',vertcat(actString,{help_comment}));
        else
            set(textBox,'String',vertcat({actString},{help_comment}));
        end
    catch
        keyboard
    end
    set(textBox,'Units','pixels');
    margin = 5;
    textExt = get(textBox,'Extent');
    pan1PosPxl =  get(get(get(textBox,'Parent'),'Parent'),'Position');
    % textExt(3:4) = textExt(3:4)+margin;
    set(get(textBox,'Parent'),'Position', [0 pan1PosPxl(4)-textExt(4)-margin textExt(3)+2*margin textExt(4)+2*margin])
    % pan2PosPxl =  get(panel2H,'Position');
    set(textBox,'Units','pixels','Position', [margin 0 textExt(3) textExt(4)])
    
    set(textBox,'Units','normalized');
    yrange = min((pan1PosPxl(4)-textExt(4)-margin),(-pan1PosPxl(4)+textExt(4)-margin));
    
    if textExt(4)+margin > pan1PosPxl(4)
        actpos = get(get(textBox,'Parent'),'Position');
        set(get(textBox,'Parent'),'Position',[actpos(1) 0 actpos(3) actpos(4)])
        sliderH = uicontrol('Style','Slider','Parent',get(get(textBox,'Parent'),'Parent'),...
            'Units','normalized','Position',[0.96 0 0.04 1],...
            'Value',0,'Callback',{@slider_callback1,get(textBox,'Parent')});
        
        set(gcf,'WindowScrollWheelFcn', {@mouseScroll,get(textBox,'Parent'),sliderH});
    end
end
disp(help_comment)

% file_act = strrep(file_act,'.m',''); fileparts


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