function [fh, textBox,textBoxAxes] = idSocial_auxiliaries_message(help_comment,textBox,textBoxAxes) 

if nargin<2
    textBox=[];
end
if nargin<3
    textBoxAxes=[];
end

if nargin<1
    help_comment = '';
    colwidth = 560;
    fh = figure('Units','pixels','Position',[400 200 colwidth 100],'Name','Progress', 'ToolBar', 'none','MenuBar', 'none','NumberTitle','off','Resize','off');
    set(fh,'Tag','idSocialMessageBoxFigure');
    % new_text = strcat(help_comment,'\n');
    % new_text = sprintf([new_text{:}]);
    panelH = uipanel(fh,'Units','normalized',...
        'Position', [0.01 0.01 .98 .98],'BackgroundColor',[1 1 1]);
    set(panelH,'Units','pixels');
    pan1PosPxl =  get(panelH,'Position');
    panel2H = uipanel(panelH,'Units','pixels',...
        'Position', [0 -100 pan1PosPxl(3) pan1PosPxl(4)+100],'BackgroundColor',[1 1 1]);
    
    
    
%     textBox = uicontrol(panel2H,'Units','normalized','Position',[0.02 0.02 0.96 0.96],'Style','text', ...
%         'FontName', 'Consolas','String',help_comment,'HorizontalAlignment','left','FontSize',12,'BackgroundColor','w');
    textBoxAxes = axes('Units','normalized','Position',[0.02 0.02 0.96 0.96], 'Parent', panel2H,'XLim',[0 1],'YLim',[0 1]);
   axis off

    textBox = text(0,0,'Test','FontSize',12,'BackgroundColor','w','HorizontalAlignment','left','VerticalAlignment','bottom');
    set(textBoxAxes,'Tag','idSocialMessageBox');
    set(textBoxAxes,'Units','pixels');
    margin = 0;
%     textExt = get(textBox,'Extent');
     set(textBox,'Units','Pixels');

    textExt = get(textBox,'Extent');
    % textExt(3:4) = textExt(3:4)+margin;
%     set(panel2H,'Position', [0 pan1PosPxl(4)-textExt(4) pan1PosPxl(3)+2*margin textExt(4)])
%     set(panel2H,'Position', [0 pan1PosPxl(4)-textExt(4) textExt(3)+2*margin textExt(4)])

    % pan2PosPxl =  get(panel2H,'Position');
    set(textBoxAxes,'Units','pixels','Position', [margin 0 textExt(3) textExt(4)])
    
    set(textBox,'Units','normalized');
%     yrange = min((pan1PosPxl(4)-textExt(4)-margin),(-pan1PosPxl(4)+textExt(4)-margin));
%     yrange = -pan1PosPxl(4)+textExt(4)-margin;% min((pan1PosPxl(4)-textExt(4)-margin),(-pan1PosPxl(4)+textExt(4)-margin));
    yrange = -(-pan1PosPxl(4)+textExt(4)+margin);% min((pan1PosPxl(4)-textExt(4)-margin),(-pan1PosPxl(4)+textExt(4)-margin));
    

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
            actString =  vertcat(actString,{help_comment});
            
        else
            actString = vertcat({actString},{help_comment});
            
        end
        set(textBox,'String',actString);
        
    catch
        keyboard
    end
    set(textBox,'Units','pixels');
    margin = 0;
    marginX = 15;
    
    corr_factor = .858; % Arial?
%     corr_factor = 1;%.95; % Consolas?
    
    textExt = get(textBox,'Extent');
    textExt(3) = textExt(3) + 10;
    textExt(4) = textExt(4)*corr_factor;
    
    %     set(textBoxAxes,'Units','pixels','Position', [margin 0 textExt(3) textExt(4)])
    
    
    
    figpos = get(gcf,'Position');
    set(gcf,'Position', [figpos(1) figpos(2) (textExt(3)+2*marginX)*(1+0.0) figpos(4)]);
    textExt = get(textBox,'Extent');
    textExt(4) = textExt(4)*corr_factor;
    
    panelH = get(get(get(textBox,'Parent'),'Parent'),'Parent');
    panel2H = get(textBoxAxes,'Parent');
    set(panelH ,'Units','normalized','Position', [0.01 0.01 .99 .99]);
    set(panelH,'Units','pixels');
    pan1PosPxl =  get(panelH,'Position');
    
%     set(panel2H,'Units','pixels','Position', [margin 0 textExt(3) textExt(4)])

    set(panel2H,'Position', [0 pan1PosPxl(4)-textExt(4) pan1PosPxl(3)+2*marginX textExt(4)])
    set(textBoxAxes,'Units','normalized','Position', [0 0 1 1])
    set(textBoxAxes,'Units','pixels')

    yrange = -(-pan1PosPxl(4)+textExt(4)+margin);
    
    if textExt(4)+margin > pan1PosPxl(4)
        set(gcf,'Position', [figpos(1) figpos(2) (textExt(3)+2*marginX)*(1+0.04) figpos(4)]);
        set(panelH ,'Units','normalized','Position', [0.01 0.01 .99 .99]);
%               set(panel2H ,'Units','normalized','Position', [0.01 0.01 .99 .99]);


        set(panelH,'Units','pixels');
        panelHpos = get(panelH,'Position');
        set(panelH,'Units','normalized');



        set(panel2H,'Units','pixels');
        actpos = get(panel2H,'Position');
                set(panel2H,'Position',[1 -textExt(4)+panelHpos(4) actpos(3) textExt(4)])

%         set(panel2H,'Position',[actpos(1) 0 actpos(3) textExt(4)])
%         set(textBox,'Units','pixels','Position', [margin -40 textExt(3) textExt(4)])
        set(textBoxAxes,'Units','normalized','Position', [0 0 1 1])
        set(textBoxAxes,'Units','pixels')
        
        sliderH = uicontrol('Style','Slider','Parent',panelH,...
            'Units','normalized','Position',[0.96 0 0.04 1],...
            'Value',0,'Callback',{@slider_callback1,panel2H});
        
        set(gcf,'WindowScrollWheelFcn', {@mouseScroll,panel2H,sliderH});
    end
    drawnow
end
disp(help_comment)



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
        dummy=get(arg1,'position');
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