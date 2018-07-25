function legend_handle = idSocial_auxiliaries_correctLegend(legend_handle,parent_ax,pos)

if nargin<3 
    pos = [ 13.35 3.55 5.3 .7];
elseif isempty(pos)
    pos = legend_handle.Position;
end

if ~verLessThan('matlab','8.5.0')
    set(legend_handle,'Units','centimeters');
    strng = legend_handle.String;
    if isempty(pos)
%     pos = legend_handle.Position;
    end
    new_lh = axes('Units','centimeters','Position',pos,'Box','on','Xtick',[],'YTick',[],'LineWidth',legend_handle.LineWidth);
%     set( new_lh,'Position',[ 13.5 3.5 5.3 pos(4)])

    set( new_lh,'Position',pos)

    hold on
    lns = NaN(1,numel(strng));
    txt = NaN(1,numel(strng));
    %     wdth = pos(3)/(numel(strng));
    %     yp = pos(2)+pos(4)/2;
    coOrder = get(parent_ax,'ColorOrder');
    set(new_lh,'XLim',[0 2*numel(strng)],'YLim',[0 1]);
    for k = 1:2:2*numel(strng)
        lns(k) = plot([k-1+.2 ,k],[.5 .5],'-','Color', coOrder(ceil(k/2),:),'LineWidth',legend_handle.LineWidth);
        set(lns(k),'Tag','legend')
    end
    for k = 2:2:2*numel(strng)
        txt(k) = text(k-.5,.5,strng{k/2},'VerticalAlignment','middle','HorizontalAlignment','center','units','data');
    end
    %     set(new_lh,'Position',pos)
    delete(legend_handle);
    %     legend_handleOld = legend_handle;
    legend_handle = new_lh;
else
    lh_pos =  get(legend_handle,'Position');
    set(legend_handle,'Position',[2.25 lh_pos(2) lh_pos(3) lh_pos(4)]);
    
    axis_width_plotUnits = diff(get(parent_ax,'XLim'));
    axis_width_cm = get(parent_ax,'Position');
    axis_width_cm = axis_width_cm(3);
    cm2plotWidth = axis_width_cm/axis_width_plotUnits;
    
    ax_x = get(parent_ax,'Position');
    ax_x = ax_x(1);
    
    lh_x_plotUnits = 4.6;
    lh_x2_plotUnits = 19.5;
    x_legend = ax_x(1) + lh_x_plotUnits*cm2plotWidth;
    width_legend = ax_x(1) + lh_x2_plotUnits*cm2plotWidth - x_legend;
    set(legend_handle,'Position',[x_legend lh_pos(2) width_legend lh_pos(4)]);
    
    % Find all lines which have no empty tags
    lh_pos = get(legend_handle,'Position');
    linesInPlot = findobj(legend_handle,'type','line','-and','-regexp','Tag','[^'']');
    textInPlot = findobj(legend_handle,'type','text');
    xd = NaN(size(linesInPlot,1),2);
    xt = NaN(size(linesInPlot,1),3);
    for k = 1:size(linesInPlot,1)
        xd(k,:) = get(linesInPlot(k),'XData');
        xt(k,:) = get(textInPlot(k),'Position');
    end
    minx = min(xd(:,1));
    len = diff(xd(1,:));
    new_len = len*.75;
    gaps = new_len + (xd(1,1) - xd(2,2));
    gap_text = xt(end,1) - xd(end,2);
    x_act = minx;
    for k = size(linesInPlot,1):-1:1
        set(linesInPlot(k),'XData',[x_act x_act+new_len]);
        set(textInPlot(k),'Position',[x_act + new_len + gap_text xt(k,2) xt(k,3)],'Units','data');
        text_extend = get(textInPlot(k),'Extent');
        x_act = x_act + new_len +  text_extend(3) + 2*gap_text;
    end
end