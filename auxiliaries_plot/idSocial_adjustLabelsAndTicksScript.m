function idSocial_adjustLabelsAndTicksScript(ax,xshift,yshift,xshift_label,yshift_label)
if nargin<1 || isempty(ax)
    ax=gca;
end
if nargin<2 || isempty(xshift)
    xshift = -.1;
end
if nargin<3 || isempty(yshift)
    yshift = -.01;
end
if nargin<4 || isempty(xshift_label)
    xshift_label = -.5;
end
if nargin<5 || isempty(yshift_label)
    yshift_label = 0;
end
% Put Ticks and Label manually (X)
label_list = get(ax,'XTickLabel');
label_pos = get(ax,'XTick');
set(ax,'XTickLabel','')
yl = get(ax,'YLim');

for k=1:numel(label_list)
    th=text(label_pos(k),yl(1)+yshift,label_list{k},'Fontsize',get(ax,'FontSize'),'VerticalAlignment','top','HorizontalAlignment','center','Units','data');
end
texExt = get(th,'Extent');
xlh=get(ax,'xlabel');
% xlabel('Surface Diameter (cm)','FontSize',fontsize_label);
% set(xlh,
xlpos=get(xlh,'Position');
set(xlh,'Position',[xlpos(1) yl(1)+yshift-texExt(4)+yshift_label])

% Put Ticks and Label manually (Y)
label_list = get(ax,'YTickLabel');
label_pos = get(ax,'YTick');
set(ax,'YTickLabel','')
xl = get(ax,'XLim');
% xshift = -.1;
for k=1:numel(label_list)
    th=text(xl(1)+xshift,label_pos(k),label_list{k},'Fontsize',get(ax,'FontSize'),'VerticalAlignment','middle','HorizontalAlignment','right','Units','data');
end
texExt = get(th,'Extent');
ylh=get(ax,'ylabel');
% xlabel('Surface Diameter (cm)','FontSize',fontsize_label);
% set(xlh,
ylpos=get(ylh,'Position');
set(ylh,'Position',[xl(1)+xshift-texExt(3)+xshift_label ylpos(2)])
