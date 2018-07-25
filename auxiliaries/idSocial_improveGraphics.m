function ax_out = idSocial_improveGraphics(fg)

fontsize_ax = 16;
linewidth_ax = 2;
fontsize_title = 20;
fontsize_labels = 20;

set(0,'DefaultPatchLineSmoothing','On')
set(0,'DefaultLineLineSmoothing','On')
% scrSize = get(0, 'MonitorPositions');
units_orig = get(fg,'Units');

set(fg, 'Units','normalized','Position',[.05 .05 .7 .7])
set(fg,'Color','w','Renderer','opengl')
% set(fg,'MenuBar', 'None');

set(fg,'Units',units_orig)


% Axes
ax = findobj(get(fg,'Children'),'Type','axes','-and','Tag','');
set(ax,'FontSize',fontsize_ax,'LineWidth',linewidth_ax)

% Axes limits
ax_ch=allchild(ax);
xdata_cell = cell(1,size(ax_ch,1));
try
for k=1:size(ax_ch,1)
    act_ch = get(ax_ch(k));
    if isfield(act_ch,'XData')
        xdata_cell{k} = get(ax_ch(k),'XData');
    end
end
all_xdata = cellfun(@(x) x(:),xdata_cell(:),'UniformOutput',false);
all_xdata = vertcat(all_xdata{:});
set(ax,'XLim',[min(all_xdata) max(all_xdata)]);
catch
end

xlim=get(gca,'XLim');
axpos = get(gca,'Position');
dat2heightRatio = diff(xlim)/axpos(4);
umargin = .05; % margin in cm
umargin2dat = umargin * dat2heightRatio;
set(gca,'XLim',[xlim(1)-umargin2dat xlim(2)+umargin2dat])

% Title
th=get(ax,'Title');
for k = 1:size(th,1)
    try
        if iscell(th)
            for k2 = 1:numel(th)
                set(th{k2},'FontSize',fontsize_title);
            end
        else
            set(th(1),'FontSize',fontsize_title);
            
        end
    catch
        keyboard
    end
end

% Axes Labels
haxes = findall(fg,'type','axes');
h_axLabels = get(haxes,{'XLabel' 'YLabel'});
set([h_axLabels{:}],'FontSize',fontsize_labels)

% box(ax,'off');
% axes('Position',get(ax,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',get(ax,'LineWidth'),'Tag','box');
% uistack(gca,'bottom')

% if nargin>1 && ~isempty(save_dir)
% %     root_dir = textscan('%s',save_dir,'Delimiter','\\')
%     export_fig(save_dir)
% end
