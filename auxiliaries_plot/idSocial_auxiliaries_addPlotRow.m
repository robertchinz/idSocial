function [ax,fh,options] = idSocial_auxiliaries_addPlotRow(fh,figops)

if nargin < 1 || isempty(fh)
    fh = gcf;
end

if nargin < 2 || isempty(figops)
    figops = [];
end
dinA4size=[210 297]/10;

if verLessThan('matlab','8.4.0')
    
    color_order = [ 0    0.4470    0.7410;
        0.8500    0.3250    0.0980;
        0.9290    0.6940    0.1250;
        0.4940    0.1840    0.5560;
        0.4660    0.6740    0.1880;
        0.3010    0.7450    0.9330;
        0.6350    0.0780    0.1840];
    n=16;
    
    def_options.color_order = [interp1(color_order(:,1), 1:(7-1)/(n-1):7)' ...
        interp1(color_order(:,2), 1:(7-1)/(n-1):7)' ...
        interp1(color_order(:,3), 1:(7-1)/(n-1):7)'];
else
    def_options.color_order = parula;
end

def_options.no_parts = 3;
def_options.act_row = 1;
def_options.act_col = 1;
def_options.act_noCol = 1;
def_options.left_margin=2;
def_options.right_margin=2;
def_options.top_margin=1;
def_options.bottom_margin=2.5;
def_options.between_axes_margin=2;
def_options.between_axes_margin_vertical=[];
def_options.dinA4size = dinA4size;
def_options.part_height = (def_options.dinA4size(2)-def_options.bottom_margin-def_options.top_margin)/def_options.no_parts;
def_options.plot_area_width=def_options.dinA4size(1)-def_options.right_margin-def_options.left_margin;

def_options.axes_linewidth = 2;
def_options.axes_fontsize = 14;



if nargin < 2 || isempty(figops)
    ud = get(fh,'Userdata');
    if isstruct(ud)
        [~, options]=idSocial_readparams([],ud,def_options);
    end
end
if isempty(figops)
    figops = def_options;
end
if ~isempty(figops) && isstruct(figops)
    [~, options]=idSocial_readparams([],figops,def_options);
    ud = get(fh,'Userdata');
    if isstruct(ud)
        options.act_row = ud.act_row;
        options.act_col = ud.act_col;
    end
end

if isempty(options.between_axes_margin_vertical)
    options.between_axes_margin_vertical = options.between_axes_margin;
end

color_order = options.color_order;
no_parts = options.no_parts;
left_margin=options.left_margin;
right_margin=options.right_margin;
top_margin=options.top_margin;
bottom_margin=options.bottom_margin;
between_axes_margin=options.between_axes_margin;
between_axes_margin_vertical=options.between_axes_margin_vertical;
part_height = options.part_height;
dinA4size = options.dinA4size;
% part_height=(dinA4size(2)-bottom_margin-top_margin)/no_parts;%*1-between_axes_margin;
plot_area_width=dinA4size(1)-right_margin-left_margin;
options.part_height=part_height;
options.plot_area_width=plot_area_width;

axes_linewidth = options.axes_linewidth;
axes_fontsize = options.axes_fontsize;
%%
if options.act_row > 1
    between_margin = between_axes_margin_vertical;
else
    between_margin = 0;
end
set(fh,'Units','centimeters');
figure_position=get(fh,'Position');
act_yposition=figure_position(4);
set(fh,'Position',[figure_position(1:3) max(act_yposition,act_yposition+part_height+between_margin)])



% Move previous axes further up
allAxesInFigure = findall(fh,'type','axes');

no_axes=numel(allAxesInFigure);
for act_axes=1:no_axes
    act_position=get(allAxesInFigure(act_axes),'Position');
    set(allAxesInFigure(act_axes),'Position',...
        [act_position(1) act_position(2)+part_height+between_margin act_position(3:4)]);
end
options.act_row = options.act_row + 1;
options.act_col = 1;

if ~verLessThan('matlab','8.5.0')
    legHandles = findall(fh,'Type','legend');
    no_axes=numel(legHandles);
    for act_axes=1:no_axes
        act_position=get(legHandles(act_axes),'Position');
        set(legHandles(act_axes),'Position',...
            [act_position(1) act_position(2)+part_height+between_margin act_position(3:4)]);
    end
    
    cbHandles = findall(fh,'Type','colorbar');
    no_axes=numel(cbHandles);
    for act_axes=1:no_axes
        act_position=get(legHandles(act_axes),'Position');
        set(cbHandles(act_axes),'Position',...
            [act_position(1) act_position(2)+part_height+between_margin act_position(3:4)]);
    end
end


% ax=axes('Unit','centimeters','Position',[left_margin   bottom_margin...
%     plot_area_width part_height]);
set(gca,'Linewidth',axes_linewidth,'FontSize',axes_fontsize);
ax=gca;

set(fh,'Userdata',options)