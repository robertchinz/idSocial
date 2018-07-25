function [fh,options] = idSocial_auxiliaries_prepareFigure(figops)

if nargin < 1 || isempty(figops)
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
%     c = parula;
    def_options.color_order = get(groot,'defaultAxesColorOrder');
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
def_options.between_axes_margin_vertical = [];
def_options.dinA4size = dinA4size;
def_options.part_height = (def_options.dinA4size(2)-def_options.bottom_margin-def_options.top_margin)/def_options.no_parts;
def_options.plot_area_width=def_options.dinA4size(1)-def_options.right_margin-def_options.left_margin;

def_options.axes_linewidth = 2;
def_options.axes_fontsize = 14;

[~, options]=idSocial_readparams([],figops,def_options);

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
dinA4size = options.dinA4size;
part_height=(dinA4size(2)-bottom_margin-top_margin)/no_parts;%*1-between_axes_margin;
plot_area_width=dinA4size(1)-right_margin-left_margin;
options.part_height=part_height;
options.plot_area_width=plot_area_width;

axes_linewidth = options.axes_linewidth;
axes_fontsize = options.axes_fontsize;
%%
set(0,'DefaultPatchLineSmoothing','On')
set(0,'DefaultLineLineSmoothing','On')

screensize_pxl=get(0,'screensize');
screensize_cm=[47 29.7];
figure_position=[screensize_cm(1)-24 0 dinA4size(1) top_margin+bottom_margin];
fh=figure('Unit','centimeters','Position',figure_position,'Color','w','Renderer','opengl');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [47 29.7]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 47 29.7]);
set(gcf, 'defaultAxesColorOrder',color_order);
% set(fh,'MenuBar', 'None');
set(fh,'Userdata',options)
set(gca,'Visible','off')
% set(gca,'Linewidth',axes_linewidth,'FontSize',axes_fontsize);
idSocial_auxiliaries_addPlotRow(fh,options);
