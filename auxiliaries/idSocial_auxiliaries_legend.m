idSocial_auxiliaries_legend(axes_handle,ax_plots,options)


% Manual legend:
leg_pos = [7 2];
fontsize_legend = 7;
linelength_legend = 1;
linewidth_legend=linewidth_axes;
horzspacing_legend = .5;
vertspacing_legend = 1;
legend_margin = .2;
string_legend = {'Data'; 'Control'};
axes_handle = modevsage_axes;
ax_plots = [data_mean_day data_mean_day_rand];

axes(modevsage_axes)
color_legend = get(axes_handle,'Color');
no_leg_entries = numel(string_legend);
legend_gridx = [leg_pos(1) leg_pos(1)+linelength_legend+horzspacing_legend];
legend_gridy = leg_pos(2):vertspacing_legend:leg_pos(2)+(no_leg_entries-1)*vertspacing_legend;
% legend_gridy = legend_gridy(end:-1:1);
% 

string_legend = string_legend(end:-1:1);
ax_plots = ax_plots(end:-1:1);
% units_orig = get(axes_handle,'Units');
% set(axes_handle,'Units','centimeter');
% pos_cm = get(axes_handle,'Position');
% set(axes_handle,'Units','pixel');
% pos_pxl = get(axes_handle,'Position');
% set(axes_handle,'Units',units_orig);

th = NaN(no_leg_entries);
max_extentX = -inf;
max_extentY = -inf;

for act_entry = 1:no_leg_entries
    act_linewidth = get(ax_plots(act_entry),'LineWidth');
    act_color = get(ax_plots(act_entry),'Color');
    act_linestyle = get(ax_plots(act_entry),'LineStyle');
    plot([legend_gridx(1) legend_gridx(1)+linelength_legend], [legend_gridy(act_entry) legend_gridy(act_entry)], ...
        'Color',act_color,'LineWidth', act_linewidth,'LineStyle',act_linestyle);
    th(act_entry) = text(legend_gridx(2),legend_gridy(act_entry),string_legend{act_entry},'FontSize',fontsize_legend,'Units','data','VerticalAlignment','middle');
    ext = get(th(act_entry),'Extent');
    max_extentX = max(max_extentX,ext(3));
    max_extentY = max(max_extentY,ext(4));

end