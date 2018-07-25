function plot_mode_def=idSocial_auxiliaries_makeDefPlotMode(plot_mode)
 
if nargin<1
   plot_mode = [];
end

dimnames = {'Group'; 'Subset'; 'Trial'; 'Time'; 'Focal'; 'Neighbor';  'EdgeX'; 'EdgeY'; 'EdgeZ'};
opnames1 = {'Mean','Median','Sum','Mode','Hist','Positive Ratio'};
opnames2 = {'Mean','Median'};

plot_mode_def.xaxis = dimnames(1:6);
plot_mode_def.yaxis = [];

plot_mode_def.statistics = opnames1;
plot_mode_def.data = opnames2;
plot_mode_def.display_mode = 'plot2d';
plot_mode_def.normalization = 'none';
plot_mode_def.extraDims = [];

plot_mode_def.xlabel = '';
plot_mode_def.ylabel = '';


if ~isempty(plot_mode)
fnames = fieldnames(plot_mode);
fnames_def = fieldnames(plot_mode_def);
for k = 1:size(fnames,1)
    idx = strcmpi(fnames{k},fnames_def);
    if any(idx)
        plot_mode_def.(fnames{k}) = plot_mode.(fnames{k});
    end
end
end

plot_mode_def2.yaxis=[];
fnames_def = fieldnames(plot_mode_def);
fnames_def2 = fieldnames(plot_mode_def2);
for k = 1:size(fnames_def2,1)
    idx = strcmpi(fnames_def2{k},fnames_def);
    if ~any(idx)
        plot_mode_def.(fnames_def2{k}) = [];
    end
end
