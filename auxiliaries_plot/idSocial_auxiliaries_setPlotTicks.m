function idSocial_auxiliaries_setPlotTicks(ax,vec_size,edges,wanted_positions)

if nargin<4 || isempty(wanted_positions)
    wanted_positions = cell(numel(edges),1);
    for k=1:numel(edges)
        wanted_positions{k} = unique(round(edges{k}));
    end
end

if numel(wanted_positions) == 1
    wanted_positions = vertcat(wanted_positions,wanted_positions);
end

siz = vec_size;
plot_edges = cell(1,numel(edges));
xtick_new = cell(1,numel(edges));
xtickpos_new = cell(1,numel(edges));
for k=1:numel(edges)
    plot_edges{k} = unique(round(edges{k}));%-4:2:4;
    xtick_new{k}=edges{k}(1):(edges{k}(end)-edges{k}(1))/(siz(k)-1):edges{k}(end);
    xtickpos_new{k} = interp1(xtick_new{k},1:siz(k),plot_edges{k});
end
labelsX = cell(size(plot_edges{2},2),1);
for k = 1:size(plot_edges{2},2)
    labelsX{k} = num2str(plot_edges{2}(k));
end
labelsY = cell(size(plot_edges{1},2),1);
for k = 1:size(plot_edges{1},2)
    labelsY{k} = num2str(plot_edges{1}(k));
end
set(ax,'XTick',xtickpos_new{1},'XTickLabel',labelsX,'YTick',xtickpos_new{2},'YTickLabel',labelsY)
