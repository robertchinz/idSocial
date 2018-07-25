function [plot_mode_idx, dim_names_act]=idSocial_function_wrapper_idxcell2combinations(plot_mode,cell_size,dim_names)
% 

if size(plot_mode,2)<2 % No indices given, size(plot_mode)=1; calculate possible combinations
    no_lgs=1;
    if ~isempty(plot_mode{1})
        no_lgs_sep=NaN(size(plot_mode{1}));
    else
        no_lgs_sep=1;
    end
    for lg=1:size(plot_mode{1},2)
        no_lgs=no_lgs*cell_size(plot_mode{1}(lg)); % {1} because if size(legend_idx,1)==1, size(yval) also ==1.
        no_lgs_sep(lg)=cell_size(plot_mode{1}(lg));
    end
    if ~isempty(no_lgs_sep)
        sets=cell(size(no_lgs_sep,1),1);
        for k=1:size(no_lgs_sep,2)
            sets{k}=1:no_lgs_sep(k);
        end
        c = cell(1, numel(sets));
        [c{:}] = ndgrid( sets{end:-1:1} );
        lg_combi = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
        lg_combi = lg_combi(:,end:-1:1);
        lg_combi2 = cell(size(lg_combi,1),1);
        for lt=1:size(lg_combi,1)
            lg_combi2{lt}=lg_combi(lt,:);
        end
        legend_idx=repmat({plot_mode{1}},size(lg_combi2,1),1);
        legend_idx(:,2)=lg_combi2;
        plot_mode_idx=legend_idx;
        dim_names_act=repmat(dim_names,size(plot_mode_idx,1),1);
    end

else
   plot_mode_idx=plot_mode;
   dim_names_act=dim_names;
end

