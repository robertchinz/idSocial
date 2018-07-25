function [H p]=idSocial_auxiliaries_plotmode2StatSign(output_plot,statistical_test_type,statistical_test_significance,bootstrap_repetitions,bootstrap_statfunc,display_mode,axes_handle,colororder,tail,compare_rand_only)

if nargin<2 || isempty(statistical_test_type)
    statistical_test_type='ranksum';
end

if nargin<3 || isempty(statistical_test_significance)
    statistical_test_significance=0.05;
end

if nargin<4 || isempty(bootstrap_repetitions)
    bootstrap_repetitions=100;
end

if nargin<5 || isempty(bootstrap_statfunc)
    bootstrap_statfunc='Mean';
end
if nargin<6 || isempty(display_mode)
    display_mode='plot2d';
end

if nargin<7 || isempty(axes_handle)
    axes_handle=[];
end

if nargin<8 || isempty(colororder)
    colororder=get(0,'defaultAxesColorOrder');
end

if nargin<9 || isempty(tail)
    tail='both';
end

if nargin<10 || isempty(compare_rand_only)
    compare_rand_only=false;
end

no_subplots = size(output_plot.data_sign,1);
for sp=1:max(1,no_subplots)
    no_legends = size(output_plot.legendstring,2);
    no_data_points=NaN(no_legends,1);
    for legendentr=1:no_legends
        if ~strfind(lower(output_plot.statistics_type{legendentr}),'hist')
            no_data_points(legendentr)=size(output_plot.data_sign{sp,legendentr},output_plot.statistics_on_idx(sp,legendentr));
        else
%             no_data_points(legendentr)=...
%                 sum(~isnan(output_plot.data_signMedian{sp,legendentr}),output_plot.statistics_on_idx(sp,legendentr));

            no_data_points(legendentr)=size(output_plot.data_signMedian{sp,legendentr},output_plot.statistics_on_idx(sp,legendentr));
        end
    end
    max_data_points=max(no_data_points);
    if strcmpi(display_mode,'hist')
        signval_all=NaN(no_legends, max_data_points);
    else
        signval_all=NaN(no_legends, max_data_points,output_plot.xaxis_length);
    end
    for legendentr=1:no_legends
        if ~isempty(output_plot.data_sign{sp,legendentr})
            if strcmpi(display_mode,'hist')
                signval=output_plot.data_signMedian{sp,legendentr};
                try
                    signval_all(legendentr,1:no_data_points(legendentr),:)=signval;%/(xticlabel(2)-xticlabel(1))+xticlabel(1);
                catch
                    keyboard
                end
            else
                signval=output_plot.data_sign{sp,legendentr};
                sign_dims=ndims(signval);
                sign_indices=[output_plot.statistics_on_idx(sp,legendentr) output_plot.xaxis_idx setdiff(1:sign_dims,[output_plot.statistics_on_idx(sp,legendentr) output_plot.xaxis_idx])];
                sign_val_perm=squeeze(permute(signval,sign_indices));
                
                signval_all(legendentr,1:no_data_points(legendentr),:)=sign_val_perm;
            end
        end
    end
    disp('Statistical significance is being calculated...');
    
    if isfield(output_plot,'dist_edges'); edges = output_plot.dist_edges; else edges = []; end;
    if isfield(output_plot,'dist_method'); method = output_plot.dist_method; else method = []; end;
    if isfield(output_plot,'dist_bandwidth'); bandwidth = output_plot.dist_bandwidth; else bandwidth = []; end;

    [H, p]=idSocial_auxiliaries_StatSign4groups(signval_all,statistical_test_type,statistical_test_significance,axes_handle,colororder,bootstrap_repetitions,display_mode,bootstrap_statfunc,tail,compare_rand_only,...
        edges,bandwidth,method);
    disp('Done.')
end