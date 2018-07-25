function output_plot_merged = idSocial_auxiliaries_mergeOutputPlot(input_data,func1,func2,func1_legid,func2_legid)
% func='dynamics_TurningBalance';
% func1='dynamics_TurningBalanceAway';
% func2='dynamics_TurningBalanceTowards';
% func1_legid='Away';
% func2_legid='Towards';

pm=input_data(1,1,1).(func1).output_plot;
pmrand=input_data(1,1,1).(func2).output_plot;
% Merge plot_mode
%     no_entries=size(pmrand.data_string,2);
    pmrand.data_string=...
        cellfun(@(x) strrep(x,x,[char(x(1)+1) x(2:end)]),pmrand.data_string,'UniformOutput',false);
   
    pmrand.legendstring=...
        cellfun(@(x) strrep(x,x,[char(x(1)+1) x(2:end)]),pmrand.legendstring,'UniformOutput',false);
   
    pmrand.legendstring=...
        cellfun(@(x) strrep(x,x,[x ', ' func2_legid]),pmrand.legendstring,'UniformOutput',false);
    
    pm.legendstring=...
        cellfun(@(x) strrep(x,x,[x ', ' func1_legid]),pm.legendstring,'UniformOutput',false);
    
    
    pm.data=[pm.data pmrand.data];
    pm.data_sign=[pm.data_sign pmrand.data_sign];
    pm.data_signMedian=[pm.data_signMedian pmrand.data_signMedian];
    pm.data_Median=[pm.data_Median pmrand.data_Median];
    pm.data_dev=[pm.data_dev pmrand.data_dev];
    pm.statistics_type=[pm.statistics_type pmrand.statistics_type];
    pm.statistics_on_idx=[pm.statistics_on_idx pmrand.statistics_on_idx];
    pm.data_no_datapoints=[pm.data_no_datapoints pmrand.data_no_datapoints];
    pm.data_string=[pm.data_string pmrand.data_string];
    pm.filterstring=[pm.filterstring pmrand.filterstring];
    pm.legendstring=[pm.legendstring pmrand.legendstring];
    pm.dim_names=[pm.dim_names pmrand.dim_names];
    pm.subplotstring=[pm.subplotstring pmrand.subplotstring];
    output_plot_merged=pm;
%     def_options.act_method='dynamics_TurningBalance';