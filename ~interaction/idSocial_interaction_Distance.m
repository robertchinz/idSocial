function input_data=idSocial_interaction_Distance(input_data,options,plot_mode)
% Calculates the inter-individual distance between animals of a group.
%
% Calling function from command line or script:
% The input parameters follow the same structure for all idSocial 
% functions:
%
% input_data: Matlab structure which is created using idSocial_loadData. 
%             It contains the trajectory locations, additional information
%             and the results from previous analyses. Results in input_data
%             are updated if the same variable input_data is also used as 
%             the output parameter (e.g., input_data =
%             idSocial_interaction_Distance(input_data, ...).
%
% options:    Matlab structure with fields corresponding to function 
%             specific options and additional options and filters. A list 
%             of the available options can be obtained by calling this 
%             function without input parameters, for example 
%                   options = idSocial_interaction_Distance;
%             options(1) contains options values, options(2) a short 
%             description of each options. 


if nargin < 3 || isempty(plot_mode)
    plot_mode = [];
end

%% Default options

def_options = idSocial_auxiliaries_createDefOptions(true);
def_options(1).act_method=strrep(mfilename,'idSocial_','');
def_options(2).act_method='';
def_options(1).limits_min=-inf;
def_options(2).limits_min='Min. neighbor distance (in BL)';
def_options(1).limits_max=inf;
def_options(2).limits_max='Max. neighbor distance (in BL)';
def_options(1).binwidth=inf;
def_options(2).binwidth='Binwidth (in BL)';
def_options(1).method={'hist','ksdensity_gauss','ksdensity_epanechnikov','ksdensity_triangular'};
def_options(2).method='Chose between histogram and kernel smoothing function';
def_options(1).kds_bandwidth=1;
def_options(2).kds_bandwidth='Bandwidth for ksdensity';


if nargin == 1 % Input: options. Output: plot_mode
    opt = input_data;
elseif nargin > 1
    opt = options;
end
if nargin >= 1
    if (isfield(opt(1),'edges') && ~isempty(opt(1).edges) && all(isfinite(opt(1).edges))) || (isfield(opt(1),'limits_min') && ~isempty(opt(1).limits_min) && isfinite(opt(1).limits_min) && ...
            isfield(opt(1),'limits_max') &&  ~isempty(opt(1).limits_max) && isfinite(opt(1).limits_max) && ...
            isfield(opt(1),'binwidth') && ~isempty(opt(1).binwidth) && isfinite(opt(1).binwidth));
        plot_mode_def.statistics = {'Hist+Mean','Hist+Median','Pool'};
        plot_mode_def.xaxis={'EdgeX'};
        plot_mode_def.display_mode='hist';
        plot_mode_def.normalization={'None','No_Frames','No_Frames_Effective','Density','Radial_Density','Density_Effective','Radial_Density_Effective'};
        plot_mode_def.xlabel = 'Distance';
        plot_mode_def.ylabel = 'Counts';
        plot_mode_def.extraDims = {'Focal','Neighbor','Distance'};

else
        plot_mode_def.statistics = {'Mean','Median','Pool'};
        plot_mode_def.display_mode='plot2d';
        def_options(1).act_method='interaction_Distance4Statistics';
        plot_mode_def.label = 'Distance';
        plot_mode_def.extraDims = {'Focal','Neighbor'};

    end
end

if nargin == 1
    input_data = idSocial_auxiliaries_makeDefPlotMode(plot_mode_def);
    return;
end
if nargin < 1 % No input. Output: Def. options.
    input_data = def_options;
    return;
end

[~, options_new]=idSocial_readparams(input_data,options,def_options,def_options.act_method);

if (isfield(opt(1),'limits_min') && ~isempty(opt(1).limits_min) && isfinite(opt(1).limits_min) && ...
            isfield(opt(1),'limits_max') && ~isempty(opt(1).limits_max) && isfinite(opt(1).limits_max) && ...
            isfield(opt(1),'binwidth') && ~isempty(opt(1).binwidth) && isfinite(opt(1).binwidth));
        options_new.edges = options_new.limits_min:options_new.binwidth:options_new.limits_max;

elseif (isfield(opt(1),'edges') && ~isempty(opt(1).edges) && all(isfinite(opt(1).edges)))
    options_new.edges = opt(1).edges;
else
    options_new.edges = [];
end
% if isfield(options_new,'normalization') && ~isempty(options_new.normalization) && strcmpi(options_new.normalization,'density')
%     plot_mode_def.statistics={'mean'};
% end

%% Information
info=               input_data(1,1).info;
% blpxl=      nanmean(info.blpxl(:)); 

plot_mode=...
    idSocial_auxiliaries_setPlotMode(plot_mode,plot_mode_def);
if ~isempty(strfind(lower(plot_mode.statistics{1}),'hist')) && ~strcmpi(plot_mode.normalization,'none')
    plot_mode.ylabel = 'P';
end
plot_mode.extraDims = plot_mode_def.extraDims;

%% Execute function

functionInfo.handle=@idSocial_distanceDistribution;
functionInfo.input_params={'trajectory';...
    options_new.edges;...
    'info.bodylength_in_pixels';...
    options_new.method;...
    options_new.kds_bandwidth;...    
    };

functionInfo.output2function={''; ...
    '';...
    '';...
    'interaction_Distance';...
    'interaction_Distance4Statistics';...
    'interaction_DistanceInfo';...
    };
%% 
plot_mode.edges=options_new.edges;
plot_mode.dist_method=options_new.method;
plot_mode.dist_bandwidth=options_new.kds_bandwidth;
plot_mode.timeintervals_in_min=options_new.timeintervals_in_min;
input_data=idSocial_function_wrapper(input_data,options_new,def_options,plot_mode,functionInfo,strrep(mfilename,'idSocial_',''));
% idSocial_plot(input_data,def_options.act_method,plot_mode);
