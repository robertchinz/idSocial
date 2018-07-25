function input_data=idSocial_kinematics_Speed(input_data,options,plot_mode)
% Calculates the speed of each individual.
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
%% Default options
def_options = idSocial_auxiliaries_createDefOptions(true);
def_options(1).act_method=strrep(mfilename,'idSocial_','');
def_options(2).act_method='';
def_options(1).limits_min=inf;
def_options(2).limits_min='Min. neighbor distance (in BL)';
def_options(1).limits_max=inf;
def_options(2).limits_max='Max. neighbor distance (in BL)';
def_options(1).binwidth=inf;
def_options(2).binwidth='Binwidth (in BL)';
def_options(1).method={'hist','ksdensity_gauss','ksdensity_epanechnikov','ksdensity_triangular'};
def_options(2).method='Chose between histogram, histogram of absolute values and kernel smoothing function';
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
        plot_mode_def.normalization={'None','No_Frames','No_frames_effective','Density'};
        plot_mode_def.xlabel = 'Speed (BL/s)';
        plot_mode_def.ylabel = 'Counts';
       
    else
        plot_mode_def.statistics = {'Mean','Median','Pool'};
        plot_mode_def.display_mode='plot2d';
        def_options(1).act_method='kinematics_Speed4Statistics';
        plot_mode_def.label = 'Speed (BL/s)';
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

[~, options_new]=idSocial_readparams(input_data,options,def_options,def_options(1).act_method);

if (isfield(opt(1),'limits_min') && ~isempty(opt(1).limits_min) && isfinite(opt(1).limits_min) && ...
            isfield(opt(1),'limits_max') && ~isempty(opt(1).limits_max) && isfinite(opt(1).limits_max) && ...
            isfield(opt(1),'binwidth') && ~isempty(opt(1).binwidth) && isfinite(opt(1).binwidth));
        options_new.edges = options_new.limits_min:options_new.binwidth:options_new.limits_max;

elseif (isfield(opt(1),'edges') && ~isempty(opt(1).edges) && all(isfinite(opt(1).edges)))
    options_new.edges = opt(1).edges;
else
    options_new.edges = [];
end
if isfield(options_new,'normalization') && ~isempty(options_new.normalization) && strcmpi(options_new.normalization,'density')
    plot_mode_def.statistics={'mean'};
end
% %% Fixed settings
% plot_mode.xaxis={'EdgesX'};
% %% Default settings
% plot_mode_def.display_mode='hist';
% plot_mode_def.statistics={'hist'};
plot_mode=...
    idSocial_auxiliaries_setPlotMode(plot_mode,plot_mode_def);


if ~isempty(strfind(lower(plot_mode.statistics{1}),'hist')) && ~strcmpi(plot_mode.normalization,'none')
    plot_mode.ylabel = 'P';
end
%% Execute function

functionInfo.handle=@idSocial_speedDistribution;
functionInfo.input_params={'trajectory';...
    options_new.edges;....
    'info.framerate';...
    'info.bodylength_in_pixels';...
    [];...
    options_new.method;...
    options_new.kds_bandwidth;...
    };
functionInfo.output2function={''; ...
    '';...
    '';...
    'kinematics_Speed';...
    'kinematics_Speed4Statistics';...
    'kinematics_SpeedInfo';...
    };


%% 
plot_mode.edges=options_new.edges; plot_mode.dist_method=options_new.method; plot_mode.dist_bandwidth=options_new.kds_bandwidth;
plot_mode.timeintervals_in_min=options_new.timeintervals_in_min;
input_data=idSocial_function_wrapper(input_data,options_new,def_options,plot_mode,functionInfo,strrep(mfilename,'idSocial_',''));
