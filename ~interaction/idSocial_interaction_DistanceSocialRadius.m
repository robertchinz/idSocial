function input_data=idSocial_interaction_DistanceSocialRadius(input_data,options,plot_mode)
% Calculates the number of neighbors within a given radius around a focal
% individual.
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

def_options(1).social_radius=5; 
def_options(2).social_radius='Limit radius within which neighbors are counted, in BL'; 


if nargin == 1 % Input: options. Output: plot_mode
    opt = input_data;
elseif nargin > 1
    opt = options;
end
if nargin >= 1
    plot_mode_def.display_mode='hist';
    plot_mode_def.xaxis={'EdgeX'};
    plot_mode_def.xlabel = 'Group Size';
    plot_mode_def.ylabel = 'Counts';
    plot_mode_def.statistics = {'hist'};
    plot_mode_def.normalization={'None','No_frames'};
    plot_mode_def.extraDims = {'Focal','','Group Size'};

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


%% Information
info=               input_data(1,1).info;


plot_mode=...
    idSocial_auxiliaries_setPlotMode(plot_mode,plot_mode_def);

if ~strcmpi(plot_mode.normalization,'none')
    plot_mode.ylabel = 'P';
end
plot_mode.extraDims = plot_mode_def.extraDims;
%% Execute function

functionInfo.handle=@idSocial_distanceRatio;
functionInfo.input_params={'trajectory';...
    options_new.social_radius;...
    'info.bodylength_in_pixels';...
    };

functionInfo.output2function={
    'interaction_DistanceSocialRadius';...
    'interaction_DistanceSocialRadius4Statistics';...
    'interaction_DistanceSocialRadiusInfo';...
    };
%% 

input_data=idSocial_function_wrapper(input_data,options_new,def_options,plot_mode,functionInfo);
