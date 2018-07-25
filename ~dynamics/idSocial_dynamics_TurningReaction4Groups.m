function input_data=idSocial_dynamics_TurningReaction4Groups(input_data,options,plot_mode)
% Calculates the turning strength (the component of its acceleration vector 
% perpendicular to the focal's direction of movement) versus the proportion
% of neighbours at the side the focal is turning to.
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
% def_options.timeintervals_in_min=[];
def_options(1).interaction_radius=[-inf inf];
def_options(2).interaction_radius='Interaction radius in BL';


if nargin >= 1
    plot_mode_def.display_mode='hist';
    plot_mode_def.xaxis={'EdgeX'};
    plot_mode_def.xlabel = 'Configuration (N_1:N_2)';
    plot_mode_def.ylabel = 'P(N_2|N_1:N_2)';
    plot_mode_def.statistics = {'Positive_Ratio','Mean','Median','Pool'};
    plot_mode_def.extraDims = {'Focal','','Configuration'};

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
% blpxl=      nanmean(info.blpxl(:)); 
%% Default settings

plot_mode=...
    idSocial_auxiliaries_setPlotMode(plot_mode,plot_mode_def);
plot_mode.extraDims = plot_mode_def.extraDims;

%% Execute function

functionInfo.handle=@idSocial_turningReaction4Groups;
functionInfo.input_params={'trajectory';...
    'info.framerate';...
    'info.bodylength_in_pixels';...
     options_new.interaction_radius;...
    };

functionInfo.output2function={''; ...
    '';...
    '';...
    'dynamics_TurningReaction4Groups';...
    'dynamics_TurningReaction4GroupsInfo'};

%% 

plot_mode.timeintervals_in_min=options_new.timeintervals_in_min;
input_data=idSocial_function_wrapper(input_data,options_new,def_options,plot_mode,functionInfo,strrep(mfilename,'idSocial_',''));
% idSocial_plot(input_data,strrep(mfilename,'idSocial_',''),plot_mode);
