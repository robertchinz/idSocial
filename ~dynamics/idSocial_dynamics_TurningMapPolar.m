function input_data=idSocial_dynamics_TurningMapPolar(input_data,options,plot_mode)
% Map of turning strength of a focal individual depending on neighbor
% positions.
%
% In a polar coordinate system in which the focal individual is
% at the origin and moving in direction of the positive y-axis, the
% value of left/right acceleration of the focal in each frame is binned 
% depending on the position of the neighbor. 
% Bins are defined by their length in radial direction 
% (<binwidthRadius>, in Body Length), and the number of slices in the 
% angular direction (<number_of_bins_2pi>). 
% All neighbors closer than <maxRadius> will be taken into account. 
%
% For example, the default values binwidthRadius=1.5, number_of_bins_2pi=8 
% and maxRadius=6 creates bins of size (1.5 BL, pi/4)
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
def_options(1).maxRadius=6;
def_options(2).maxRadius='Max. radius (in BL)';
def_options(1).binwidthRadius=1.5;
def_options(2).binwidthRadius='Radial intervals';
def_options(1).number_of_bins_2pi=8;
def_options(2).number_of_bins_2pi='Angular binwidth is calculated as pi/x';


if nargin >= 1
    plot_mode_def.statistics = {'Positive_Ratio','Mean','Median','Pool'};
    plot_mode_def.display_mode='dynamicMapPolar';
    plot_mode_def.xaxis={'EdgeX'};
    plot_mode_def.yaxis={'EdgeY'};
    plot_mode_def.xlabel = 'Left-Right';
    plot_mode_def.ylabel = 'Behind-Front';
    
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
options_new.edges{1} = 0:options_new.binwidthRadius:options_new.maxRadius+options_new.binwidthRadius;
options_new.edges{2} = 0:2*pi/options_new.number_of_bins_2pi:2*pi+2*pi/options_new.number_of_bins_2pi;

options_new.method = 'hist';
options_new.spacing = 1;
options_new.discard_percentage_outliers=0;
options_new.rotate_z=false;
options_new.normalization='density'; 

plot_mode=...
    idSocial_auxiliaries_setPlotMode(plot_mode,plot_mode_def);


%% Execute function

functionInfo.handle=@idSocial_dynamicMapsPolar;
functionInfo.input_params={'trajectory';...
    options_new.edges;...
    options_new.spacing;...
    'info.framerate';...
    'info.bodylength_in_pixels';...
    options_new.discard_percentage_outliers;...
    options_new.rotate_z;...
    options_new.method;...
    options_new.normalization;...
    };

functionInfo.output2function={'';...
    '';...
    '';...
    'dynamics_PositionMapPolar';...
    'dynamics_TurningMapPolar';...
    'dynamics_AccelerationMapPolar';...
    };

    
%% 
plot_mode.edges=options_new.edges;
plot_mode.dist_method=options_new.method; 
plot_mode.timeintervals_in_min=options_new.timeintervals_in_min;
plot_mode.normalization=options_new.normalization;
input_data=idSocial_function_wrapper(input_data,options,def_options,plot_mode,functionInfo);