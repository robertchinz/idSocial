%% idSocial MANUAL AND EXAMPLES 
% This file contains some help and examples for the usage of
% idSocial.
% You can copy this script and adapt it for your own 
% purposes.
% If you find bugs or have suggestions, please contact
% bugs.idSocial@gmail.com.

%% 0. FIRST STEPS
% To be able to use idSocial, simply copy the folder
% 'idSocial/' to a location of your choice and include the
% folder and all its subfolders in the Matlab(R) path
% by choosing 'File'->'Set Path'->'Add with Subfolders' in
% the Matlab(R) command window or editor.
%
% The basic steps of how to use idSocial to analyse 
% experimental data taken from various trials will be
% explained in Chapter 1 and 2.
% In order to run an easy example or create your own script 
% you can copy the code from sections 1.4 
% (Defining data location), 1.6 (Data preprocessing) and the
% specific function you want to use from Chapter 3 and paste
% it to a new Matlab(R) script, or see  
% <matlab:edit(fullfile(fileparts(which('idSocial_manual.m')),'idSocial_example.m')) idSocial_example.m> 

%% 1. DATA PREPROCESSING: Set data location and basic options
% 
%% 1.1 Input trajectory location
% 
% You have to tell idSocial the location of your data, i.e.,
% the trajectory files. If you are using idTracker, the 
% files are called 'trajectories.mat' and can be found in 
% the subdirectory '..\segm' within the directory of your
% videos.
% Each trajectory file corresponds to a single trial and has
% the format of idTracker output (a numerical Matlab(R) 
% array of size (# of frames) -by- (# of individuals) -by- 2).
% For a single trial, simply define a variable containing
% the full path to the trajectory file:
%
%   trajectory_location = ...
%           'path\trajectory_1.mat';
%
% To create a list of various trials use
% 
%   trajectory_location = ...
%           {'path\trajectory_1.mat' ...    % Trial 1
%            'path\trajectory_2.mat' ...    % Trial 2
%            ...
%            'path\trajectory_N.mat' ...    % Trial N
%           };
% 
% Note that the actual name of the trajectory file is not
% important, and different files (in different locations) 
% can have the same name. 
%% 1.2 Defining groups and subgroups of data
% 
% idSocial lets you compare data taken for 
% 
% # various groups (e.g., trials for test group and control 
%   group) and lets you divide your data in 
% # different subgroups (e.g., trials from different days, 
%   weeks).
% 
% Example:
% During two days, we record three videos each day for one 
% test group and a second control group.   
% In order to tell idSocial which trials belong to which 
% group or subgroup, we create a list of the following form:
% 
% Group 1 (Test), Day 1:  
%   (In the following, Group_1 = 1 stands for 'Test group',
%   Day_1 = 1 for 'Day 1')
%
%   trajectory_location{Group_1}{Day_1} = ...
%          {'trajectory_1.mat' ...  % Gr. 1, Day 1, Tr. 1
%           'trajectory_2.mat' ...  % Gr. 1, Day 1, Tr. 2
%           'trajectory_3.mat' ...  % Gr. 1, Day 1, Tr. 3
%           };
% 
% Group 1 (Test), Day 2: 
%
%   trajectory_location{Group_1}{Day_2} = ...
%          {'trajectory_1.mat' ...  % Gr. 1, Day 2, Tr. 1
%           'trajectory_2.mat' ...  % Gr. 1, Day 2, Tr. 2
%           'trajectory_3.mat' ...  % Gr. 1, Day 2, Tr. 3
%           };
%
% Group 2 (Control), Day 1:  
%
%   trajectory_location{Group_2}{Day_1} = ...
%          {'trajectory_1.mat' ...  % Gr. 2, Day 1, Tr. 1
%           'trajectory_2.mat' ...  % Gr. 2, Day 1, Tr. 2
%           'trajectory_3.mat' ...  % Gr. 2, Day 1, Tr. 3
%           };
%
% Group 2 (Control), Day 2: 
%
%   trajectory_location{Group_2}{Day_2} = ...
%          {'trajectory_1.mat' ...  % Gr. 2, Day 2, Tr. 1
%           'trajectory_2.mat' ...  % Gr. 2, Day 2, Tr. 2
%           'trajectory_3.mat' ...  % Gr. 2, Day 2, Tr. 3
%           };
%
% Note that subgroups of different groups do not have to
% contain the same number of trials, and different groups do
% not have to contain the same number of subgroups.
% In our example, the number of experiments per day could
% vary, and the number of days did not have to be the same
% for test and control group.
%% 1.3 Optional information: |datosegm| (for idTracker users):
% 
% Analogously, if you work with idTracker, you can create
% lists for the locations of the 'datosegm.mat'-files, which
% contain additional data like the framerate, body length 
% and setup dimensions. For each 'trajectories.mat'-file 
% there is one corresponding 'datosegm.mat'.
% If you do not tell idSocial the location, it will 
% automatically look for it in the folder of the 
% corresponding trajectory-file.
% idSocial will use information about the frame rate 
% and the body length of the individuals for each video if 
% the corresponding data is present in 'datosegm.mat'.
% Otherwise, you have to set frame rate and body length
% manually (see 1.6 Example: Data preprocessing,
% options.framerate and options.blpxl).
%% 1.4 Example: Defining data location

trajectory_location=[];datosegm_location=[];

trajectory_location{1}{1}={'trayectoria_trial_orden_1.mat'...
    'trayectoria_trial_orden_2.mat',...
    'trayectoria_trial_orden_3.mat',...
   };
datosegm_location{1}{1}={'datosegm_trial_orden_1.mat'...
    'datosegm_trial_orden_2.mat',...
    'datosegm_trial_orden_3.mat',...
    };
%% 1.5 Setting options for data preprocessing
% idSocial preprocesses the data before you can start the 
% actual analysis. 
% The following parameters determine which preprocessing is 
% applied. They are stored in the variable 'options' which
% is of the form
% 
%       options.start_min = ...;
%       options.project_path = ...;
%       ...
%
% It is also possible to set seperate options for each trial 
% In this case, options must have a structure similar to 
% 'trajectory_location' in Section 1.2:
% 
% Group 1, Day 1, Trial 1
%       options{Group_1}{Day_1}{Trial_1}.start_min = ...;
%       options{Group_1}{Day_1}{Trial_1}.blpxl = ...; 
%       ...
% Group 1, Day 1, Trial 2
%       options{Group_1}{Day_1}{Trial_2}.start_min = ...;
%       options{Group_1}{Day_1}{Trial_2}.blpxl = ...; 
%       ...
%
% Furthermore, if you want to set the same options for all
% trials belonging to same Group/Subgroup, you can use the
% following syntax:
%
% Same options for all trials of the same group:
% 
%       options{Group_1} = options1;    
%       options{Group_2} = options2; 
%       ...
% 
% Same options for all trials of the same subgroup (day):
% 
%       options{Group_1}{Day_1} = options1;    
%       options{Group_1}{Day_2} = options2; 
%       ...
%
% where options1 and options2 are Matlab(R) structures with 
% possible fieldnames and values shown in 1.6. In general,
% not all options have to be set by the user. Options which
% are not set manually will take on default values. 

%% 1.6 Example: Data preprocessing  
% Not all options have to be set by the user. Options which
% are not set manually will take on default values. An
% exception is 'project_path'. In case it is undefined when 
% starting data preprocessing, a dialog box will ask the 
% user to chose a storage location for temporary files.   

clear options           % This deletes any existing variable
                        % 'options' 
options.start_min = 0;  % Only frames from a later time (in
                        % minutes) will be analysed.
                        % Default: 0.
options.end_min = 180;  % Frames from later times (in
                        % minutes) will be discarded.
                        % Default: inf.
options.project_path = ''; % If this is empty, no 
                        % figures will be saved. If it is a
                        % valid path, figures will be saved 
                        % to path\figures and 
                        % path\figures_latex.
                        % Default: ''.
options.blpxl = 50;     % Body length of the individuals in 
                        % pixels.
                        % If set to 'idTracker' and 
                        % idTracker has been used to 
                        % generate the trajectory-files
                        % and idSocial finds the 
                        % corresponding 'datosegm'-files, 
                        % body length will be taken directly
                        % from 'datosegm.mat' for each 
                        % trial.
                        % blpxl is used to present results 
                        % in units of body length.
                        % Default: 1.
options.framerate = 30; % Framerate of the videos from which
                        % trajectories have been obtained.
                        % If idTracker is used to generate 
                        % the trajectory-files, and idSocial
                        % finds the corresponding 
                        % 'datosegm'-files, framerate 
                        % will be taken directly from 
                        % 'datosegm.mat' for each trial. 
                        % framerate is used to present 
                        % results in units of 1 min/sec. 
                        % Default: 1.
options.interpolate_trajectories = false; % If set 'true', 
                        % missing data points will be 
                        % interpolated. NOTE: Only use 
                        % interpolation if you are sure it 
                        % will not create 'fictional' data.
                        % Default: false.
options.interpolation_mode = 'linear'; % Set interpolation
                        % method. Possible options:
                        % 'linear' (default)
                        % 'spline' 
                        % Default: 'linear'.
options.smooth_method = 'moving';% Set smoothing method. Is 
                        % left empty no smoothing will be 
                        % performed. Possible options:
                        % 'moving' (default)
                        % 'lowess'
                        % 'loess'
                        % 'rlowess'
                        % 'rloess'
                        % 'spline'
                        % Default: 'moving'.
                        % See also: MATLAB(R) 'smooth'
options.smooth_degree = 5; % Size of the moving window used for
                        % smoothing, i.e., number of 
                        % neighboring data points used for
                        % calculation of adjusted central 
                        % data point.
                        % Default: 1.
options.filter_minProbIdentityAssignment = 1.0;
                        % Lower limit for the reliability of
                        % identities assigned by idTracker.
                        % This option does only apply if
                        % idTracker has been used to extract
                        % trajectories.
                        % For each frame and each
                        % individual, idTracker returns a
                        % value between 0 and 1 which 
                        % reflects the reliability of the 
                        % assignment (0: Not reliable, 
                        % 1: Identities are completely 
                        % reliable).
                        % In addition, values of -1 or -2 
                        % reflect the way idTracker has 
                        % dealt with crossings of two or 
                        % more individuals, where identities
                        % cannot be determined: 
                        % Value of -1 means there is a 
                        % crossing, but blobs can be 
                        % recovered by image erosion, and 
                        % the center of mass coordinates of 
                        % those blobs are interpolated in 
                        % order to get approximate 
                        % coordinates of the overlapping 
                        % individuals during the crossing.
                        % A value of -2 has a similar
                        % meaning to -1, but image erosion 
                        % does not recover seperate blobs. 
                        % Interoplation is done under the 
                        % only condition that the 
                        % interpolated trajectory 
                        % coordinates lie within the blob 
                        % in a certain frame.
                        % Dafault: 1.
options.filter_focal_speedlimits_bl_per_s = [-inf inf];
                        % Lower and upper limit of the
                        % focal individual's speed. Data 
                        % corresponding to the focal will be
                        % neglected in frames in which the 
                        % focal moves slower than the lower
                        % or faster than the upper limit.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [-inf inf].
options.filter_neighbor_speedlimits_bl_per_s = [-inf inf];
                        % Lower and upper limit of the
                        % neighbor individual's speed. Data 
                        % corresponding to the focal will be
                        % neglected in frames in which the 
                        % neighbor moves slower than the 
                        % lower or faster than the upper 
                        % limit.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [-inf inf].
options.filter_focal_accelerationlimits_bl_per_s2 = [-inf inf];
                        % Lower and upper limit of the
                        % focal individual's acceleration. 
                        % Data corresponding to the focal 
                        % will be neglected in frames in 
                        % which the focal accelerates slower
                        % than the lower or faster than the 
                        % upper limit.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [-inf inf].
options.filter_neighbor_accelerationlimits_bl_per_s2 = [-inf inf];
                        % Lower and upper limit of the
                        % neighbor individual's acceleration. 
                        % Data corresponding to the focal 
                        % will be neglected in frames in 
                        % which the focal accelerates slower
                        % than the lower or faster than the 
                        % upper limit.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [-inf inf].
options.filter_distancelimits_bl = [-inf inf];
                        % Lower and upper limit of the
                        % distance between the focal 
                        % individual and its neighbor. Data 
                        % corresponding to the focal will be
                        % neglected in frames in which the
                        % distance between the both is
                        % smaller than the lower or greater
                        % than the upper limit.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [-inf inf].
options.filter_focal_rectangularROI = [-inf inf; ...
                                            -inf inf;...
                                            -inf inf];
                        % Lower and upper limit for x-, y- 
                        % and z-coordinates (in pixels). 
                        % Focal coordinates exceeding these 
                        % limits will be discarded. This 
                        % lets you set a rectangular region 
                        % of interest for focal coordinates.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [-inf inf; -inf inf;...
                        %                       -inf inf];
options.filter_neighbor_rectangularROI = [-inf inf; ...
                                            -inf inf;...
                                            -inf inf];
                        % Lower and upper limit for x-, y- 
                        % and z-coordinates (in pixels).
                        % Neighbor coordinates exceeding 
                        % these limits will be discarded. 
                        % This lets you set a rectangular 
                        % region of interest for neighbor 
                        % coordinates.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [-inf inf; -inf inf;...
                        %                       -inf inf];
options.filter_focal_circularROI = [0,0,0,inf];
                        % x-, y-, z-coordinates of the 
                        % center and the radius (in pixels) 
                        % of a circular region of interest. 
                        % Focal coordinates outside this 
                        % region will be discarded. 
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [0,0,0,inf].
options.filter_neighbor_circularROI = [0,0,0,inf];
                        % x-, y-, z-coordinates of the 
                        % center and the radius (in pixels) 
                        % of a circular region of interest. 
                        % Neighbor coordinates outside this 
                        % region will be discarded.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: [0,0,0,inf].
options.order_neighbors = false;
                        % If set 'true', neighbors will be 
                        % sorted according to their distance 
                        % to the focal. The neighbor with 
                        % index 1 will be the nearest 
                        % neighbor, neighbor 2 the second 
                        % nearest etc., up to the neighbor 
                        % with index '# of individuals - 1' 
                        % (the neighbor with index 
                        % '# of individuals' corresponds to 
                        % the focal itself).
                        % NOTE: Sorting neighbors at the
                        % preprocessing stage will fix the 
                        % order for all later analysis, even 
                        % if the options.order_neighbors for 
                        % a particular function is not
                        % activated. Setting this option at 
                        % the preprocessing stage lets you 
                        % apply filters to selected 
                        % groups/subgroups/etc.
                        % Default: false.
options.significance_between_groups=false;
                        % If set 'true', color markers
                        % will indicate if plotted data sets
                        % are significantly different.
                        % NOTE: The same filter options
                        % can also be applied at a later
                        % point when calling a particular
                        % function. Filtering at the
                        % preprocessing stage will fix the
                        % used filters for all later usage
                        % but lets you apply filters to
                        % selected groups/subgroups/etc.
                        % Default: false.
options.group_names = {'Group 1'; 'Group 2';}; 
                        % If trials  are sorted into 
                        % different groups, names for each 
                        % group can be set (e.g., 
                        % {'Test'; 'Control'};
                        % Default: ''.
options.temp_savepath = 'D:\~idSocialTemp\'; % idSocial will
                        % save results temporarily at the
                        % given location. If left empty, 
                        % you will be asked to choose the 
                        % location when executing any 
                        % idSocial method for the first time 
                        % with the given data.
                        % Default: If not set, user will be
                        % asked to enter a location.
options.random_data = true; % Create artifical trajectories 
                        % by randomizing the temporal order
                        % of the original trajectory data. 
                        % These artificial trajectories
                        % can later be used to compare the
                        % original data to the randomized
                        % controls. 
                        % See also options.random_data 
                        % below.
                        % Default: false.

[data, ~, info]=idSocial_loadData(trajectory_location, datosegm_location, options);
%% ---------------------------------------------------------
%% 2. DATA ANALYSIS
% After preprocessing the data, you are ready to analyse it. 
% idSocial consists of a collection of methods for the 
% analysis of a number of different aspects of animal
% behavior. 
% The available methods can be divided into different fields 
% of analysis: Kinematics (basic movement measures like 
% speed, acceleration etc.), Interaction (basic measures of 
% sociability like inter-individual distance), dynamics
% (methods which analyse individual behavior connected with 
% social interaction), hierarchies (analysis of aggression, 
% territoriality and leader-follower relationship), 
% collective (advanced analysis of group behavior) and 
% network (focusing on network properties of social 
% interactions).
%
% Each idSocial-function takes three input parameters: 
%
% # |data|, which is the structure returned by the data 
% preprocessing step executed by |idSocial_loadData|, and
% also the output returned by every idSocial-function. It 
% stores all the information and calculations of your 
% recent idSocial calculations and therefore should be 
% updated by every function call. This update step is done
% by storing the output in |data|.  
% # |options|, which contains parameters for the specific 
% function.
% # |plot_mode|, which lets you choose how you want to slice 
% or pool your data.
% 
% Example:
%
%  data = ...
%   idSocial_interaction_Distance(data, options, plot_mode)
%% 2.1 Input parameter |data|
% The output structure returned by |idSocial_loadData|, 
% containing, above all, the trajectory data and, possibly
% the results from previous calls of any idSocial method. 
%% 2.2 Input parameter |options|
% A structure containing parameters necessary for
% the respective method, as well as additional parameters
% regarding the presentation of the results.
% Parameters common to all methods are:

options.timeintervals_in_min = 2; 
                % idSocial divides each 
                % trajectory into shorter pieces
                % corresponding to periods of time set by  
                % timeintervals_in_min.
                % The average for each of these periods will
                % be calculated (and, in case of a basic 
                % line plot, presented with the respective 
                % standard error).
                % Note that idSocial uses overlapping
                % periods: For timeintervals_in_min = 2,
                % it will bin data from 0-2 min, 1-3 min,
                % 2-4 min, etc.
                % If options.timeintervals_in_min = [],
                % trajectories are divided in three periods,
                % where the 1st one includes times from 0
                % to half the total length of the
                % trajectory, the 2nd one includes the whole
                % trajectory, and the 3rd all the times from
                % the middle of the trajectory to the end.
options.random_data = true; 
                % Include results obtained by
                % randomizing the temporal order of the 
                % original trajectory data to the plot
                % in order to create a control of the 
                % meaningfulness of the results.
                % To be able to use this option, 
                % idSocial_loadData must also have been 
                % executed with the variable 
                % options.random_data = true.
options.plot_axislabelstring =  {'x-label'; 'y-label'};
                % Lets you choose the label of the x- and 
                % y-axis.
                % If empty, most of the functions
                % will choose appropriate labels.
options.plot_titlestring = {'Title'};
                % Lets you choose the title of the final
                % plot.
options.filter_focal_speedlimits_bl_per_s = [-inf inf];
                % Lower and upper limit of the focal 
                % individual's speed. Data corresponding to 
                % the focal will be neglected in frames in 
                % which the focal moves slower than the 
                % lower or faster than the upper limit. 
                % NOTE: The same filter options can also be
                % applied at the preprocessing stage using 
                % |idSocial_loadData|. Filtering at the 
                % preprocessing stage will fix the used 
                % filters for all later usage but lets you 
                % apply filters to selected 
                % groups/subgroups/etc., while setting 
                % filters in the options for a particular 
                % function ill affect all the data (all
                % groups/subgroups/trials), but only for the
                % particular function.
                % Default: [-inf inf].
options.filter_neighbor_speedlimits_bl_per_s = [-inf inf];
                % Lower and upper limit of the neighbor 
                % individual's speed. Data corresponding to 
                % the focal will be neglected in frames in 
                % which the neighbor moves slower than the 
                % lower or faster than the upper limit.
                % NOTE: The same filter options can also be 
                % applied at the preprocessing stage using 
                % |idSocial_loadData|. Filtering at the 
                % preprocessing stage will fix the used 
                % filters for all later usage but lets you 
                % apply filters to selected 
                % groups/subgroups/etc., while setting 
                % filters in the options for a particular 
                % function will affect all the data (all
                % groups/subgroups/trials), but only for the
                % particular function.
                % Default: [-inf inf].
options.filter_distancelimits_bl = [-inf inf];
                % Lower and upper limit of the distance 
                % between the focal individual and its 
                % neighbor. Data  corresponding to the focal
                % will be neglected in frames in which the 
                % distance between the both is smaller than 
                % the lower or greater than the upper limit.
                % NOTE: The same filter options can also be 
                % applied at the preprocessing stage using 
                % |idSocial_loadData|. Filtering at the 
                % preprocessing stage will fix the used 
                % filters for all later usage but lets you 
                % apply filters to selected 
                % groups/subgroups/etc., while setting 
                % filters in the options for a particular 
                % function will affect all the data (all 
                % groups/subgroups/trials), but only for the
                % particular function.
                % Default: [-inf inf].
options.filter_focal_rectangularROI = [-inf inf; ...
                                            -inf inf;...
                                            -inf inf];
                % Lower and upper limit for x-, y- and z-
                % coordinates (in pixels). Focal coordinates
                % exceeding these limits will be discarded. 
                % This lets you set a rectangular region of 
                % interest for focal coordinates.
                % NOTE: The same filter options can also be 
                % applied at the preprocessing stage using 
                % |idSocial_loadData|. Filtering at the 
                % preprocessing stage will fix the used 
                % filters for all later usage but lets you 
                % apply filters to selected 
                % groups/subgroups/etc., while setting 
                % filters in the options for a particular 
                % function will affect all the data (all 
                % groups/subgroups/trials), but only for the
                % particular function.
                % Default: [-inf inf; -inf inf;...
                %                       -inf inf];
options.filter_neighbor_rectangularROI = [-inf inf; ...
                                            -inf inf;...
                                            -inf inf];
                % Lower and upper limit for x-, y- and z-
                % coordinates (in pixels). Neighbor 
                % coordinates exceeding these limits will be
                % discarded. This lets you set a rectangular
                % region of interest for neighbor 
                % coordinates.
                % NOTE: The same filter options can also be
                % applied at the preprocessing stage using
                % |idSocial_loadData|. Filtering at the 
                % preprocessing stage will fix the used 
                % filters for all later usage but lets you 
                % apply filters to selected 
                % groups/subgroups/etc., while setting 
                % filters in the options for a particular 
                % function will affect all the data (all 
                % groups/subgroups/trials), but only for the
                % particular function.
                % Default: [-inf inf; -inf inf;...
                %                       -inf inf];
options.filter_focal_circularROI = [0,0,0,inf];
                % x-, y-, z-coordinates of the center and 
                % the radius (in pixels) of a circular 
                % region of interest.  Focal coordinates 
                % outside this region will be discarded. 
                % NOTE: The same filter options can also be 
                % applied at the preprocessing stage using 
                % |idSocial_loadData|. Filtering at the 
                % preprocessing stage will fix the used 
                % filters for all later usage but lets you 
                % apply filters to 
                % selected groups/subgroups/etc., while 
                % setting filters in the options for a 
                % particular function will affect all the 
                % data (all groups/subgroups/trials), but 
                % only for the particular function.
                % Default: [0,0,0,inf].
options.filter_neighbor_circularROI = [0,0,0,inf];
                % x-, y-, z-coordinates of the center and 
                % the radius (in pixels) of a circular 
                % region of interest. Neighbor coordinates 
                % outside this region will be discarded.
                % NOTE: The same filter options can also be 
                % applied at the preprocessing stage using 
                % |idSocial_loadData|. Filtering at the 
                % preprocessing stage will fix the used 
                % filters for all later usage but lets you 
                % apply filters to selected 
                % groups/subgroups/etc., while setting 
                % filters in the options for a particular 
                % function will affect all the data (all 
                % groups/subgroups/trials), but only for the
                % particular function.
                % Default: [0,0,0,inf].
options.order_neighbors = false;
                % If set 'true', neighbors will be sorted 
                % according to their distance to the focal.
                % The neighbor with index 1 will be the 
                % nearest neighbor, neighbor 2 the second 
                % nearest etc., up to the neighbor with 
                % index '# of individuals - 1' (the neighbor
                % with index '# of individuals' corresponds 
                % to  the focal itself).
                % NOTE: Sorting neighbors at the 
                % preprocessing stage will fix the  order 
                % for all later analysis, even if the 
                % options.order_neighbors for a particular 
                % function is not activated. Setting this 
                % option at the preprocessing stage lets you 
                % apply filters to selected 
                % groups/subgroups/etc., while setting it 
                % for a particular function will only affect
                % this particular function, but all the data
                % (all groups/subgroups/trials).
                % Default: false.
options.statistical_test_type = 'ranksum';
                % Statistical test used to determine if 
                % displayed data sets are statistically    
                % different. Options are: 
                % ranksum: two-sided Wilcoxon rank sum
                %          test; see Matlab(R) ranksum.
                % ttest: two-sample t-test
                %         See Matlab(R) ttest2.
                % bootstrap: Bootstrapping.
                %           See idSocial_bootstrap.
options.statistical_test_significance = 0.05; 
                % Significance level. If the p-Value is less
                % than the chosen significance level,
                % according pairs of data points will be 
                % marked as statistically different.
options.plot_labelfontsize = 14;
                % Title font size. 
options.plot_axislabelsize = 16;
                % Font size of axis labels.
options.plot_tickfontsize = 14;
                % Font size of x- and y-axis tick labels  
options.plot_bootstrap_statfunc = 'Mean';
                % Only applicable if 
                % statistical_test_type = 'bootstrap'.
                % Determines on which statistical measure
                % will be used for the bootstrapping.
                % Options are 'Mean' and 'Median'.
options.plot_bootstrap_repetitions = 10000;
                % Only applicable if 
                % statistical_test_type = 'bootstrap'.
                % Determines the number of repetitions,
                % i.e., of bootrap samples.
options.plot_deviation = true; 
                % If 'true', error bars will be plotted
                % showing the standard error (if the 
                % displayed values are means) or
                % confidence interval (for medians).
% @TODO: There are options missing here!@
%% 2.3 Input parameter |plot_mode|
% Defines how the results will be displayed.
%
% In general, we want to plot and compare different pieces 
% of all the data. This is done by setting |plot_mode|,
% which uses string specifiers to distinguish subsets 
% of data:
% 'Group', 'Day', 'Trial', 'Time', 'Focal', 'Neighbor',
%                                           'EdgeX', 'EdgeY'
%
% |plot_mode| is a structure with the following fields: 
% 
%  plot_mode.xaxis: Specifies the x-axis: 
%               plot_mode.xaxis = {'Time'}
%               will plot the data versus time.
%               Some functions do not permit changing the 
%               axes, for example in case of distributions/
%               histograms, the x-axis is set to 'EdgeX' 
%               automatically, corresponding to the bin 
%               edges set before running the function (and,
%               analogously, the y-axis is set to 'EdgeY'
%               in case of bivariate histograms/maps).
%
%  plot_mode.yaxis: see XAXIS. Only applicable to bivariate 
%               histograms/maps.
%
%  plot_mode.filter: Specifies/filters data for calculations. 
%               Filter is a cell array with as many 
%               rows as the number of data subsets you want 
%               to compare (e.g., the number of legend 
%               entries in a plot). Each row has the form
%
%                    {'Specifier 1' [indices] ...
%                          'Specifier 2' [indices] ...}
%               where 'Specifier' is a string specifier 
%               from the list above, and 'indices' determine
%               the data indices which will be used for the
%               calculations. E.g., 
%
%                    plot_mode.filter =
%                           {'Group' 1 'Trial' 1; ...
%                           'Group' 2 'Trial' 3 }
%
%               corresponds to data from Group 1, 
%               Trial 1 in the first plot, and Group 2, 
%               Trial 3 in the second.
%               
%               If the field filter does not exist or is 
%               empty, all the data will be used.
%
%  plot_mode.data: Specifies which statistical operation will 
%               be applied to which dimension of OUTPUT.
%               Like FILTER, it is a cell array with as many 
%               rows as the number of data subsets you want 
%               to compare.
%               Each row has the form
%                   {'Stat. op.' 'Specifier 1' ...
%                          'Stat. op.' 'Specifier 2' ...}
%               Statistical operation may be 'Mean',
%               'Median', 'Min', 'Max', 'Hist' or 'Pool'.
%               E.g., 
%
%                    plot_mode.data =
%                           {'Mean' 'Trial'; ...
%                           'Median' 'Trial'}
%
%               will calculate the average over all trials 
%               (which are left after filtering) for the 
%               first plot, and the median of all trials for
%               the second.
%               
%               If the field data does not exist or is 
%               empty, all the data will be pooled together
%               and the mean is calculated.
%
%  plot_mode.autocombs: Instead of selecting each plot entry 
%               seperately by defining various rows in 
%               FILTER, AUTOCOMBS can be used to 
%               automatically generate all the possible 
%               index combinations for a list of string 
%               specifiers. 
%               E.g., if #Groups = 2 and #Trials = 2 
%
%                    plot_mode.autocombs =
%                           {'Group Trial'}
%
%               will generate combinations 
%               (Group 1, Trial 1), (Group 1, Trial 2),  
%               (Group 2, Trial 1), (Group 2, Trial 2) 
%               and plot the corresponding data separately 
%               (in the same figure window).
%
%               Note that if FILTER is defined, combinations 
%               are only generated from the indices present
%               in FILTER. In this case, FILTER and DATA can
%               only have one row, and the statistical 
%               operation(s) in data will be applied to all 
%               combinations.   
%
%  plot_mode.legendstring: Can be used to manually set the 
%               legend string. For each row in FILTER/each 
%               combination generated by autocombs, there
%               must be one row in LEGENDSTRING:
%                   plot_mode.legendstring =
%                           {'Mean over trials'; ...
%                           'Median of trials'}               
%
%  plot_mode.statistics does not apply to all functions, 
%               since for many functions it is 
%               fixed. It determines the statistical 
%               operation which is applied directly to the 
%               output of the core function. This only makes
%               sense if the output is a numeric vector or a 
%               n -by- m cell of numeric vectors (e.g., a 1 
%               -by- m cell where each element of the cell  
%               corresponds to the bin of a histogram and 
%               contains all the values falling into this 
%               bin). 
%               The statistical operation can be 'Mean', 
%               'Median' and 'Hist' (the latter simply 
%               counts the number of values).
%               E.g., 
%
%                    plot_mode.statistics =
%                           'Mean';           
%% ---------------------------------------------------------
%% 3. FUNCTION OVERVIEW
%% 3.1 KINEMATICS - basic movement measures
%% Speed
% |idSocial_kinematics_Speed| calculates the average 
% absolute velocity over time (i.e., the mean over all
% frames in each time interval defined by 
% |options.timeintervals_in_min|). Use |plot_mode.xaxis| to
% choose units on the X-axis (Group, Day, Trial, Time, 
% Focal, Neighbor).
 
% Clear previous options
options = []; plot_mode = [];
% Set new options
options.random_data = false;
options.timeintervals_in_min = 2; 
% Choose what to plot 
plot_mode.xaxis={'Time'};
plot_mode.legend= {'Group'};
% plot_mode.filter={'Group' 1 'Trial' 1 'Focal' 1};
% Execute function
data=idSocial_kinematics_Speed(data,options,plot_mode);
%% Speed distribution
% |idSocial_kinematics_SpeedDistribution| plots the
% velocity distribution. It also shows the median(s), 
% and in case of multiple data sets, tests if the medians 
% are significantly different. _Additional options_:
% |options.edges| sets the edges of the bins for the
% distribution. 

% Clear previous options
options = [];
% Set new options
options.random_data = false;
options.timeintervals_in_min = 2;
% Extra options for idSocial_kinematics_SpeedDistribution 
options.edges = 0:.5:40; 
% Choose what to plot 
plot_mode.legend= {'Group'};
% Execute function
data=idSocial_kinematics_SpeedDistribution(data,options,plot_mode);
%% Acceleration
% |idSocial_kinematics_Acceleration| calculates the average 
% absolute acceleration over time (i.e., the mean over all
% frames in each time interval). Use |plot_mode.xaxis| to
% choose units on the X-axis (Group, Day, Trial, Time, 
% Focal, Neighbor).
 
% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = false;
options.timeintervals_in_min = 2; 
% Choose what to plot 
plot_mode.xaxis={'Time'};
plot_mode.legend= {'Group'};
% Execute function
data=idSocial_kinematics_Acceleration(data,options,plot_mode);
%% Acceleration distribution
% |idSocial_kinematics_AccelerationDistribution| plots the
% acceleration distribution. It also shows the median(s), 
% and in case of multiple data sets, tests if the medians 
% are significantly different. _Additional options_:
% |options.edges| sets the edges of the bins for the
% distribution. 

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = false;
options.timeintervals_in_min = 2;
% Extra options for idSocial_kinematics_SpeedDistribution 
options.edges = 0:.5:40; 
% Choose what to plot 
plot_mode.legend= {'Group'};
% Execute function
data=idSocial_kinematics_AccelerationDistribution(data,options,plot_mode);
%% Average time to reach a given point in the focal's coordinate system
% |idSocial_kinematics_ApproachTimeMap| calculates the
% average time an individual needs to reach a certain point
% in the coordinate system defined by the individuals
% direction of movement (y-axis) and its position (origin).
% This is done by mapping the x-y-position of the focal at
% time T + dt onto the the focal's coordinate system at time
% T. 
% _Additional options_: |options.edges| sets the edges of the
% two dimensional bins into which the x-y-coordinates are 
% sorted. |options.spacing| sets the spacing for 
% overlapping bins which can be applied to smooth the 
% result. For |spacing|=2 (N), contiguous bins previously 
% defined by |edges| will be seperated by one (N-1) 
% additional bins whos edges are situated at the centers of
% the original bins. |options.interval_in_seconds| 
% sets the maximum time dt in seconds which is taken 
% into account. Note that the longer the time, the less 
% meaningful will be the result for close distances of the
% focal at time T because the individual will have time to 
% visit the same positions various times. 
% |options.frame_step| determines the temporal distance
% between frames for which the focal's cordinate system is 
% calculated and the mapping is done. Higher values reduce 
% the processing time at the expense of higher precision.   

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = false;
options.timeintervals_in_min = [];
% Extra options for idSocial_kinematics_ApproachTimeMap
options.edges = {-10:1:10; -10:1:10};
options.spacing = 3;
options.interval_in_seconds = 10;
options.frame_step = 30;
% Choose what to plot 
plot_mode.legend= {'Group'};
plot_mode.filter={'Group' 1;};
% Execute function
data=idSocial_kinematics_ApproachTimeMap(data,options,plot_mode);
%% 3.2 INTERACTION - basic measures of sociability 
%
%% Inter-individual distance (Mean)
% |idSocial_interaction_Distance| calculates the mean 
% distance over time (i.e., the mean over all
% frames in each time interval defined by 
% |options.timeintervals_in_min|) between all possible pairs
% of individuals.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=5;
options.random_data = false;
% Choose what to plot
plot_mode.xaxis={'Time'};
plot_mode.legendstring={'Group 1';};
plot_mode.filter={'Group' 1;};
% Execute function
data=idSocial_interaction_Distance(data,options,plot_mode);
%% Inter-individual distance (Median)
% |idSocial_interaction_DistanceMedian| calculates the 
% median distance over time (i.e., the median of distance 
% distributions of all frames in each time interval defined 
% by |options.timeintervals_in_min|), where distances are
% calculated for all possible pairs of individuals.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=5;
options.random_data = true;
% Choose what to plot
plot_mode.xaxis={'Time'};
plot_mode.legendstring={'Group'};
% Execute function
data=idSocial_interaction_DistanceMedian(data,options,plot_mode);
%% Distance distributions
% |idSocial_kinematics_DistanceDistribution| plots the
% distribution of distances between all possible pairs of 
% individuals. It also shows the median(s), 
% and in case of multiple data sets, tests if the medians 
% are significantly different. _Additional options_:
% |options.edges| sets the edges of the bins for the
% distribution. 

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = true;
options.timeintervals_in_min = [];
% Extra options for idSocial_interaction_DistanceDistribution
options.edges = 0:.5:40; 
% Choose what to plot 
plot_mode.legend= {'Group'};
% Execute function
data=idSocial_interaction_DistanceDistribution(data,options,plot_mode);
%% Group polarization
% |idSocial_interaction_Polarization| calculates the 
% mean group polarization over time (i.e., the mean
% polarization over all frames in each time interval defined 
% by |options.timeintervals_in_min|).

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=5;
options.random_data = true;
% Choose what to plot
plot_mode.xaxis={'Time'};
plot_mode.legendstring={'Group'};
% Execute function
data=idSocial_interaction_Polarization(data,options,plot_mode);
%% Map of probabilities for relative orientation
% |idSocial_dynamics_AlignmentMap| depicts the probability of
% relative orientations of focal and neighbor individual at
% a given inter-individual distance.
% Distances and the angle between focal movement direction 
% and neighbor movement direction are binned into spatial intervals with
% edges given by |options.edges|. 
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. ||options.spacing| sets the span of
% bin overlap to smooth the resulting plot.
%
% <<focal_system_angles.png>>

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_AlignmentMap
options.edges={0:.5:5; 0:.5:2*pi}; 
options.spacing=1;
options.normalize_distances=false;
% Execute function
data=idSocial_interaction_AlignmentMap(data,options,plot_mode);
%% 3.3 DYNAMICS - Social influence on individual behavior 
%% Turning strength depending on relative lateral neighbor distance
% |idSocial_dynamics_TurningReaction| calculates the average 
% strength of left/right turns of a focal individual
% depending on the lateral distance of the neighbor. This is 
% visualized by mapping the focal's turning acceleration 
% onto a discrete spatial interval corresponding to the
% lateral position of the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span for a
% moving average filter to smooth the resulting plot.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = true;
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_TurningReaction
options.edges=-11.25:.5:11.25; 
options.spacing=3;
% Execute function
data=idSocial_dynamics_TurningReaction(data,options,plot_mode);
%% Turning strength depending on neighbor presence
% |idSocial_dynamics_TurningBalance| calculates and plots
% the average strength of turns of a focal individual
% towards the side where there is a neighbor present, and 
% the averge strength of turns of a focal individual
% towards the side where there is no neighbor present. This 
% is visualized by mapping the focal's turning acceleration 
% onto a discrete spatial interval corresponding to the
% lateral position of the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span for a
% moving average filter to smooth the resulting plot.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = true;
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_TurningReaction
options.edges=-11.25:.5:11.25; 
options.spacing=2;
% Execute function
data=idSocial_dynamics_TurningBalance(data,options,plot_mode);
%% Turning probability depending on relative lateral neighbor distance
% |idSocial_dynamics_TurningReaction| calculates the 
% probability of a focal individual turning left/right 
% given a certain lateral position of the neighbor. This is 
% visualized by mapping the probability onto a discrete
% spatial interval corresponding to the lateral position of 
% the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span for a
% moving average filter to smooth the resulting plot.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = true;
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_TurningReactionProbability 
options.edges=-11.25:.5:11.25; 
options.spacing=3;
% Execute function
data=idSocial_dynamics_TurningReactionProbability(data,options,plot_mode);
%% Forward acceleration strength depending on relative neighbor front/back distance 
% |idSocial_dynamics_AccelerationReaction| calculates the 
% average strength of acceleration/deceleration of a focal
% individual depending on the distance of the neighbor in 
% front/behind of the focal. This is visualized by mapping the 
% focal's acceleration onto a discrete spatial interval 
% corresponding to the position of the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span for a
% moving average filter to smooth the resulting plot.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = true;
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_AccelerationReaction
options.edges=-11.25:.5:11.25; 
options.spacing=3;
% Execute function
data=idSocial_dynamics_AccelerationReaction(data,options,plot_mode);
%% Probability of forward acceleration depending on relative neighbor front/back distance 
% |idSocial_dynamics_AccelerationReactionProbability| calculates the 
% probability of a focal individual accelerating/decelerating 
% given distance of the neighbor in front/behind of the focal. 
% This is visualized by mapping the probability onto a 
% discrete spatial interval corresponding to the position of
% the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span for a
% moving average filter to smooth the resulting plot.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = true;
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_AccelerationReactionProbability
options.edges=-11.25:.5:11.25; 
options.spacing=3;
% Execute function
data=idSocial_dynamics_AccelerationReactionProbability(data,options,plot_mode);
%% Map of probabilities for relative neighbor positions
% |idSocial_dynamics_PositionMap| depicts the probabiliy of
% finding a neighbor individual at a certain position in the
% coordinate system defined by the focal's position
% (origin) and movement direction (y-axis).
% Neighbor positions are binned into spatial intervals with
% edges given by |options.edges|. 
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. ||options.spacing| sets the span of
% bin overlap to smooth the resulting plot.
%
% <<focal_system_cartesian.png>>

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_PositionMap
options.edges={-10.5:3.5:10.5; -10.5:3.5:10.5}; 
options.spacing=1;
options.order_neighbors = true;
% Choose what to plot 
plot_mode.filter= {'Neighbor' 1; 'Neighbor' 2; 'Neighbor' 3};
% Execute function
data=idSocial_dynamics_PositionMap(data,options,plot_mode);
%% Map of turning strength depending on relative neighbor position
% |idSocial_dynamics_TurningMap| calculates the average 
% strength of left/right turns of a focal individual
% depending on the position of the neighbor in the
% coordinate system defined by the focal's position
% (origin) and movement direction (y-axis). This is 
% visualized by mapping the focal's turning acceleration 
% onto a discrete spatial interval corresponding to the
% relative position of the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span of
% bin overlap to smooth the resulting plot.
%
% <<focal_system_cartesian.png>>

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_TurningMap
options.edges={-10.5:3.5:10.5; -10.5:3.5:10.5}; 
options.spacing=1;
% Choose what to plot 
plot_mode.filter= {'Trial' 1; 'Trial' 2; 'Trial' 3};
% Execute function
data=idSocial_dynamics_TurningMap(data,options,plot_mode);
%% Map of turning probability depending on relative neighbor position
% |idSocial_dynamics_TurningProbabilityMap| calculates the 
% probability of a focal individual turning left/right 
% given a certain position of the neighbor in the
% coordinate system defined by the focal's position
% (origin) and movement direction (y-axis). This is 
% visualized by mapping the probability onto discrete 
% spatial intervals corresponding to the relative position
% of the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span of
% bin overlap to smooth the resulting plot.
%
% <<focal_system_cartesian.png>>

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_TurningMap
options.edges={-10.5:3.5:10.5; -10.5:3.5:10.5};  
options.spacing=1;
% Choose what to plot 
plot_mode.filter= {'Trial' 1; 'Trial' 2; 'Trial' 3};
% Execute function
data=idSocial_dynamics_TurningProbabilityMap(data,options,plot_mode);
%% Map of forward acceleration probability depending on relative neighbor position
% |idSocial_dynamics_AccelerationMap| calculates the 
% probability of a focal individual accelerating/decelerating 
% given the position of the neighbor in the
% coordinate system defined by the focal's position
% (origin) and movement direction (y-axis). This is 
% visualized by mapping the probability onto discrete
% spatial intervals corresponding to the relative position 
% of the neighbor.
% _Additional options_: |options.edges| sets the edges of the
% spatial intervals. |options.spacing| sets the span of
% bin overlap to smooth the resulting plot.
%
% <<focal_system_cartesian.png>>

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min = [];
% Extra options for idSocial_dynamics_AccelerationMap
options.edges={-10.5:3.5:10.5; -10.5:3.5:10.5};  
options.spacing=1;
% Choose what to plot 
plot_mode.filter= {'Trial' 1; 'Trial' 2; 'Trial' 3};
% Execute function
data=idSocial_dynamics_AccelerationMap(data,options,plot_mode);
%% 3.4 HIERARCHIES - Aggression and Territoriality 
%% Aggression: Chasing behavior 
% |idSocial_hierarchies_Chasing| calculates the ratio 
% of frames in which the focal individual chases a neighbor
% during time intervals of width 
% |options.timeintervals_in_min|.
% A chase is detected when the focal's speed is greater than
% an adjustable lower limit 
% |options.velocity_bottom_limit_bl_per_s| (in units of Body
% Length/second), the neighbor it is chasing is located
% within an adjustable angle |options.front_angle_deg| in 
% front of the focal, and the distance 
% |options.distance_upper_limit_bl| between focal and 
% neighbor is less than an adjustable limit distance. 
% _Additional options_: (see above) 
% |options.front_angle_deg|, 
% |options.velocity_bottom_limit_bl_per_s|,
% |options.distance_upper_limit_bl|.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=2;
% Extra options for idSocial_hierarchies_Chasing
options.velocity_bottom_limit_bl_per_s = .1; 
options.distance_upper_limit_bl = 3; 
options.front_angle_deg = 30; 
% Choose what to plot
plot_mode.xaxis={'Focal'};
plot_mode.legendstring={'Group'};
% Execute function
data=idSocial_hierarchies_Chasing(data,options,plot_mode);
%% Aggression: Hierarchies from chasing behavior 
% |idSocial_hierarchies_ChasingGraph| draws a graph 
% representing the hierarchy calculated from the ratio 
% of frames in which the focal individual chases a neighbor
% during time intervals of width 
% |options.timeintervals_in_min|.
% In the resulting graph, the the arrows point from focal to
% neighbor if (# of chases(focal->neighbor)) > 
% (# of chases(neighbor->focal)), and the vertical position
% of node for focal i depends on the mean difference 
% (# of chases focal i -> neighbor j) - 
% (# of chases neighbor j -> focal i) for all 
% (M-1)-neighbors, where M is the total number of 
% individuals.
% A chase is detected when the focal's speed is greater than
% an adjustable lower limit 
% |options.velocity_bottom_limit_bl_per_s| (in units of Body
% Length/second), the neighbor it is chasing is located
% within an adjustable angle |options.front_angle_deg| in 
% front of the focal, and the distance 
% |options.distance_upper_limit_bl| between focal and 
% neighbor is less than an adjustable limit distance.
% _Additional options_: (see above) 
% |options.front_angle_deg|, 
% |options.velocity_bottom_limit_bl_per_s|,
% |options.distance_upper_limit_bl|.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=[];
% Extra options for idSocial_hierarchies_ChasingGraph
options.velocity_bottom_limit_bl_per_s = .1; 
options.distance_upper_limit_bl = 3; 
options.front_angle_deg = 30; 
% Choose what to plot
plot_mode.autocombs={'Trial'};
% Execute function
data=idSocial_hierarchies_ChasingGraph(data,options,plot_mode);
%% Aggression: Symmetry of Chasing behavior
% |idSocial_hierarchies_ChasingSymmetry| calculates the ratio 
% of frames in which the focal individual chases a neighbor
% to his right over the total number of chases during time
% intervals of width |options.timeintervals_in_min|.
% A chase is detected when the focal's speed is greater than
% an adjustable lower limit 
% |options.velocity_bottom_limit_bl_per_s| (in units of Body
% Length/second), the neighbor it is chasing is located
% within an adjustable angle |options.front_angle_deg| in 
% front of the focal, and the distance 
% |options.distance_upper_limit_bl| between focal and 
% neighbor is less than an adjustable limit distance. 
% _Additional options_: (see above) 
% |options.front_angle_deg|, 
% |options.velocity_bottom_limit_bl_per_s|,
% |options.distance_upper_limit_bl|.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=4;
% Extra options for idSocial_hierarchies_ChasingSymmetry
options.velocity_bottom_limit_bl_per_s = .1; 
options.distance_upper_limit_bl = 4; 
options.front_angle_deg = 30; 
% Choose what to plot
plot_mode.xaxis={'Time'};
plot_mode.legendstring={'Group'};
plot_mode.filter={'Trial' 3};
% Execute function
data=idSocial_hierarchies_ChasingSymmetry(data,options,plot_mode);
%% Territories
% |idSocial_hierarchies_Territories| visualizes the 
% spatial distribution of individuals in the arena. Regions
% with higher probability of finding an individual are more
% opaque. Borders mark the borders of regions where 
% individuals spend a certain ratio of time, adjustable by 
% the input parameter % |options.discard_bottom_percent|, 
% which sets the percentage of time which is to be 
% neglected. Thus, borders mark regions where individuals 
% spend 100-|options.discard_bottom_percent| percent of 
% their time.
% _Additional options_: |options.discard_bottom_percent|
% (see above), |options.no_bins| is a vector with two 
% elements which determines the resolution of the 
% map in x- and y-direction, or a scalar value, in which case 
% the same resolution is assumed in both directions.
% |options.interpolation_factor| sets the interpolation
% factor for smoother maps.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=5;
% Extra options for idSocial_hierarchies_Territories
options.discard_bottom_percent=9;
options.no_bins=[];
options.interpolation_factor=4;
% Choose what to plot
plot_mode.subplots = {'Trials'};
plot_mode.autocombs = {'Focal'};
plot_mode.legendstring = {'Indiv. 1'; 'Indiv. 2'; 'Indiv. 3'; 'Indiv. 4'; 'Indiv. 5';};
% Execute function
data=idSocial_hierarchies_Territories(data,options,plot_mode);
%% Territory Overlap
% |idSocial_hierarchies_TerritoryOverlap| shows the 
% spatial overlap between pairs of territories. It first 
% counts the occurences of an individual passing 
% through a position of the arena and normalizes the result
% for each position by the total number of counts (the 
% number of frames of the video). Thus, for each position it 
% gets the probability of finding an individual at a certain
% position during the time of the video. 
% To calculate the overlap between pairs of individuals, the
% square root of the product of the two results
% corresponding to the two individuals is calculated for
% each position. Finally, the sum over all positions is 
% taken as a measure for the overlap.
% _Additional options_ are the same as for 
% |idSocial_hierarchies_Territories|, since both measures 
% are calculated simultaneously: |options.no_bins| is a 
% vector with two elements which determines the resolution 
% of the grid which is used to define finite positions in 
% the arena in x- and y-direction, or a scalar value, in
% which case the same resolution is assumed in both 
% directions.
% |options.interpolation_factor| sets the interpolation
% factor for smoother maps.
% |options.discard_bottom_percent| is needed by 
% |idSocial_hierarchies_Territory| (see above) and has no
% effect on |idSocial_hierarchies_TerritoryOverlap|.

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=5;
% Extra options for idSocial_hierarchies_TerritoryOverlap
options.discard_bottom_percent=9;
options.no_bins=[6 6];
options.interpolation_factor=3;
%Choose what to plot
plot_mode.autocombs = {'Focal'};
plot_mode.xaxis={'Time'};
% Execute function
data=idSocial_hierarchies_TerritoryOverlap(data,options,plot_mode);
%% 3.5 NETWORKS - Leader-Follower Relationship
%% Movement correlation: Hierarchy from delays in leading and following
% |idSocial_networks_movementCorrelationDelayGraph|
% calculates a leader-follower hierarchy from the delay with
% which an individual copies the movement direction of
% another. The hierarchy for the whole group of individuals 
% is calculated from the pairwise interactions as follows:
% The vertical position of an individual is determined by
% the mean delay over all neighbors.
% The arrows depict the pairwise relation between two
% individuals and point from the leader towards the neighbor
% which is copying the leader's direction.  
% _Additional options_: |options.delay_seconds| sets the
% maximum delay time in seconds for which copying behavior 
% is counted as leading/following. 
% |options.max_distance_pxl| sets the maximum 
% inter-individual distance in pixel for which copying 
% behavior is counted as leading/following. 
% |delay_correlation_lim| sets the bottom limit of the
% correlation coefficient for which copying behavior is 
% counted as leading/following. The correlation coefficient
% reflects the angle between movement directions of two
% individuals (the correlation coefficient is proportional 
% to the cosine of mutual angle). It corresponds to the 
% minimum angle within the range of time defined by 
% |options.delay_seconds|.
%
% <<movement_correlation.png>>

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=[];
% Extra options for idSocial_networks_movementCorrelationDelayGraph
options.delay_seconds = 2;
options.delay_correlation_lim = .0;
%Choose what to plot
plot_mode.autocombs={'Trial'};
% Execute function
data=idSocial_networks_movementCorrelationDelayGraph(data,options,plot_mode);
%% Movement correlation: Hierarchy from sign of delays in leading and following

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=[];
% Extra options for idSocial_hierarchies_Territories
options.delay_seconds = 2;
options.delay_correlation_lim = .0;
%Choose what to plot
plot_mode.autocombs={'Trial'};
% Execute function
data=idSocial_networks_movementCorrelationSignDelayGraph(data,options,plot_mode);
%% DELETE? Movement correlation: Hierarchy from correlation coefficient in leading and following
% |idSocial_networks_movementCorrelationDelayGraph|
% calculates a leader-follower hierarchy from the 
% correlation coefficient which reflects the angle between 
% the directions of movement of two individuals (the
% correlation coefficient equals the cosine of 
% the angle between the directions of movement). The
% correlation coefficient is calculated as the cosine of the 
% minimum angle between the directions of movement of two 
% individuals within a range of time defined by 
% |options.delay_seconds|.
% The hierarchy for the whole group of individuals 
% is calculated from the pairwise interactions as follows:
% The vertical position of an individual is determined by
% the mean correlation coefficient over all neighbors.
% The arrows depict the pairwise relation between two
% individuals and point from the leader towards the neighbor
% which is copying the leader's direction.  
% _Additional options_: |options.delay_seconds| sets the
% maximum delay time in seconds for which copying behavior 
% is counted as leading/following. 
% |options.max_distance_pxl| sets the maximum 
% inter-individual distance in pixel for which copying 
% behavior is counted as leading/following. 
% |delay_correlation_lim| sets the bottom limit of the
% correlation coefficient for which copying behavior is 
% counted as leading/following. The correlation coefficient
% reflects the angle between movement directions of two
% % individuals (the correlation coefficient is proportional 
% to the cosine of mutual angle). It corresponds to the 
% minimum angle within the range of time defined by 
% |options.delay_seconds|.


% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=[];
% Extra options for idSocial_hierarchies_Territories
options.delay_seconds = 2;
options.delay_correlation_lim = .0;
%Choose what to plot
plot_mode.autocombs={'Trial'};
% Execute function
data=idSocial_networks_movementCorrelationGraph(data,options,plot_mode);
%% Movement correlation: Delays in leading and following
% |idSocial_networks_movementCorrelationDelayGraph|

% Clear previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.timeintervals_in_min=[];
% Extra options for idSocial_hierarchies_Territories
options.delay_seconds = 2;
options.delay_correlation_lim = .0;
%Choose what to plot
plot_mode.xaxis={'Time'};
plot_mode.autocombs={'Focal'};
plot_mode.subplots={'Trial'};
% Execute function
data=idSocial_networks_movementCorrelationDelay(data,options,plot_mode);
%% 4. SUPPELMENTARY INFORMATION
%% 4.1 DATA PREPROCESSING: Trajectory filtering and smoothing
%% 4.2 RANDOM CONTROLS: Randomization of temporal order of trajectories
%% 4.3 PLOT FUNCTIONS
