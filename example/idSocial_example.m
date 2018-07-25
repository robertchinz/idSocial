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

%% 1.6 Example: Data preprocessing  

clear options           % This delet    es any existing variable
                        % 'options' 
options.start_min = 0;  % Only frames from a later time (in
                        % minutes) will be analysed.
options.end_min = 180;  % Frames from later times (in
                        % minutes) will be discarded.
options.project_path = ''; % If this is empty, no 
                        % figures will be saved. If it is a
                        % valid path, figures will be saved 
                        % to path\figures and 
                        % path\figures_latex.
options.blpxl = 50;     % Body length of the individuals in 
                        % pixels.
                        % If idTracker is used to generate 
                        % the trajectory-files, and idSocial
                        % finds the corresponding 
                        % 'datosegm'-files, body length 
                        % will be taken directly from 
                        % 'datosegm.mat' for each trial.
                        % blpxl is used to present results 
                        % in units of body length.
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
options.interpolate_trajectories = false; % If set 'true', 
                        % missing data points will be 
                        % interpolated. NOTE: Only use 
                        % interpolation if you are sure it 
                        % will not create 'fictional' data.
options.interpolation_mode = 'linear'; % Set interpolation
                        % method. Possible options:
                        % 'linear' (default)
                        % 'spline' 
options.smooth_method = 'moving';% Set smoothing method. Is 
                        % left empty no smoothing will be 
                        % performed. Possible options:
                        % 'moving' (default)
                        % 'lowess'
                        % 'loess'
                        % 'rlowess'
                        % 'rloess'
                        % 'spline'
                        % See also: MATLAB(R) 'smooth'
options.smooth_degree = 5; % Size of the moving window used for
                        % smoothing, i.e., number of 
                        % neighboring data points used for
                        % calculation of adjusted central 
                        % data point.
options.filter_minProbIdentityAssignment = .7;
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
options.group_names = {'Group 1'; 'Group 2';}; 
                        % If trials  are sorted into 
                        % different groups, names for each 
                        % group can be set (e.g., 
                        % {'Test'; 'Control'};
options.temp_savepath = 'D:\~idSocialTemp\'; % idSocial will
                        % save results temporarily at the
                        % given location. If left empty, 
                        % you will be asked to choose the 
                        % location when executing any 
                        % idSocial method for the first time 
                        % with the given data.
options.random_data = true; % Create artifical trajectories 
                        % by randomizing the temporal order
                        % of the original trajectory data. 
                        % These artificial trajectories
                        % can later be used to compare the
                        % original data to the randomized
                        % controls. 
                        % See also options.random_data 
                        % below.

[data, ~, info]=idSocial_loadData(trajectory_location, datosegm_location, options);
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
options.timeintervals_in_min=2.5;
options.random_data = false;
options.filter_focal_speedlimits_bl_per_s = [0 10];
options.random_data = true;
% Choose what to plot
plot_mode.xaxis={'Time'};
plot_mode.legendstring={'Group 1';};
plot_mode.filter={'Group' 1;};
plot_mode.filter={'Trial' 1; 'Trial' 2; 'Trial' 3};

% Execute function
data=idSocial_interaction_Distance(data,options,plot_mode);