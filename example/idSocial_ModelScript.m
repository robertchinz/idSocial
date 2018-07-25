%% idSocial example script

root_dir = 'C:\Users\hinz\Documents\Syncplicity Folders\Syncplicity\git_idSocial\';
output_directory = [root_dir 'example\Development\results\'];
temporary_dir = [root_dir '~temp'];
trajectory_dir=[root_dir 'example\Development\N2\day21'];


% Create list of trajectory files
trajectory_location={...
    [trajectory_dir '\trial45_8mm\trajectories.mat'];...
    [trajectory_dir '\trial46_8mm\trajectories.mat'];...
    [trajectory_dir '\trial47_7mm\trajectories.mat'];...
    [trajectory_dir '\trial48_6mm\trajectories.mat'];...
};


% Display options for pre-processing by calling idSocial_loadData without 
% input parameters
disp(idSocial_loadData)

% Set parameters for data preprocessing
clear options 
options.project_path = output_directory;  % Path where to save results.   
options.blpxl = 50; % Body length of animals in pixel.
options.framerate = 30; % Frame rate of the original video.
options.smooth_method = 'moving'; % Method for trajectory smoothing.
options.smooth_degree = 30; % Window width for moving average smoothing.
options.filter_focal_speedlimits_bl_per_s = [.1 inf]; % In all following 
% analysis, data of a focal individual will only be cosidered if it is 
% speed exceeds 0.1 BL.
options.temp_savepath = temporary_dir; % Path for temporary files (idSocial 
% sometimes needs to buffer data when there is not sufficient memory)
options.random_data = true; % Control data from randomization will be 
% calculated.
options.no_RandomNeighbors = 4; % Number of virtual neighbors obtained 
% from randomization.

% Run idSocial_loadData:
data=idSocial_loadData(trajectory_location, [], options);

%% Speed Distributions
% Reset previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.edges=0:.1:4;
% Run analysis
data=idSocial_kinematics_SpeedDistribution(data,options,plot_mode);
% Display results
disp(data.kinematics_SpeedDistribution.results)
figure('Color','w','Units','normalized','Position',[.1 .1 .5 .5]);
plot(options.edges(1:end-1),data.kinematics_SpeedDistribution.results{3,2})
xlabel('Speed [BL/s]')
ylabel('Counts')
%% Acceleration Distributions
% Reset previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.edges=0:.25:7;
% Run analysis
data=idSocial_kinematics_AccelerationDistribution(data,options,plot_mode);
% Display results
disp(data.kinematics_AccelerationDistribution.results)
figure('Color','w','Units','normalized','Position',[.1 .1 .5 .5]);
plot(options.edges(1:end-1),data.kinematics_AccelerationDistribution.results{3,2})
xlabel('Acc. [BL/s^2]')
ylabel('Counts')
%% Distance Distributions 
% using kernel density smoothing 
% Reset previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = true;
options.edges = 0:.25:21;
options.method='ksdensity_epanechnikov';
options.kds_bandwidth = 2;
plot_mode.statistics = {'hist+1stMode'};
plot_mode.pool_data = true;
% Run analysis
data=idSocial_interaction_DistanceDistribution(data,options,plot_mode);
disp(data.interaction_DistanceDistribution.results)
figure('Color','w','Units','normalized','Position',[.1 .1 .5 .5]); 
plot(options.edges(1:end-1),data.interaction_DistanceDistribution.results{3,2})
xlabel('Distance [BL]')
ylabel('Counts')

%% Distance Distributions 
% using kernel density smoothing 
% Reset previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.random_data = false;
options.edges = {-6:1:6;-6:1:6};
options.spacing =5;
% Run analysis
data=idSocial_dynamics_PositionMap(data,options,plot_mode);


