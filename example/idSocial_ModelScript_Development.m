%% idSocial example script

root_dir = pwd; %'C:\Users\hinz\Documents\Syncplicity Folders\Syncplicity\git_idSocial';
output_directory = [root_dir '\example\Development\results\'];
temporary_dir = [root_dir '~temp'];
trajectory_dir=[root_dir '\example\Development\N2'];

trajectory_location=[];
% Create list of trajectory files
trajectory_location{1}{7}={...
    [trajectory_dir '\day7\trial9_4mm\trajectories.mat'];...
    [trajectory_dir '\day7\trial10_4mm\trajectories.mat'];...
    [trajectory_dir '\day7\trial11_4mm\trajectories.mat'];...
    [trajectory_dir '\day7\trial12_4mm\trajectories.mat'];...
};

trajectory_location{1}{15}={...
    [trajectory_dir '\day15\trial13_5mm\trajectories.mat'];...
    [trajectory_dir '\day15\trial14_5mm\trajectories.mat'];...
    [trajectory_dir '\day15\trial16_4mm\trajectories.mat'];...
    [trajectory_dir '\day15\trial17_5mm\trajectories.mat'];...
};

trajectory_location{1}{24}={...
    [trajectory_dir '\day24\trial50_7mm\trajectories.mat'];...
    [trajectory_dir '\day24\trial51_7mm\trajectories.mat'];...
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
data=idSocial_kinematics_Speed(data,options,plot_mode);
% Display results
disp(data.kinematics_SpeedDistribution.results)
figure('Color','w','Units','normalized','Position',[.1 .1 .5 .5]);
plot(options.edges(1:end-1),data.kinematics_Speed.results{3,2})
xlabel('Speed [BL/s]')
ylabel('Counts')
%% Acceleration Distributions
% Reset previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.edges=0:.25:7;
% Run analysis
data=idSocial_kinematics_Acceleration(data,options,plot_mode);
% Display results
disp(data.kinematics_AccelerationDistribution.results)
figure('Color','w','Units','normalized','Position',[.1 .1 .5 .5]);
plot(options.edges(1:end-1),data.kinematics_Acceleration.results{3,2})
xlabel('Acc. [BL/s^2]')
ylabel('Counts')
%% Distance Distributions 
% using kernel density smoothing 
% Reset previous options and plot_mode
options = []; plot_mode = [];
% Set new options
options.edges = 0:.25:21;
options.method='hist';
plot_mode.filter = {'Day' 7; 'Day' 15; 'Day' 24};
% Run analysis
data=idSocial_interaction_Distance(data,options,plot_mode);
disp(data.interaction_Distance.results)
figure('Color','w','Units','normalized','Position',[.1 .1 .5 .5]); 
plot(options.edges(1:end-1),data.interaction_Distance.results{3,2})
hold on
plot(options.edges(1:end-1),data.interaction_Distance.results{3,3})
plot(options.edges(1:end-1),data.interaction_Distance.results{3,4})

xlabel('Distance [BL]')
ylabel('Counts')

%% Distance vs. age, compare with randomized controls
% using kernel density smoothing 
% Reset previous options and plot_mode

options = []; plot_mode = [];
% Set new options
options.random_data = true;
plot_mode.statistics = {'mean'};
plot_mode.xaxis = {'Subset'};

data=idSocial_interaction_Distance(data,options,plot_mode);

disp(data.interaction_Distance.results)
figure('Color','w','Units','normalized','Position',[.1 .1 .5 .5]); 
% Plot results
plot([7 15 24],data.interaction_Distance.results{3,2}([7 15 24]),'b-o')
hold on
% Plot controls
plot([7 15 24],data.interaction_Distance.results{7,2}([7 15 24]),'k-+')

xlabel('Day')
ylabel('Distance [BL]')
