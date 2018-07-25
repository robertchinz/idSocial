% Load trajectory to workspace
% load('E:\idSocialUI\example\Development\N2\day21\trial45_8mm\trajectories.mat')
load('C:\Users\hinz\Documents\Syncplicity Folders\Syncplicity\git_idSocial\example\Development\N2\day21\trial45_8mm\trajectories.mat')

% Size of variable 'trajectory'
[no_frames, no_fish, no_dimensions] = size(trajectories);

% Plot trajectories for fish 1
figure('Name','Trajectory'); 
plot3(trajectories(:,1,1),trajectories(:,1,2),1:no_frames)
xlabel('X'); ylabel('Y'); zlabel('Time');

% Calculate stuff, for example inter-indiviudal distance

dist = sqrt((trajectories(:,1,1)-trajectories(:,2,1)).^2 + ...
    (trajectories(:,1,2)-trajectories(:,2,2)).^2);

figure('Name','Distance'); 
plot(dist)
xlabel('Time'); ylabel('Distance'); 

%% More examples
% Speed for fish 1 from difference between coordinates in adjacent frames
speed = sqrt(diff(trajectories(:,1,1),[],1).^2 + ...
    diff(trajectories(:,1,2),[],1).^2);

% Mean speed:
nanmean(speed)

% Acceleration for fish 1 from difference between speed in adjacent
% frames
acc = sqrt(diff(trajectories(:,1,1),[],1).^2 + ...
    diff(trajectories(:,1,2),[],1).^2);

% Acceleration median:
nanmedian(acc)
