function [filter_all filter_trajectory filter_velocity filter_acceleration]=strucdyn_filters(trayectorias,prob,avge,options)
avg=[];
blpxl=[];
framerate=[];
max_dist_in_bl=[];
min_dist_in_bl=[];
max_prob=[];
speedlim_bl_per_s=[];
maxspeed_bl_per_s=[];
neighbour_speedlim_bl_per_s=[];
neighbour_maxspeed_bl_per_s=[];
roi4boundary=[];
% thisfile=[mfilename('fullpath') '.m'];
% system(['copy ' thisfile '+test.m MyBigFat.m']);
if nargin<2 || isempty(prob)
    prob=ones(size(trayectorias,1),size(trayectorias,2));
    avg=1;
end
if nargin<3 || isempty(avge)
    avg=1;
else
    avg=avge;
end
if nargin>=4 && ~isempty(options)
    if isfield(options,'avg')
        avg=options.avg;
    elseif ~isfield(options,'avg') && isfield(options,'avg_preparedata')
        avg=options.avg_preparedata;
    end
    if isfield(options,'blpxl')
        blpxl=options.blpxl;
    end
    if isfield(options,'framerate')
        framerate=options.framerate;
    end
    if isfield(options,'max_dist_in_bl')
        max_dist_in_bl=options.max_dist_in_bl;
    end
    if isfield(options,'dist_limits_in_bl')
        max_dist_in_bl=options.dist_limits_in_bl(2);
        min_dist_in_bl=options.dist_limits_in_bl(1);
    end
    if isfield(options,'max_prob')
        max_prob=options.min_prob;
    end
    if isfield(options,'speedlim_bl_per_s')
        speedlim_bl_per_s=options.speedlim_bl_per_s;
    end
    if isfield(options,'neighbour_speedlim_bl_per_s')
        neighbour_speedlim_bl_per_s=options.neighbour_speedlim_bl_per_s;
    end
    
    if isfield(options,'maxspeed_bl_per_s')
        maxspeed_bl_per_s=options.maxspeed_bl_per_s;
    end
    if isfield(options,'neighbour_maxspeed_bl_per_s')
        neighbour_maxspeed_bl_per_s=options.neighbour_maxspeed_bl_per_s;
    end
    if isfield(options,'roi4boundary')
        roi4boundary=options.roi4boundary;
    end
    
end
if isempty(avg); avg=1; end;
if isempty(blpxl); blpxl=90;  end;
if isempty(framerate); framerate=30; end;
if isempty(max_dist_in_bl); max_dist_in_bl=inf; end;
if isempty(min_dist_in_bl); min_dist_in_bl=0; end;
if isempty(max_prob); max_prob=.5; end;
if isempty(speedlim_bl_per_s); speedlim_bl_per_s=0; end;
if isempty(neighbour_speedlim_bl_per_s); neighbour_speedlim_bl_per_s=0; end;
if isempty(maxspeed_bl_per_s); maxspeed_bl_per_s=Inf; end;
if isempty(neighbour_maxspeed_bl_per_s); neighbour_maxspeed_bl_per_s=Inf; end;
if isempty(roi4boundary); roi4boundary=[-Inf -Inf Inf Inf]; end;




% Parameters for filtering etc.


% For boundary_filter [SO FAR NOT APPLIED, ADJUST FOR REAL DATA IF
% NECESSARY]:

%     blcm=2.5;                       % Body length in cm;


% For velocity filter. Framerate and max. allowed velocity in (body length)/s:
% framerate=30;
% speedlim_bl_per_s=0;

% For distance filter. Max. allowed distance in body length:
% max_dist_in_bl=2.5;


% max_prob=.4; % All data with a probability for being wrong lower than this will be filtered (probability is calculated by Alfonso's tracking and stored within the trayectorias.mat)

% blpxl=150;                      % Body length in pixels, 20120220_hungryfish;
%     blpxl=90;                      % Body length in pixels;
filter_all=cell(2,1);
for fish=[1 2; 2 1]
    
    prob_filter=false(size(trayectorias,1),1);
    velocity_filter=false(1,size(trayectorias,1));
    neighbour_velocity_filter=false(1,size(trayectorias,1));
    x=[trayectorias(:,fish(1),1) trayectorias(:,fish(2),1)]';
    y=[trayectorias(:,fish(1),2) trayectorias(:,fish(2),2)]';
    %     prob=prob(1:9000);
    
    xavg=tsmovavg(x,'s',avg); % Smoothing trajectories with moving window average, according to Katz11 (window size 10)
    yavg=tsmovavg(y,'s',avg);
    
    isfinite_filter = isfinite(x(1,:)) & isfinite(y(1,:));
    isfinite_filter_nb = isfinite(x(2,:)) & isfinite(y(2,:));
    isfinite_total=isfinite_filter & isfinite_filter_nb;
    
    vx=diff(xavg,1,2);
    vy=diff(yavg,1,2);
    
    acceleration_x=diff(vx,1,2);
    acceleration_y=diff(vy,1,2);
    
    
    %% Filtering experimental data
    
    
    % Filter trajectories 1: Ommit trajectories for which the tracking has
    % determined a low probability:
    
    prob_filter((prob(:,fish(2))>max_prob)&(prob(:,fish(1))>max_prob)&~isnan(prob(:,fish(2)))&~isnan(prob(:,fish(1))))=true;
    
    % Filter trajectories 2: Ommit frames in which any fish is closer than bl
    % body lengths to boundary (Katz11: 2.5 BL)
    scaling=1/blpxl;                % Set to (length of fish in BL)/(length of fish in pixels)
    marg=2.5/scaling;               % Margin in pixels
    %     boundary_filter=((x(1,:)>roi4boundary(1))&(x(2,:)>roi4boundary(2))&(x(1,:)<roi4boundary(3))&(x(2,:)<roi4boundary(4)));
    boundary_filter=((x(1,:)>roi4boundary(1))&(y(1,:)>roi4boundary(2))&(x(1,:)<roi4boundary(3))&(y(1,:)<roi4boundary(4))) & ...
        ((x(2,:)>roi4boundary(1))&(y(2,:)>roi4boundary(2))&(x(2,:)<roi4boundary(3))&(y(2,:)<roi4boundary(4))) & isfinite_total;
    
    % Filter trajectories 3: Ommit trajectories where velocity is low, <fs
    % (Katz11: 0.5 BL/s)
    neighbour_speedlim_pxl_per_frame=neighbour_speedlim_bl_per_s*blpxl/framerate;
    neighbour_maxspeed_pxl_per_frame=neighbour_maxspeed_bl_per_s*blpxl/framerate;
    speedlim_pxl_per_frame=speedlim_bl_per_s*blpxl/framerate;
    maxspeed_pxl_per_frame=maxspeed_bl_per_s*blpxl/framerate;
    focal_vel=(vx(1,:).^2+vy(1,:).^2).^.5;
    neighbour_vel=(vx(2,:).^2+vy(2,:).^2).^.5;
    %     velocity_filter=[(focal_vel>speedlim_pxl_per_frame).*(~isnan(focal_vel))  0];
    velocity_filter((focal_vel>speedlim_pxl_per_frame)&(~isnan(focal_vel))& focal_vel<maxspeed_pxl_per_frame )=true;
    %     velocity_filter((focal_vel<maxspeed_pxl_per_frame)&(~isnan(focal_vel)))=true;
    
    %     neighbour_velocity_filter(( neighbour_vel> neighbour_speedlim_pxl_per_frame)&(~isnan( neighbour_vel)))=true;
    if neighbour_maxspeed_pxl_per_frame==inf && neighbour_speedlim_pxl_per_frame<=0 % if neighbour velocity does not matter
        neighbour_velocity_filter=true(1,size(trayectorias,1));
    else
        neighbour_velocity_filter(( neighbour_vel< neighbour_maxspeed_pxl_per_frame)&...
            (~isnan( neighbour_vel)) & ...
            ( neighbour_vel> neighbour_speedlim_pxl_per_frame)...
            )=true;
    end
    % Filter trajectories 4: Ommit trajectories for which distance between
    % individuals is big
    max_dist_in_pxl=max_dist_in_bl*blpxl;
    min_dist_in_pxl=min_dist_in_bl*blpxl;
    distance_focal_to_neighbour=((xavg(1,:)-xavg(2,:)).^2+(yavg(1,:)-yavg(2,:)).^2).^.5;
    distance_filter=distance_focal_to_neighbour<max_dist_in_pxl & distance_focal_to_neighbour>min_dist_in_pxl;
    
    %     velocity_filter=ones(size(velocity_filter));
    %     velocity_filter(end)=0;
    
    % Combine filters: For N trajectory points, there are N-1 data points for
    % velocity and N-2 for acceleration. If the kth trajectory point is filtered, so is the kth of velocity and acceleration. In other words: The same filter is applied to all three vectors in the end.
    all_filter_trajectory=logical(distance_filter&prob_filter'); %
    %     filter_trajectory{fish(1)}=all_filter_trajectory;
    % In terms of velocity, if the kth element of the trajectory is filtered,
    % so is the kth and (k-1)-element of velocity, whose calculation involves the filtered kth trajectory point. Similarly, if velocity k is filtered because of exceeding the speed limit, (k+1) of the trajectory is filtered.
    all_filter_velocity=all_filter_trajectory & [all_filter_trajectory(2:end) true]& velocity_filter & [true velocity_filter(1:end-1)];
    %     filter_velocity{fish(1)}=all_filter_velocity;
    all_filter_neighbour_velocity=all_filter_velocity & neighbour_velocity_filter & [true neighbour_velocity_filter(1:end-1)];
    % Again, having filtered one velocity data point entails filtering of two
    % acceleration data points which are calculated from it.
    all_filter_acceleration=all_filter_neighbour_velocity & [all_filter_neighbour_velocity(2:end) true] & isfinite_total;
    %     filter_acceleration{fish(1)}=all_filter_acceleration;
    
    % All filters are the same in the end:
    disp(['[' mfilename '] ' 'ATTENTION: Frame rate is set to ' num2str(framerate) ' fps. Please check if this is correct.'])
    
    disp(['[' mfilename '] ' 'ATTENTION: Body length in pixels is set to ' num2str(blpxl) ' pxl/Bl. Please check if this is correct.'])
    disp(['[' mfilename '] ' 'strucdyn_filters is averaging over ' num2str(avg) ' frames.'])
    
    disp(['[' mfilename '] ' 'Number of frames: ' num2str(size(x,2)) ' = ' num2str(size(x,2)/framerate/60) ' min'])
    disp(['[' mfilename '] ' 'Number of frames (not NaN or Inf): ' num2str(sum(isfinite_total)) ' = ' num2str(sum(isfinite_total)/framerate/60) ' min'])
    %             if filter==1
    all_filter_trajectory=all_filter_acceleration;
    all_filter_velocity=all_filter_acceleration;
    filter_all{fish(1)}=all_filter_acceleration;
    
    disp(['[' mfilename '] ' 'Number of good frames after filtering: ' num2str(sum(all_filter_acceleration & isfinite_total)) '=' num2str(sum(all_filter_acceleration(isfinite_total))/sum(isfinite_total)*100) '% = ' num2str(sum(all_filter_acceleration(isfinite_total))/framerate/60) ' min'])
    disp(['[' mfilename '] ' 'Passed probability filter (min. probability=' num2str(max_prob) '): ' num2str(sum(prob_filter(isfinite_filter))) '=' num2str(sum(prob_filter(isfinite_filter))/sum(isfinite_filter)*100) '%'])
    disp(['[' mfilename '] ' 'Passed distance filter (distance=[' num2str( min_dist_in_bl ) ' ' num2str( max_dist_in_bl ) '] BL): ' num2str(sum(distance_filter(isfinite_total))) '=' num2str(sum(distance_filter(isfinite_total))/sum(isfinite_total)*100) '%'])
    disp(['[' mfilename '] ' 'Passed velocity filter (vel=' num2str(speedlim_bl_per_s) '-'  num2str(maxspeed_bl_per_s) ' BL/s): '  num2str(sum(velocity_filter(isfinite_total))) '=' num2str(sum(velocity_filter(isfinite_total))/sum(isfinite_total)*100) '%'])
    disp(['[' mfilename '] ' 'Passed neighbour velocity filter (vel=' num2str(neighbour_speedlim_bl_per_s) '-'  num2str(neighbour_maxspeed_bl_per_s) ' BL/s): '  num2str(sum(neighbour_velocity_filter(isfinite_total))) '=' num2str(sum(neighbour_velocity_filter(isfinite_total))/sum(isfinite_total)*100) '%'])
    disp(['[' mfilename '] ' 'Passed boundary filter (roi=[' num2str(roi4boundary/blpxl) '] BL): '  num2str(sum(boundary_filter(isfinite_total))) '=' num2str(sum(boundary_filter(isfinite_total))/sum(isfinite_total)*100) '%'])
    
    %             else
    %         %         keyboard
    %                 all_filter_acceleration=true(size(all_filter_acceleration));
    %                 all_filter_acceleration(end-1:end)=false(1,2);
    %                 all_filter_trajectory=all_filter_acceleration;
    %                 all_filter_velocity=all_filter_acceleration;
    %         %         keyboard
    %
    %             end
    %     if filter==1
    
    
    %     end
    
    
    
    
end
