function [md simplestats]=idSocial_preparedata(varargin)
% function [md simplestats]=strucdyn_preparedata(trayectorias,prob,filters_on,avg,options)
% [md simplestats]=strucdyn_preparedata(trayectorias,prob,filters_on,avg,options)
% 'focal_r',[],'focal_v',[],'focal_a','focal_v_magnitude',[],'focal_a_magnitude',[],'neighbour_r',[],'neighbour_v',[],'neighbour_a',[],...
%                 'distance_focal_neighbour',[],'temporal_distance_focal_neighbour',[],'focal_movement_direction',[],'focal_to_neighbour_direction',[],'focal_to_neighbour_vector',[],...
%                 'neighbour_movement_direction',[],...
%                 'angle_focal_movement_direction_and_y_axis',[],'angle_focal_movement_direction_and_neighbour_movement_direction',[],...
%                 'angle_focal_movement_direction_and_neighbour_position',[],...
%                 'angle_neighbour_movement_direction_and_focal_position',[],...
%                 'angle_focal_movement_direction_t1_and_t2_integrated',[],...
%                 'angle_focal_movement_direction_t1_and_t2',[],...
%                 'angle_focal_movement_direction_t1_and_t2_magnitude',[],...
%                 'angle_neighbour_movement_direction_t1_and_t2_magnitude',[],...
%                 'angle_acceleration_focal_movement_direction_t1_and_t2',[],...
%                 'focal_to_neighbour_vector_rotated',[],...
%                 'focal_to_neighbour_direction_rotated',[],...
%                 'focal_v_rotated',[],'focal_a_rotated',[],...
%                 'neighbour_v_rotated',[],...
%                 'neighbour_v_transformed',[],...
%                 'neighbour_movement_direction_rotated',[],'neighbour_a_rotated',[],...
%                 'neighbour_a_transformed',[],...
%                 'angle_focal_acceleration_and_neighbour_position',[],...
%                 'focal_acceleration_direction',[],...
%                 'angle_focal_acceleration_and_movement_direction',[],...
%                 'angle_focal_movement_direction_and_focal_acceleration',[],...
%                 'focal_position_direction',[],...
%                 'angle_focal_position_and_y_axis',[],...
%                 'focal_distance_to_origin',[],...
%                 'focal_polar_radius',[],...
%                 'focal_polar_angle',[],...
%                 'neighbour_polar_radius',[],...
%                 'neighbour_polar_angle',[],...
%                 'focal_polar_radius_velocity',[],...
%                 'focal_polar_angle_velocity',[],...
%                 'neighbour_polar_radius_velocity',[],...
%                 'neighbour_polar_angle_velocity',[],...
%                 'polar_angle_distance_focal_to_neighbour',[]);

% Getting a whole set of measures from trajectories.
% RCH 04/10/2010


noprobs=false;
[cax,args,nargs] = axescheck(varargin{:});


inpchartest=NaN(1,nargs-1);
for k=1:nargs-1
    inpchartest(k)=ischar(args{k+1});
end
if ~isempty(find(inpchartest,1,'first'))
    charidx=min(find(inpchartest,1,'first'),nargs);
else
    charidx=nargs;
end

tic;
if charidx<5 || (charidx>=5 && isempty(args{5}))
    options=[];
elseif (charidx>=5 && ~isempty(args{5}))
    options=args{5};
end
if ischar(args(1))
    if exist(args(1),'file')==2
        load(args(1));
        trayectorias=idSocial_interpolateTrajectories(trayectorias,interp_mode);
    end
else
    trayectorias=args{1};
end

if charidx<2 || (charidx>=2 && isempty(args{2}))
    prob=ones(size(trayectorias(:,2,1)));
    noprobs=true;
elseif (charidx>=2 && ~isempty(args{2}))
    prob=args{2};
end

if charidx<3 || (charidx>=3 && isempty(args{3}))
    filters_on=true;
elseif (charidx>=3 && ~isempty(args{3}))
    filters_on=args{3};
end
if (charidx<4  &&  (isempty(options) || ~isfield(options,'avg_preparedata'))) || (charidx>=4 && isempty(args{4}))
    avg=3;
elseif (~isempty(options) && isfield(options,'avg_preparedata'))
    avg=options.avg_preparedata;
elseif (charidx>=4 && ~isempty(args{4}))
    avg=args{4};
end
if ~isempty(options) && isfield(options,'interpolation_mode') && ~isempty(options.interpolation_mode)
    interp_mode=options.interpolation_mode;
else
    interp_mode='linear';
end

% if nargs > 1 && ~ischar(args{2})
%     prob = args{2};
%     args = args(3:end); % strip off x and nbins/ctrs
% else
%     prob=zeros(size(args(1)));
%     noprobs=true;
%     args = args(2:end); % strip off x
% end
no_fish=size(trayectorias,2);
args=args(charidx+1:end);
params_temporal_distance=[1 1 1];
% pnames = {'prob','filters_on','avg','options','polar_center','polar_axis','exclusion_radius','params_temporal_distance'};
% dflts =  { prob,     filters_on,       avg,    options, [0 0], [1 0], 0 , 1:no_fish,[1 1 1]};
% [errid,errmsg,prob,filters_on,avg,options,polar_center,polar_axis,exclusion_radius] = internal.stats.getargs(pnames, dflts, args{:});




if isempty(options) || options.interpolate_trajectories
    trayectorias=idSocial_interpolateTrajectories(trayectorias,interp_mode,60);
    disp('Interpolating trajectory directly.')
end

if isempty(options) || (isfield(options,'interpolate_polar_angle_and_radius') && options.interpolate_polar_angle_and_radius)
    % Interpolation of polar angle will be done after angles have been
    % calculated.
end

if ~isempty(options) && (~isfield(options,'focal') || isempty(options.focal))
    options.focal=1:no_fish;
    focal=1:no_fish;
elseif ~isempty(options) && isfield(options,'focal') && ~isempty(options.focal)
    focal=options.focal; %% OVERWRITES EXTRA-INPUT!
elseif isempty(options)
    focal=1:no_fish;
end

if ~isempty(options) && ~isfield(options,'polar_center')
    options.polar_center=[0 0];
    polar_center=[0 0];
elseif ~isempty(options) && isfield(options,'polar_center')
    polar_center=options.polar_center; %% OVERWRITES EXTRA-INPUT!
end
if ~isempty(options) && ~isfield(options,'polar_axis')
    options.polar_axis=[1 0];
    polar_axis=[1 0];
elseif ~isempty(options) && isfield(options,'polar_axis')
    polar_axis=options.polar_axis; %% OVERWRITES EXTRA-INPUT!
end
if ~isempty(options) && ~isfield(options,'exclusion_radius')
    options.exclusion_radius=0;
    exclusion_radius=0;
elseif ~isempty(options) && isfield(options,'exclusion_radius')
    exclusion_radius=options.exclusion_radius; %% OVERWRITES EXTRA-INPUT!
end
if ~isempty(options) && isfield(options,'avg_preparedata')
    avg=options.avg_preparedata;
end
if ~isempty(options) && ~isfield(options,'params_temporal_distance')
    options.params_temporal_distance=[1 1 1];
    params_temporal_distance=[1 1 1];
elseif ~isempty(options) && isfield(options,'params_temporal_distance')
    params_temporal_distance=options.params_temporal_distance; %% OVERWRITES EXTRA-INPUT!
end


if ~isempty(options) && isfield(options,'preparedata_centerofmass') && options.preparedata_centerofmass==1 % Adds center of mass as additional fish.
    temptray=NaN(size(trayectorias)+[0 1 0]);
    temptray(:,1:no_fish,:)=trayectorias;
    temptray(:,no_fish+1,:)=(mean(squeeze(trayectorias),2));
    trayectorias=temptray;
end
no_fish=size(trayectorias,2);
polar_axis=polar_axis';
polar_center=polar_center';

disp(options)

md=cell(no_fish,no_fish);
for fish1=focal
    
    for fish2=(fish1+1):no_fish
        if fish1~=fish2
            
            
            
            movementdata=cell(2,1);
            simplestats=cell(2,1);
            for k=1:2
                movementdata{k}=struct('focal_r',[],'focal_v',[],'focal_a',[],'focal_v_magnitude',[],'focal_a_magnitude',[],'neighbour_r',[],'neighbour_v',[],'neighbour_a',[],...
                    'distance_focal_neighbour',[],'distance_focal_nextframe_to_neighbour',[],'distance_neighbour_nextframe_to_focal',[],...
                    'distance_change_focal_nextframe_to_neighbour',[],'distance_change_acceleration_focal_nextframe_to_neighbour',[],...
                    'distance_change_neighbour_nextframe_to_focal',[],'distance_change_acceleration_neighbour_nextframe_to_focal',[],...
                    'temporal_distance_focal_neighbour',[],'focal_movement_direction',[],'focal_to_neighbour_direction',[],'focal_to_neighbour_vector',[],...
                    'neighbour_movement_direction',[],...
                    'angle_focal_movement_direction_and_y_axis',[],'angle_focal_movement_direction_and_neighbour_movement_direction',[],...
                    'angle_focal_movement_direction_and_neighbour_position',[],...
                    'angle_neighbour_movement_direction_and_focal_position',[],...
                    'angle_focal_movement_direction_t1_and_t2',[],...
                    'angle_focal_movement_direction_t1_and_t2_magnitude',[],...
                    'angle_neighbour_movement_direction_t1_and_t2',[],...
                    'angle_neighbour_movement_direction_t1_and_t2_magnitude',[],...
                    'angle_acceleration_focal_movement_direction_t1_and_t2',[],...
                    'angle_acceleration_neighbour_movement_direction_t1_and_t2',[],...
                    'angle_focal_movement_direction_t1_and_t2_integrated',[],...
                    'focal_to_neighbour_vector_rotated',[],...
                    'focal_to_neighbour_direction_rotated',[],...
                    'focal_v_rotated',[],'focal_a_rotated',[],...
                    'neighbour_v_rotated',[],...
                    'neighbour_v_transformed',[],...
                    'neighbour_movement_direction_rotated',[],'neighbour_a_rotated',[],...
                    'neighbour_a_transformed',[],...
                    'angle_focal_acceleration_and_neighbour_position',[],...
                    'focal_acceleration_direction',[],...
                    'angle_focal_acceleration_and_movement_direction',[],...
                    'angle_focal_movement_direction_and_focal_acceleration',[],...
                    'focal_position_direction',[],...
                    'angle_focal_position_and_y_axis',[],...
                    'focal_distance_to_origin',[],...
                    'focal_polar_radius',[],...
                    'focal_polar_angle',[],...
                    'neighbour_polar_radius',[],...
                    'neighbour_polar_angle',[],...
                    'focal_polar_radius_velocity',[],...
                    'focal_polar_angle_velocity',[],...
                    'focal_polar_angle_acceleration',[],...
                    'neighbour_polar_radius_velocity',[],...
                    'neighbour_polar_angle_velocity',[],...
                    'neighbour_polar_angle_acceleration',[],...
                    'focal_exclusion_distance',[],...
                    'focal_exclusion_angle',[],...
                    'angle_focal_to_polar_center_and_focal_neighbour_direction',[],...
                    'focal_x_to_polarcenter_direction',[],...
                    'neighbour_x_to_polarcenter_direction',[],...
                    'focal_x_to_polarcenter',[],...
                    'neighbour_x_to_polarcenter',[],...
                    'polar_angle_distance_focal_to_neighbour',[]);
            end
            
            
            
            tray_pair=trayectorias(:,[fish1 fish2],:);
            
            if length(focal)==no_fish
                nb_switch=[1 2; 2 1];
            else
                nb_switch=[1 2]';
            end
            
            for fish=nb_switch
                
                
                
                x=[tray_pair(:,fish(1),1) tray_pair(:,fish(2),1)]';
                y=[tray_pair(:,fish(1),2) tray_pair(:,fish(2),2)]';
                %     prob=prob(1:9000);
                disp(['[' mfilename '] ' 'strucdyn_preparedata is averaging over ' num2str(avg) ' frames'])
                
                
                
                
                %% For ringarena/polar coordinates with origin polar_center
                
                focal_r4polar=vertcat(x(1,:),y(1,:));
                nb_r4polar=vertcat(x(2,:),y(2,:));
                
                try
                    focal_x_to_polarcenter=bsxfun(@minus,polar_center,focal_r4polar);
                catch
                    keyboard
                end
                focal_x_to_polarcenter_direction=bsxfun(@rdivide,focal_x_to_polarcenter,sqrt(sum(focal_x_to_polarcenter.^2)));
                
                focal_polar_radius=sqrt(sum(focal_x_to_polarcenter.^2));
                focal_polar_angle=...
                    atan2(-focal_x_to_polarcenter_direction(1,:).*polar_axis(2,:)-polar_axis(1,:).*(-focal_x_to_polarcenter_direction(2,:)),...
                    (-focal_x_to_polarcenter_direction(1,:)).*polar_axis(1,:)+(-focal_x_to_polarcenter_direction(2,:)).*polar_axis(2,:)); % "-" because x_to_center points towards the center, now I want it point outwards.
                %             mod(atan2(t1(1).*t2(2)-t2(1).*t1(2),t1(1).*t2(1)+t1(2).*t2(2)),2*pi)
                neighbour_x_to_polarcenter=bsxfun(@minus,polar_center,nb_r4polar);
                neighbour_x_to_polarcenter_direction=bsxfun(@rdivide,neighbour_x_to_polarcenter,sqrt(sum(neighbour_x_to_polarcenter.^2)));
                
                neighbour_polar_radius=sqrt(sum(neighbour_x_to_polarcenter.^2));
                neighbour_polar_angle=...
                    atan2(-neighbour_x_to_polarcenter_direction(1,:).*polar_axis(2,:)-polar_axis(1,:).*(-neighbour_x_to_polarcenter_direction(2,:)),...
                    (-neighbour_x_to_polarcenter_direction(1,:)).*polar_axis(1,:)+(-neighbour_x_to_polarcenter_direction(2,:)).*polar_axis(2,:)); % "-" because x_to_center points towards the center, now I want it point outwards.
                
                
                
                if ~isempty(options) && isfield(options,'interpolate_polar_angle_and_radius') && options.interpolate_polar_angle_and_radius % First interpolate polar coordinates, then derive cartesian coordinates from it.
                    %                 dd_x=find(~isnan(focal_polar_angle));
                    %                 dd_y=focal_polar_angle(~isnan(focal_polar_angle));
                    %                 focal_polar_angle=interp1(dd_x,dd_y,1:length(focal_polar_angle));
                    
                    focal_polar_angle=interp_angles(focal_polar_angle);
                    
                    
                    dd_x=find(~isnan(focal_polar_radius));
                    dd_y=focal_polar_radius(~isnan(focal_polar_radius));
                    focal_polar_radius=interp1(dd_x,dd_y,1:length(focal_polar_radius));
                    
                    
                    
                    %                 dd_x=find(~isnan(neighbour_polar_angle));
                    %                 dd_y=neighbour_polar_angle(~isnan(neighbour_polar_angle));
                    %                 neighbour_polar_angle=interp1(dd_x,dd_y,1:length(neighbour_polar_angle));
                    
                    neighbour_polar_angle=interp_angles(neighbour_polar_angle);
                    
                    dd_x=find(~isnan(neighbour_polar_radius));
                    dd_y=neighbour_polar_radius(~isnan(neighbour_polar_radius));
                    neighbour_polar_radius=interp1(dd_x,dd_y,1:length(neighbour_polar_radius));
                    
                    disp('¡¡¡ATTENTION!!! No moving averages in polar coordinates. Function mediamovil4angles is faulty. Fix it!')
                    %                 focal_polar_angle=mediamovil4angles(focal_polar_angle,avg);
                    %                 neighbour_polar_angle=mediamovil4angles(neighbour_polar_angle,avg);
                    %                 focal_polar_radius=tsmovavg(focal_polar_radius,'s',avg);
                    %                 neighbour_polar_radius=tsmovavg(neighbour_polar_radius,'s',avg);
                    %
                    %                 focal_polar_angle(~filters{fish(1)})=NaN;
                    %                 focal_polar_radius(~filters{fish(1)})=NaN;
                    %                 neighbour_polar_angle(~filters{fish(1)})=NaN;
                    %                 neighbour_polar_radius(~filters{fish(1)})=NaN;
                    
                    
                    
                    
                    % Re-do carthesian coordinates from interpolated polars
                    %                 keyboard
                    focal_r= vertcat(focal_polar_radius.*cos(focal_polar_angle), -focal_polar_radius.*sin(focal_polar_angle));
                    focal_r(1,:)=focal_r(1,:)+polar_center(1); focal_r(2,:)=focal_r(2,:)+polar_center(2);
                    neighbour_r= vertcat(neighbour_polar_radius.*cos(neighbour_polar_angle), -neighbour_polar_radius.*sin(neighbour_polar_angle));
                    neighbour_r(1,:)=neighbour_r(1,:)+polar_center(1); neighbour_r(2,:)=neighbour_r(2,:)+polar_center(2);
                    focal_x_to_polarcenter=bsxfun(@minus,polar_center,focal_r);
                    %                 keyboard
                    xavg=vertcat(focal_r(1,:), neighbour_r(1,:));
                    yavg=vertcat(focal_r(2,:), neighbour_r(2,:));
                    
                    disp('Interpolation using polar coordinates.')
                else % if isempty(options) || options.interpolate_polar_angle_and_radius: Do not calculate focal_r etc. from the interpolated polar coordinates.
                    focal_polar_angle=mediamovil4angles(focal_polar_angle,avg);
                    neighbour_polar_angle=mediamovil4angles(neighbour_polar_angle,avg);
%                     focal_polar_radius=tsmovavg(focal_polar_radius,'s',avg);
                    focal_polar_radius= moving_average(avg,focal_polar_radius,[],2);
%                     neighbour_polar_radius=tsmovavg(neighbour_polar_radius,'s',avg);
                    neighbour_polar_radius= moving_average(avg,neighbour_polar_radius,[],2);
                    
                    %                  focal_v=vertcat(vx(1,:),vy(1,:));
                    
%                     xavg=tsmovavg(x,'s',avg); % Smoothing trajectories with moving window average, according to Katz11 (window size 10)
%                     yavg=tsmovavg(y,'s',avg);
                    xavg= moving_average(avg,x,[],2);
                    yavg= moving_average(avg,y,[],2);
                    focal_r=vertcat(xavg(1,:),yavg(1,:));
                    neighbour_r=vertcat(xavg(2,:),yavg(2,:));
                    
                end
                
                focal_to_neighbour_direction=bsxfun(@rdivide,neighbour_r-focal_r,sqrt(sum((neighbour_r-focal_r).^2)));
                
                
                focal_exclusion_distance=sqrt(focal_polar_radius.^2-exclusion_radius^2);
                focal_exclusion_angle=asin(bsxfun(@rdivide,exclusion_radius,focal_polar_radius));
                angle_focal_to_polar_center_and_focal_neighbour_direction=...
                    atan2(focal_x_to_polarcenter_direction(1,:).*focal_to_neighbour_direction(2,:)-focal_to_neighbour_direction(1,:).*focal_x_to_polarcenter_direction(2,:),...
                    focal_x_to_polarcenter_direction(1,:).*focal_to_neighbour_direction(1,:)+focal_x_to_polarcenter_direction(2,:).*focal_to_neighbour_direction(2,:));
                
                
                focal_polar_angle_velocity=[strucdyn_angle_difference(focal_polar_angle(1:end-1),focal_polar_angle(2:end)) NaN];
                focal_polar_angle_acceleration=[diff(focal_polar_angle_velocity) NaN];
                focal_polar_radius_velocity=[diff(focal_polar_radius) NaN];
                
                
                neighbour_polar_angle_velocity=[strucdyn_angle_difference(neighbour_polar_angle(1:end-1),neighbour_polar_angle(2:end)) NaN];
                neighbour_polar_angle_acceleration=[diff(neighbour_polar_angle_velocity) NaN];
                neighbour_polar_radius_velocity=[diff(neighbour_polar_radius) NaN];
                
                polar_angle_distance_focal_to_neighbour=strucdyn_angle_difference(focal_polar_angle,neighbour_polar_angle);
                
                %% Prepare filters.
                
                tray_pair_avg=NaN(size(trayectorias,1),2,2);%trayectorias(:,[fish1 fish2],:);
                tray_pair_avg(:,1,:)=[xavg(1,:)' yavg(1,:)'];
                tray_pair_avg(:,2,:)=[xavg(2,:)' yavg(2,:)'];
                if filters_on
                    if ~noprobs
                        filters=strucdyn_filters(tray_pair_avg,prob,avg,options);
                    else
                        filters=strucdyn_filters(tray_pair_avg,[],avg,options);
                    end
                else
                    filters=cell(2,1);
                    for k=1:2
                        filters{k}=true(size(tray_pair_avg(:,2,1)));
                        filters{k}(end-1:end)=false(1,2);
                    end
                end
                
                %%
                
                
                vx=diff(xavg,1,2);
                vy=diff(yavg,1,2);
                
                acceleration_x=diff(vx,1,2);
                acceleration_y=diff(vy,1,2);
                % Filtering the vecors:
                
                
                
                xavg(:,~filters{fish(1)})=NaN;
                yavg(:,~filters{fish(1)})=NaN;
                
                vx(:,~filters{fish(1)}(1:end-1))=NaN;
                vy(:,~filters{fish(1)}(1:end-1))=NaN;
                
                acceleration_x(:,~filters{fish(1)}(1:end-2))=NaN;
                acceleration_y(:,~filters{fish(1)}(1:end-2))=NaN;
                
                %             xavg=xavg(:,1:end-2);
                %             yavg=yavg(:,1:end-2);
                %             vx=vx(:,1:end-1);
                %             vy=vy(:,1:end-1);
                
                vxt=vx;
                vx=NaN(size(vx,1),size(xavg,2));
                vx(:,1:end-1)=vxt;
                vyt=vy;
                vy=NaN(size(vy,1),size(yavg,2));
                vy(:,1:end-1)=vyt;
                
                axt=acceleration_x;
                acceleration_x=NaN(size(vx,1),size(xavg,2));
                acceleration_x(:,1:end-2)=axt;
                ayt=acceleration_y;
                acceleration_y=NaN(size(vy,1),size(xavg,2));
                acceleration_y(:,1:end-2)=ayt;
                
                focal_r=vertcat(xavg(1,:),yavg(1,:));
                focal_v=vertcat(vx(1,:),vy(1,:));
                
                neighbour_r=vertcat(xavg(2,:),yavg(2,:));
                neighbour_v=vertcat(vx(2,:),vy(2,:));
             
                
                focal_polar_angle(~filters{fish(1)})=NaN;
                focal_polar_radius(~filters{fish(1)})=NaN;
                neighbour_polar_angle(~filters{fish(1)})=NaN;
                neighbour_polar_radius(~filters{fish(1)})=NaN;
                
                %% Transformation to the coordinate system of the focal fish
                
                %             distance_focal_neighbour=sqrt((xavg(1,:)-xavg(2,:)).^2+(yavg(1,:)-yavg(2,:)).^2);
                
                %             focal_r=vertcat(xavg(1,:),yavg(1,:));
                %             focal_v=vertcat(vx(1,:),vy(1,:));
                focal_v_magnitude=sqrt(focal_v(1,:).^2+focal_v(2,:).^2);
                focal_acceleration=vertcat(acceleration_x(1,:),acceleration_y(1,:));
                focal_a_magnitude=sqrt(focal_acceleration(1,:).^2+focal_acceleration(2,:).^2);
                %     keyboard
                %             neighbour_r=vertcat(xavg(2,:),yavg(2,:));
                %             neighbour_v=vertcat(vx(2,:),vy(2,:));
                neighbour_acceleration=vertcat(acceleration_x(2,:),acceleration_y(2,:));
                neighbour_v_magnitude=sqrt(neighbour_v(1,:).^2+neighbour_v(2,:).^2);
                neighbour_a_magnitude=sqrt(neighbour_acceleration(1,:).^2+neighbour_acceleration(2,:).^2);
                
                focal_position_direction=bsxfun(@rdivide,focal_r,sqrt(sum(focal_r.^2)));
                focal_distance_to_origin=sqrt(sum(focal_r.^2));
                focal_movement_direction=bsxfun(@rdivide,focal_v,sqrt(sum(focal_v.^2)));
                focal_to_neighbour_direction=bsxfun(@rdivide,neighbour_r-focal_r,sqrt(sum((neighbour_r-focal_r).^2)));
                focal_to_neighbour_vector=neighbour_r-focal_r;
                focal_acceleration_direction=bsxfun(@rdivide,focal_acceleration,sqrt(sum((focal_acceleration).^2)));
                
                neighbour_movement_direction=bsxfun(@rdivide,neighbour_v,sqrt(sum(neighbour_v.^2)));
                
                distance_focal_neighbour=sqrt((focal_r(1,:)-neighbour_r(1,:)).^2+(focal_r(2,:)-neighbour_r(2,:)).^2);
                
                distance_focal_nextframe_to_neighbour=[sqrt((focal_r(1,2:end)-neighbour_r(1,1:end-1)).^2+(focal_r(1,2:end)-neighbour_r(2,1:end-1)).^2) NaN];
                distance_neighbour_nextframe_to_focal=[sqrt((neighbour_r(1,2:end)-focal_r(1,1:end-1)).^2+(neighbour_r(1,2:end)-focal_r(2,1:end-1)).^2) NaN];
                
                distance_change_focal_nextframe_to_neighbour=distance_focal_nextframe_to_neighbour-distance_focal_neighbour;
                distance_change_acceleration_focal_nextframe_to_neighbour=[diff(distance_change_focal_nextframe_to_neighbour) NaN];
                
                distance_change_neighbour_nextframe_to_focal=distance_neighbour_nextframe_to_focal-distance_focal_neighbour;
                distance_change_acceleration_neighbour_nextframe_to_focal=[diff(distance_change_neighbour_nextframe_to_focal) NaN];
                
                % Angle between the direction of movement of the focal fish and the y-axis
                %
                %             If you'd like to test the sense:
                %             t1=[1 0];
                %             t2=[0 1];
                %             mod(atan2(t1(1).*t2(2)-t2(1).*t1(2),t1(1).*t2(1)+t1(2).*t2(2)),2*pi)
                %             /pi*180
                %             >0! Therefore: counterclockwise gives positive sign.
                %
                angle_focal_movement_direction_and_y_axis=mod(atan2(focal_movement_direction(1,:)*1-0*focal_movement_direction(2,:),focal_movement_direction(1,:)*0+focal_movement_direction(2,:)*1),2*pi);
                angle_focal_movement_direction_and_neighbour_movement_direction=mod(atan2(neighbour_movement_direction(1,:).*focal_movement_direction(2,:)-focal_movement_direction(1,:).*neighbour_movement_direction(2,:),focal_movement_direction(1,:).*neighbour_movement_direction(1,:)+focal_movement_direction(2,:).*neighbour_movement_direction(2,:)),2*pi);
%                 angle_focal_movement_direction_and_neighbour_position=mod(atan2(focal_movement_direction(1,:).*focal_to_neighbour_direction(2,:)-focal_to_neighbour_direction(1,:).*focal_movement_direction(2,:),focal_movement_direction(1,:).*focal_to_neighbour_direction(1,:)+focal_movement_direction(2,:).*focal_to_neighbour_direction(2,:)),2*pi);
                angle_focal_movement_direction_and_neighbour_position=mod(atan2(focal_to_neighbour_direction(1,:).*focal_movement_direction(2,:)-focal_movement_direction(1,:).*focal_to_neighbour_direction(2,:),focal_movement_direction(1,:).*focal_to_neighbour_direction(1,:)+focal_movement_direction(2,:).*focal_to_neighbour_direction(2,:)),2*pi);
                
angle_neighbour_movement_direction_and_focal_position=mod(atan2(neighbour_movement_direction(1,:).*(-focal_to_neighbour_direction(2,:))-(-focal_to_neighbour_direction(1,:)).*neighbour_movement_direction(2,:),neighbour_movement_direction(1,:).*(-focal_to_neighbour_direction(1,:))+neighbour_movement_direction(2,:).*(-focal_to_neighbour_direction(2,:))),2*pi);
                angle_focal_acceleration_and_neighbour_position=atan2(focal_acceleration_direction(1,:).*(focal_to_neighbour_direction(2,:))-(focal_to_neighbour_direction(1,:)).*focal_acceleration_direction(2,:),focal_acceleration_direction(1,:).*(focal_to_neighbour_direction(1,:))+focal_acceleration_direction(2,:).*(focal_to_neighbour_direction(2,:)));
                %             angle_focal_acceleration_and_movement_direction=mod(atan2(focal_acceleration_direction(1,:).*(focal_movement_direction(2,:))-(focal_movement_direction(1,:)).*focal_acceleration_direction(2,:),focal_acceleration_direction(1,:).*(focal_movement_direction(1,:))+focal_acceleration_direction(2,:).*(focal_movement_direction(2,:))),2*pi);
                angle_focal_acceleration_and_movement_direction=atan2(focal_acceleration_direction(1,:).*(focal_movement_direction(2,:))-(focal_movement_direction(1,:)).*focal_acceleration_direction(2,:),focal_acceleration_direction(1,:).*(focal_movement_direction(1,:))+focal_acceleration_direction(2,:).*(focal_movement_direction(2,:)));
                
                angle_focal_movement_direction_and_focal_acceleration=mod(atan2(focal_acceleration_direction(2,:).*(focal_movement_direction(1,:))-(focal_movement_direction(2,:)).*focal_acceleration_direction(1,:),focal_acceleration_direction(1,:).*(focal_movement_direction(1,:))+focal_acceleration_direction(2,:).*(focal_movement_direction(2,:))),2*pi);
                
                angle_focal_position_and_y_axis=mod(atan2(focal_position_direction(1,:)*1-0*focal_position_direction(2,:),focal_position_direction(1,:)*0+focal_position_direction(2,:)*1),2*pi);
                %             angle_focal_acceleration_and_neighbour_position=mod(atan2(focal_acceleration_direction(1,:).*(focal_movement_direction(2,:))-(focal_movement_direction(1,:)).*focal_acceleration_direction(2,:),focal_acceleration_direction(1,:).*(focal_movement_direction(1,:))+focal_acceleration_direction(2,:).*(focal_movement_direction(2,:))),2*pi);
                angle_focal_movement_direction_t1_and_t2=[atan2(focal_movement_direction(1,2:end).*focal_movement_direction(2,1:end-1)-focal_movement_direction(1,1:end-1).*focal_movement_direction(2,2:end),focal_movement_direction(1,2:end).*focal_movement_direction(1,1:end-1)+focal_movement_direction(2,2:end).*focal_movement_direction(2,1:end-1)) NaN];
               angle_focal_movement_direction_t1_and_t2_magnitude=abs(angle_focal_movement_direction_t1_and_t2);
                angle_acceleration_focal_movement_direction_t1_and_t2=[diff(angle_focal_movement_direction_t1_and_t2,1,2) NaN];
                
                angle_neighbour_movement_direction_t1_and_t2=[atan2(neighbour_movement_direction(1,2:end).*neighbour_movement_direction(2,1:end-1)-neighbour_movement_direction(1,1:end-1).*neighbour_movement_direction(2,2:end),neighbour_movement_direction(1,2:end).*neighbour_movement_direction(1,1:end-1)+neighbour_movement_direction(2,2:end).*neighbour_movement_direction(2,1:end-1)) NaN];
                angle_neighbour_movement_direction_t1_and_t2_magnitude=abs(angle_neighbour_movement_direction_t1_and_t2);
                angle_acceleration_neighbour_movement_direction_t1_and_t2=[diff(angle_neighbour_movement_direction_t1_and_t2,1,2) NaN];
                
                
                %             angle_focal_movement_direction_t1_and_t2=[angle_focal_movement_direction_t1_and_t2 NaN];
                
                %----------------Integrating over
                %angle_focal_movement_direction_t1_and_t2 until direction
                %changes---------------------------------------------------
                
                da=angle_focal_movement_direction_t1_and_t2(end:-1:1);
                da(isnan(angle_focal_movement_direction_t1_and_t2(end:-1:1)))=0;
                da(da<0)=0;
                
                cx1=cumsum(da);
                
                cx2=zeros(size(cx1));
                i=strfind([0 da>0.'],[0 1]);
                cx2(i)=diff([0; cx1(i)'-da(i)']);
                x3=cx1-cumsum(cx2);
                % x3(~x) = NaN;
                x3(da<=0)=0;
                
                da=angle_focal_movement_direction_t1_and_t2(end:-1:1);
                da(isnan(angle_focal_movement_direction_t1_and_t2(end:-1:1)))=0;
                da(da>0)=0;
                cx1=cumsum(da);
                
                cx3=zeros(size(cx1));
                i=strfind([0 da<0.'],[0 1]);
                cx3(i)=diff([0; cx1(i)'-da(i)']);
                x4=cx1-cumsum(cx3);
                % x4(~x(:,1)) = NaN;
                
                x4(da>=0)=0;
                angle_focal_movement_direction_t1_and_t2_integrated=x3;
                angle_focal_movement_direction_t1_and_t2_integrated(x3==0)=x4(x3==0);
                angle_focal_movement_direction_t1_and_t2_integrated=angle_focal_movement_direction_t1_and_t2_integrated(end:-1:1);
                
                temporal_distance_focal_neighbour=distance_focal_neighbour/params_temporal_distance(1) + ...
                    min(angle_focal_movement_direction_and_neighbour_position,2*pi-angle_focal_movement_direction_and_neighbour_position)*180/pi/params_temporal_distance(2)+ ...
                    params_temporal_distance(3);
%                 keyboard
                %----------------------------------------------------------
                
                
                %Preallocating
              
                neighbour_v_transformed=NaN(2,length(angle_focal_movement_direction_and_y_axis));
               %             neighbour_acceleration_transformed=NaN(2,length(angle_focal_movement_direction_and_y_axis));
                
                
                
                a=-angle_focal_movement_direction_and_y_axis';%horzcat(angle_focal_movement_direction_and_y_axis, NaN)';
                turnmatr=NaN(2,2,length(a));
                turnmatr(1,1,:)=cos(a);
                turnmatr(1,2,:)=-sin(a);
                turnmatr(2,1,:)=sin(a);
                turnmatr(2,2,:)=cos(a);
                
                focal_to_neighbour_vector_rotated=[sum(focal_to_neighbour_vector.*squeeze(turnmatr(:,1,:)))' sum(focal_to_neighbour_vector.*squeeze(turnmatr(:,2,:)))']'; %bsxfun(@times,focal_to_neighbour_vector,squeeze(turnmatr(:,2,:)))];
                focal_to_neighbour_direction_rotated=[sum(focal_to_neighbour_direction.*squeeze(turnmatr(:,1,:)))' sum(focal_to_neighbour_direction.*squeeze(turnmatr(:,2,:)))']'; %bsxfun(@times,focal_to_neighbour_vector,squeeze(turnmatr(:,2,:)))];
                
                focal_v_rotated=[sum(focal_v.*squeeze(turnmatr(:,1,:)))' sum(focal_v.*squeeze(turnmatr(:,2,:)))']'; %bsxfun(@times,focal_to_neighbour_vector,squeeze(turnmatr(:,2,:)))];
                neighbour_v_rotated=[sum(neighbour_v.*squeeze(turnmatr(:,1,:)))' sum(neighbour_v.*squeeze(turnmatr(:,2,:)))']'; %bsxfun(@times,focal_to_neighbour_vector,squeeze(turnmatr(:,2,:)))];
                neighbour_movement_direction_rotated=[sum(neighbour_movement_direction.*squeeze(turnmatr(:,1,:)))' sum(neighbour_movement_direction.*squeeze(turnmatr(:,2,:)))']'; %bsxfun(@times,focal_to_neighbour_vector,squeeze(turnmatr(:,2,:)))];

                focal_acceleration_rotated=[sum(focal_acceleration.*squeeze(turnmatr(:,1,:)))' sum(focal_acceleration.*squeeze(turnmatr(:,2,:)))']'; %bsxfun(@times,focal_to_neighbour_vector,squeeze(turnmatr(:,2,:)))];
                neighbour_acceleration_rotated=[sum(neighbour_acceleration.*squeeze(turnmatr(:,1,:)))' sum(neighbour_acceleration.*squeeze(turnmatr(:,2,:)))']'; %bsxfun(@times,focal_to_neighbour_vector,squeeze(turnmatr(:,2,:)))];

                             
                            
                
                
              
                neighbour_v_transformed(1,:)=neighbour_v_rotated(1,:);
                neighbour_v_transformed(2,:)=neighbour_v_rotated(2,:)-focal_v_rotated(2,:);
                neighbour_acceleration_transformed=horzcat(diff(neighbour_v_transformed,1,2), NaN(2,1));
                
              
                
                %             angle_focal_to_neighbour_direction_t1_and_t2=angle
                
                clear alpha
                %             keyboard
                
                
                %             %% Transformation to polar coordinates (polar_radius,polar_angle)
                %
                %
                %             x=[tray_pair(:,fish(1),1) tray_pair(:,fish(2),1)]';
                %             y=[tray_pair(:,fish(1),2) tray_pair(:,fish(2),2)]';
                %
                %             xavg=tsmovavg(x,'s',avg); % Smoothing trajectories with moving window average, according to Katz11 (window size 10)
                %             yavg=tsmovavg(y,'s',avg);
                %
                % %             xavg(:,~filters{fish(1)})=NaN;
                % %             yavg(:,~filters{fish(1)})=NaN;
                
                
                
                %%
                simplestats{fish(1)}.focal_mean_velocity=nanmean(focal_v_magnitude);
                simplestats{fish(1)}.focal_mean_acceleration=nanmean(focal_a_magnitude);
                simplestats{fish(1)}.neighbour_mean_velocity=nanmean(neighbour_v_magnitude);
                simplestats{fish(1)}.neighbour_mean_acceleration=nanmean(neighbour_a_magnitude);
                simplestats{fish(1)}.focal_std_velocity=nanstd(focal_v_magnitude);
                simplestats{fish(1)}.focal_std_acceleration=nanstd(focal_a_magnitude);
                simplestats{fish(1)}.neighbour_std_velocity=nanstd(neighbour_v_magnitude);
                simplestats{fish(1)}.neighbour_std_acceleration=nanstd(neighbour_a_magnitude);
                
                movementdata{fish(1)}.angle_focal_movement_direction_t1_and_t2=angle_focal_movement_direction_t1_and_t2;
                movementdata{fish(1)}.angle_focal_movement_direction_t1_and_t2_magnitude=angle_focal_movement_direction_t1_and_t2_magnitude;

                movementdata{fish(1)}.angle_neighbour_movement_direction_t1_and_t2=angle_neighbour_movement_direction_t1_and_t2;
                movementdata{fish(1)}.angle_neighbour_movement_direction_t1_and_t2_magnitude=angle_neighbour_movement_direction_t1_and_t2_magnitude;
                movementdata{fish(1)}.angle_acceleration_focal_movement_direction_t1_and_t2=angle_acceleration_focal_movement_direction_t1_and_t2;
                movementdata{fish(1)}.angle_acceleration_neighbour_movement_direction_t1_and_t2=angle_acceleration_neighbour_movement_direction_t1_and_t2;
                movementdata{fish(1)}.angle_focal_movement_direction_t1_and_t2_integrated=angle_focal_movement_direction_t1_and_t2_integrated;
                movementdata{fish(1)}.neighbour_movement_direction_rotated=neighbour_movement_direction_rotated;
                movementdata{fish(1)}.focal_r=focal_r;
                movementdata{fish(1)}.focal_v=focal_v;
                movementdata{fish(1)}.focal_v_magnitude=focal_v_magnitude;
                movementdata{fish(1)}.focal_a_magnitude=focal_a_magnitude;
                movementdata{fish(1)}.focal_a=focal_acceleration;
                movementdata{fish(1)}.neighbour_r=neighbour_r;
                movementdata{fish(1)}.neighbour_v=neighbour_v;
                movementdata{fish(1)}.neighbour_a=neighbour_acceleration;
                movementdata{fish(1)}.neighbour_a_magnitude=neighbour_a_magnitude;
                movementdata{fish(1)}.neighbour_v_magnitude=neighbour_v_magnitude;
                movementdata{fish(1)}.neighbour_movement_direction=neighbour_movement_direction;
                movementdata{fish(1)}.neighbour_movement_direction_rotated=neighbour_movement_direction_rotated;
                movementdata{fish(1)}.angle_focal_movement_direction_and_y_axis=angle_focal_movement_direction_and_y_axis;
                movementdata{fish(1)}.angle_focal_movement_direction_and_neighbour_movement_direction=angle_focal_movement_direction_and_neighbour_movement_direction;
                movementdata{fish(1)}.angle_focal_movement_direction_and_neighbour_position=angle_focal_movement_direction_and_neighbour_position;
                movementdata{fish(1)}.angle_neighbour_movement_direction_and_focal_position=angle_neighbour_movement_direction_and_focal_position;
                movementdata{fish(1)}.focal_to_neighbour_vector_rotated=focal_to_neighbour_vector_rotated;
                movementdata{fish(1)}.focal_to_neighbour_direction_rotated=focal_to_neighbour_direction_rotated;
                movementdata{fish(1)}.focal_v_rotated=focal_v_rotated;
                movementdata{fish(1)}.neighbour_v_rotated=neighbour_v_rotated;
                movementdata{fish(1)}.neighbour_v_transformed=neighbour_v_transformed;
                movementdata{fish(1)}.neighbour_a_transformed=neighbour_acceleration_transformed;
                movementdata{fish(1)}.focal_a_rotated=focal_acceleration_rotated;
                movementdata{fish(1)}.neighbour_a_rotated=neighbour_acceleration_rotated;
                movementdata{fish(1)}.distance_focal_neighbour=distance_focal_neighbour;
                movementdata{fish(1)}.distance_focal_nextframe_to_neighbour=distance_focal_nextframe_to_neighbour;
                movementdata{fish(1)}.distance_neighbour_nextframe_to_focal=distance_neighbour_nextframe_to_focal;
                movementdata{fish(1)}.distance_change_focal_nextframe_to_neighbour=distance_change_focal_nextframe_to_neighbour;
                movementdata{fish(1)}.distance_change_acceleration_focal_nextframe_to_neighbour=distance_change_acceleration_focal_nextframe_to_neighbour;
                movementdata{fish(1)}.distance_change_neighbour_nextframe_to_focal=distance_change_neighbour_nextframe_to_focal;
                movementdata{fish(1)}.distance_change_acceleration_neighbour_nextframe_to_focal=distance_change_acceleration_neighbour_nextframe_to_focal;
                movementdata{fish(1)}.temporal_distance_focal_neighbour=temporal_distance_focal_neighbour;
                movementdata{fish(1)}.focal_movement_direction=focal_movement_direction;
                movementdata{fish(1)}.focal_to_neighbour_direction=focal_to_neighbour_direction;
                movementdata{fish(1)}.focal_to_neighbour_vector=focal_to_neighbour_vector;
                movementdata{fish(1)}.angle_focal_acceleration_and_neighbour_position=angle_focal_acceleration_and_neighbour_position;
                movementdata{fish(1)}.focal_acceleration_direction=focal_acceleration_direction;
                movementdata{fish(1)}.angle_focal_acceleration_and_movement_direction=angle_focal_acceleration_and_movement_direction;
                movementdata{fish(1)}.angle_focal_movement_direction_and_focal_acceleration=angle_focal_movement_direction_and_focal_acceleration;
                movementdata{fish(1)}.angle_focal_position_and_y_axis=angle_focal_position_and_y_axis;
                movementdata{fish(1)}.focal_position_direction=focal_position_direction;
                movementdata{fish(1)}.focal_distance_to_origin=focal_distance_to_origin;
                movementdata{fish(1)}.focal_polar_radius=focal_polar_radius;
                movementdata{fish(1)}.focal_polar_angle=focal_polar_angle;
                movementdata{fish(1)}.focal_polar_angle_velocity=focal_polar_angle_velocity;
                movementdata{fish(1)}.focal_polar_angle_acceleration=focal_polar_angle_acceleration;
                movementdata{fish(1)}.focal_polar_radius_velocity=focal_polar_radius_velocity;
                movementdata{fish(1)}.neighbour_polar_radius=neighbour_polar_radius;
                movementdata{fish(1)}.focal_x_to_polarcenter_direction=focal_x_to_polarcenter_direction;
                movementdata{fish(1)}.neighbour_x_to_polarcenter_direction=neighbour_x_to_polarcenter_direction;
                movementdata{fish(1)}.focal_x_to_polarcenter=focal_x_to_polarcenter;
                movementdata{fish(1)}.neighbour_x_to_polarcenter=neighbour_x_to_polarcenter;
                movementdata{fish(1)}.neighbour_polar_angle=neighbour_polar_angle;
                movementdata{fish(1)}.neighbour_polar_angle_velocity=neighbour_polar_angle_velocity;
                movementdata{fish(1)}.neighbour_polar_angle_acceleration=neighbour_polar_angle_acceleration;
                movementdata{fish(1)}.neighbour_polar_radius_velocity=neighbour_polar_radius_velocity;
                movementdata{fish(1)}.focal_exclusion_distance=focal_exclusion_distance;
                movementdata{fish(1)}.focal_exclusion_angle=focal_exclusion_angle;
                movementdata{fish(1)}.angle_focal_to_polar_center_and_focal_neighbour_direction=angle_focal_to_polar_center_and_focal_neighbour_direction;
                movementdata{fish(1)}.polar_angle_distance_focal_to_neighbour=polar_angle_distance_focal_to_neighbour;
                %     distance_rotated{fish(1)}=sqrt((focal_to_neighbour_vector_rotated(1,:)).^2+(focal_to_neighbour_vector_rotated(2,:)).^2);
            end
            if no_fish>2 || isfield(options,'movementdata_newformat') && options.movementdata_newformat %&& length(focal)==no_fish
                md{fish1,fish2}=movementdata{1};
                md{fish2,fish1}=movementdata{2};
                if ~isempty(options)
                    md{1,1}=options;
                end
                
            end
        end
    end
end
if no_fish>2 && length(focal)~=no_fish || isfield(options,'movementdata_newformat') && options.movementdata_newformat && length(focal)~=no_fish
    
    md2=cell(length(focal),no_fish);
    for ffish=focal
        for nfish=1:no_fish
            if nfish~=ffish
                md2{ffish,nfish}=md{ffish,nfish};
            end
        end
        
    end
    md=md2;
    if ~isempty(options)
        md{min(focal),min(focal)}=options;
    end
elseif no_fish==2 && ~(isfield(options,'movementdata_newformat') && options.movementdata_newformat)
    md=movementdata;
end

disp(['[' mfilename '] ' 'Movementdata is ready! (It took ' num2str(toc) ' seconds to prepare it)'])
