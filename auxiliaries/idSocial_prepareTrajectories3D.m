function [movementdata, options]=idSocial_prepareTrajectories3D(trayectorias,probtrayectorias,options,bodylength,disp_options)
% Calculate dynamic data from trajectory.
% !! Trajectory from #Frames-by-#Individuals-by-#Dimensions to #Frames-by-#Individuals-by-#Individuals-by-#Dimensions
% def_options.smooth_method='moving';
% def_options.smooth_degree=3;
% def_options.interpolation_mode='spline';
% def_options.interp_maxlen=inf;
% def_options.interpolate_polar_angle_and_radius=false;
% def_options.polar_center=[0 0 0];
% def_options.polar_axis=[xax; yax; zax];
% def_options.exclusion_radius=0;
% def_options.params_temporal_distance=[1 1 1];
% def_options.preparedata_centerofmass=false;
% def_options.arena_limits=[];

% Output: structure MOVEMENTDATA with fields:
if nargin <5 || isempty(disp_options)
    disp_options = true;
end
% if nargin<6 || isempty(structOut)
%     structOut = false;
% end

movementdata = [];
xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation
no_dim=3;
%%
act_method=mfilename;
act_method=act_method(10:end);

def_options.act_method=act_method;%(10:end);


% Def. options for dynamic data
def_options.smooth_method='moving';
def_options.smooth_degree=3;
def_options.median_filter=false;
def_options.median_filter_order=1;
def_options.interpolation_mode='spline';
def_options.interpolate_trajectories=false;
def_options.interp_maxlen=inf;
def_options.interpolate_polar_angle_and_radius=false;
def_options.preparedata_centerofmass=false;
def_options.smooth_max_deviation = 5;
def_options.smooth_spline_degree = 2;
def_options.smooth_adaptive_noise = 10;
def_options.filter_minProbIdentityAssignment=-inf;
% def_options.polar_center=[0 0 0];
% def_options.polar_axis=[xax; yax; zax];
% def_options.exclusion_radius=0;
% def_options.params_temporal_distance=[1 1 1];
% def_options.arena_limits=[];
% def_options.rotate_z=true;
% def_options.start_min=1;
% def_options.end_min=inf;

if nargin < 1
    movementdata = def_options;
    return;
end

if nargin<2 || isempty(options)
    options=[];
end
if nargin< 4 || isempty(bodylength)
    bodylength=inf;
end

defoptnames=fieldnames(def_options);
no_defoptions=size(defoptnames,1);
for optfield=1:no_defoptions
    if ~isfield(options,defoptnames{optfield}) %% || isempty(options.(defoptnames{optfield}))
        options.(defoptnames{optfield})=def_options.(defoptnames{optfield});
    end
end
if disp_options
    disp(options)
end
tic;


focalReconstructionFlag = isfield(options,'focalReconstruction_minProbIdentity4VelCalcs') && ~isempty(options.focalReconstruction_minProbIdentity4VelCalcs) && ...
        isfield(options,'focalReconstruction_minProbIdentityFocal') && ~isempty(options.focalReconstruction_minProbIdentityFocal) && ...
        isfield(options,'focalReconstruction_minProbIdentityNeighbor') && ~isempty(options.focalReconstruction_minProbIdentityNeighbor);
%% Options for dynamics:
smooth_spline_degree=options.smooth_spline_degree;
smooth_degree=options.smooth_degree;
smooth_method=options.smooth_method;
median_filter=options.median_filter;
median_filter_order=options.median_filter_order;
interp_mode=options.interpolation_mode;
interpolate_trajectories=options.interpolate_trajectories;
interp_maxlen=options.interp_maxlen;
% polar_center=options.polar_center;
interpolate_polar_angle_and_radius=options.interpolate_polar_angle_and_radius;
% polar_axis=options.polar_axis;
% exclusion_radius=options.exclusion_radius;
% params_temporal_distance=options.params_temporal_distance;
preparedata_centerofmass=options.preparedata_centerofmass;
% arena_limits=options.arena_limits;
% rotate_z=options.rotate_z;
% idces=[options.start_min min(options.end_min,size(trayectorias,1))];
idces=[1 size(trayectorias,1)];
max_deviation = options.smooth_max_deviation;
noise_mean = options.smooth_adaptive_noise;
%%
trayectorias=trayectorias(idces(1):idces(2),:,:);
no_frames=size(trayectorias,1);
no_fish=size(trayectorias,2);
no_dim_orig=size(trayectorias,3);

% [tr,vel,acc,no_frames,no_fish,no_dim] = ...
%     idSocial_auxiliaries_formatInputTrajectory(trayectorias);
% no_vertices=8;


% if no_dim_orig==2
%     trtemp=zeros(no_frames,no_fish,3);
%     trtemp(:,:,1:no_dim_orig)=trayectorias;
%     trayectorias=trtemp;
% end
clear trtemp

%% For center of mass calculations: Add c.o.m. coordinate as individual no_fish+1.
tr=trayectorias;
if preparedata_centerofmass
    temptray=NaN(size(tr2)+[0 1 0]);
    temptray(:,1:no_fish,:)=tr2;
    temptray(:,no_fish+1,:)=(mean(squeeze(tr2),2));
    tr=temptray;
    no_fish=no_fish+1;
end
clear trayectorias

% Add neighbour dimension
% tr=zeros(no_frames,no_fish,no_fish,no_dim);
% tr(:,:,:,1:no_dim)=reshape(repmat(tr2,[1,no_fish]),no_frames,no_fish,no_fish,no_dim);

%% Filter
% Apply filter of assignment
% probabilites
if ~isempty(options.filter_minProbIdentityAssignment) && ~focalReconstructionFlag
    tr(...
        repmat(probtrayectorias,...
        [1 1 size(tr,3)])<...
        options.filter_minProbIdentityAssignment) = NaN;
end

if median_filter
    for ff=1:no_fish
        for d=1:no_dim
            tr(:,ff,d)=medfilt1(tr(:,ff,d),median_filter_order);
        end
    end
end

if no_fish==1 && ndims(tr)<3
    tr=reshape(tr,[no_frames no_fish no_dim]);   
end
if ~isempty(smooth_method) && ~isempty(smooth_degree) && isnumeric(smooth_degree) && ...
        ~strcmpi(smooth_method,'none')
    [tr_temp, tol] = idSocial_auxiliaries_smoothTrajectories(tr,smooth_method,smooth_degree,max_deviation,noise_mean);
end


%% Generate structure for output
movementdata.options=options;
if isnumeric(tr_temp)
    vel=NaN(size(tr_temp));
    acc=NaN(size(tr_temp));
    vel(1:no_frames-1,:,:,:)=diff(tr_temp,1,1);
    acc(1:no_frames-2,:,:,:)=diff(tr_temp,2,1);
    movementdata.Tr=tr_temp(:,:,1:no_dim_orig);
    movementdata.Vel = vel;
    movementdata.Acc = acc;
else
     movementdata.Tr=tr_temp;
end

if strcmpi(smooth_method,'adaptive_moving') || strcmpi(smooth_method,'adaptive_moving_acc')
    movementdata.moving_average_width = mavgwidth;
end
if disp_options
    disp(['[' mfilename '] ' 'Trajectories are ready! (It took ' num2str(toc) ' seconds to prepare it)'])
end
