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
% no_vertices=8;


if no_dim_orig==2
    trtemp=zeros(no_frames,no_fish,3);
    trtemp(:,:,1:no_dim_orig)=trayectorias;
    trayectorias=trtemp;
end
clear trtemp

% if isempty(arena_limits)
%     minx=floor(min(min(min(trayectorias(:,:,1))))*1);
%     miny=floor(min(min(min(trayectorias(:,:,2))))*1);
%     minz=floor(min(min(min(trayectorias(:,:,3))))*1);
%     maxx=floor(max(max(max(trayectorias(:,:,1))))*1);
%     maxy=floor(max(max(max(trayectorias(:,:,2))))*1);
%     maxz=floor(max(max(max(trayectorias(:,:,3))))*1);
%     arena_limits=[minx miny minz; ...
%         minx maxy minz; ...
%         maxx maxy minz; ...
%         maxx miny minz; ...
%         minx miny maxz; ...
%         minx maxy maxz; ...
%         maxx maxy maxz; ...
%         maxx miny maxz];
%     options.arena_limits=arena_limits;
%
%
% end

%% For center of mass calculations: Add c.o.m. coordinate as individual no_fish+1.
tr2=trayectorias;
if preparedata_centerofmass
    temptray=NaN(size(tr2)+[0 1 0]);
    temptray(:,1:no_fish,:)=tr2;
    temptray(:,no_fish+1,:)=(mean(squeeze(tr2),2));
    tr2=temptray;
    no_fish=no_fish+1;
end
clear trayectorias
%% Trajectories with separate probtrajectories limits for focal and neighbors
% (Useful for example in situations where only neighbor positions are
% necessary, but not neighbor identities and velocity/acceleration)
% if isfield(options,'focalReconstruction_minProbIdentity4VelCalcs') && ~isempty(options.focalReconstruction_minProbIdentity4VelCalcs) && ...
%         isfield(options,'focalReconstruction_minProbIdentityFocal') && ~isempty(options.focalReconstruction_minProbIdentityFocal) && ...
%         isfield(options,'focalReconstruction_minProbIdentityNeighbor') && ~isempty(options.focalReconstruction_minProbIdentityNeighbor)
%     [tr, good_frames]= ...
%         idSocial_auxiliaries_idTrackerNoIdentities(tr2,probtrayectorias,options.focalReconstruction_minProbIdentity4VelCalcs, ...
%         options.focalReconstruction_minProbIdentityFocal, ...
%         options.focalReconstruction_minProbIdentityNeighbor);
%     duration_total = sum(good_frames);
%     duration_trial = no_frames;
%     if duration_total<duration_trial
%         warning([mfilename ': % of frames for all focal after reconstruction: ' num2str(duration_total)])
%     else
%         disp([mfilename ': total number of frames combining all focals after reconstruction: ' num2str(duration_total)])
%     end
% else
% Add neighbour dimension
tr=zeros(no_frames,no_fish,no_fish,no_dim);
tr(:,:,:,1:no_dim)=reshape(repmat(tr2,[1,no_fish]),no_frames,no_fish,no_fish,no_dim);
%% Filter
% Apply filter of assignment
% probabilites
if ~isempty(options.filter_minProbIdentityAssignment) && ~focalReconstructionFlag
    tr(...
        repmat(probtrayectorias,...
        [1 1 size(tr,2) size(tr,4)])<...
        options.filter_minProbIdentityAssignment) = NaN;
end
% end





%% Interpolation and smoothing either from polar or cartesian coordinates
if interpolate_polar_angle_and_radius % Transform to polar/cylindrical coordinates, then interpolate (and smooth, but smoothing is disabled right now...disculpen las molestias) radius and angle.
    
    foc_to_polarcenter=cat(4,polar_center(1)-tr(:,:,:,1), polar_center(2)-tr(:,:,:,2), polar_center(3)-tr(:,:,:,3));
    foc_to_polarcenter_dir=bsxfun(@rdivide,foc_to_polarcenter,sqrt(sum(foc_to_polarcenter.^2,4)));
    foc_polar_radius=sqrt(sum(foc_to_polarcenter.^2,4));
    
    foc_polar_angle_xz=...
        atan2(polar_axis(1,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(1,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))),...
        (polar_axis(1,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(1,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))))...
        );
    foc_polar_angle_z=...
        atan2(-foc_to_polarcenter_dir(:,:,:,1)*zax(2)-zax(1)*(-foc_to_polarcenter_dir(:,:,:,2)),(-foc_to_polarcenter_dir(:,:,:,1))*zax(1)+(-foc_to_polarcenter_dir(:,:,:,2))*zax(2));
    
    foc_polar_angle_xy=NaN(size(foc_polar_angle_z));
    foc_polar_angle_xy(foc_polar_angle_z>=0)=...
        pi/2-foc_polar_angle_z(foc_polar_angle_z>=0);
    foc_polar_angle_xy(foc_polar_angle_z<0)=...
        -pi/2-foc_polar_angle_z(foc_polar_angle_z<0);
    % The following would calculate the angle in "real" 3D
    % foc_polar_angle_xz=...
    %     atan2(sqrt(...
    %         (polar_axis(2).*(-squeeze(foc_to_polarcenter_dir(:,ff,3))) - polar_axis(3).*(-squeeze(foc_to_polarcenter_dir(:,ff,2)))).^2+...
    %         (polar_axis(3).*(-squeeze(foc_to_polarcenter_dir(:,ff,1))) - polar_axis(1).*(-squeeze(foc_to_polarcenter_dir(:,ff,3)))).^2+...
    %         (polar_axis(1).*(-squeeze(foc_to_polarcenter_dir(:,ff,2))) - polar_axis(2).*(-squeeze(foc_to_polarcenter_dir(:,ff,1)))).^2)',...
    %         polar_axis*(-squeeze(foc_to_polarcenter_dir(:,ff,:))')...
    %     );
    
    if ~isempty(interp_mode)
        disp([act_method ': Interpolation using polar coordinates.'])
        for ff=1:no_fish
            foc_polar_angle_xz(:,ff,:)=repmat(interp_angles(foc_polar_angle_xz(:,ff,ff)),[1,1,no_fish]);
            foc_polar_angle_xy(:,ff,:)=repmat(interp_angles(foc_polar_angle_xy(:,ff,ff)),[1,1,no_fish]);
        end
        
        foc_polar_radius=idSocial_interpolateTrajectories(foc_polar_radius,interp_mode,interp_maxlen);
    end
    disp('¡¡¡ATTENTION!!! No moving averages in polar coordinates. Function mediamovil4angles is faulty. Fix it!')
    
    tr= cat(4,foc_polar_radius.*cos(foc_polar_angle_xz), -foc_polar_radius.*sin(foc_polar_angle_xz),squeeze(tr(:,:,:,3)));
    foc_to_polarcenter=tr;
    tr(:,:,:,1)=tr(:,:,:,1)+polar_center(1); tr(:,:,:,2)=tr(:,:,:,2)+polar_center(2); tr(:,:,:,3)=tr(:,:,:,3)+polar_center(3);
    
else
    
    tr_temp1=squeeze(tr(:,:,1,:));
    if median_filter
        for ff=1:no_fish
            for d=1:no_dim
                tr_temp1(:,ff,d)=medfilt1(tr(:,ff,ff,d),median_filter_order);
            end
        end
    end
    
    if no_fish==1
        tr_temp=reshape(tr,[no_frames no_fish no_dim]);
    else
        tr_temp=squeeze(tr(:,:,1,:));
    end
    probIdCheck = any(any(squeeze(all(isnan(tr_temp))),2));
    if ~isempty(smooth_method) && ~isempty(smooth_degree) && isnumeric(smooth_degree) && ~probIdCheck && ...
            ~strcmpi(smooth_method,'none')
        [tr_temp, tol] = idSocial_auxiliaries_smoothTrajectories(tr,smooth_method,smooth_degree,max_deviation,noise_mean);

%         switch smooth_method
%             case {'moving', 'lowess','loess','sgolay','rlowess','rloess'}
%                 if strcmpi('smooth_method','sgolay')
%                     keyboard
%                 end
%                 for ff=1:no_fish
%                     for d=1:no_dim
%                         tr_temp(:,ff,d)=smooth(tr(:,ff,ff,d),smooth_degree,smooth_method);
%                     end
%                 end
%             case 'adaptive_moving'
%                 [tr_temp, mavgwidth] = idSocial_auxiliaries_smoothTrajectoryAdaptiveMovAvg(tr,smooth_degree);
%             case 'adaptive_moving_acc'
% %                 keyboard
%                 
%                 [tr_temp, mavgwidth] = idSocial_auxiliaries_smoothAcc2Vel2TrajectoryAdaptiveMovAvg(tr,smooth_degree);
%             case 'smoothing_spline'
% %                 keyboard
%                 tr_temp = idSocial_auxiliaries_splineSmoothing(tr,max_deviation,[],smooth_spline_degree);
%             case 'smoothing_spline_acc'
% %                 keyboard
%                 tr_temp = idSocial_auxiliaries_smoothAcc2Vel2TrajectoryAdaptiveSpline(tr,max_deviation);
%             case 'smoothing_spline_adaptive'
%                 [tr_temp, tol] = ...
%                     idSocial_adaptiveTrajectorySmoothing(tr,[],noise_mean,max_deviation,smooth_method);
% 
%             case {'DouglasPeucker','iterative_end_point','split_and_merge'}
%                 tr_temp = idSocial_auxiliaries_smoothTrajectoryDouglasPeucker(tr,smooth_degree);
%             case 'moving_separate'
%                  tr_temp = idSocial_auxiliaries_smoothTrVelAccSep(tr,smooth_degree,'moving');
%         end
    elseif ~isempty(smooth_method) && ~isempty(smooth_degree) && ischar(smooth_degree) && ...
            strcmpi(smooth_degree,'auto') || strcmpi(smooth_degree,'automatic')
    end
    if probIdCheck
        warning([mfilename ': Empty trajectoy, might be due to filtering (assignment probabilities etc.)'])
    end
end

if focalReconstructionFlag
    
    [tr_temp, good_frames]= ...
        idSocial_auxiliaries_idTrackerNoIdentities(tr_temp,probtrayectorias,options.focalReconstruction_minProbIdentity4VelCalcs, ...
        options.focalReconstruction_minProbIdentityFocal, ...
        options.focalReconstruction_minProbIdentityNeighbor);
    duration_total = sum(good_frames);
    duration_trial = no_frames;
    if duration_total<duration_trial
        warning([mfilename ': % of frames for all focal after reconstruction: ' num2str(duration_total)])
    else
        disp([mfilename ': total number of frames combining all focals after reconstruction: ' num2str(duration_total)])
    end
end

if ~isempty(interp_mode) && interpolate_trajectories
    tr=zeros(no_frames,no_fish,no_fish,no_dim);
    try
        tr(:,:,:,1:no_dim)=reshape(repmat(idSocial_interpolateTrajectories(tr_temp,interp_mode),[1,no_fish,1]),no_frames,no_fish,no_fish,no_dim);
    catch
        keyboard
    end
else
    if isnumeric(tr_temp)
        tr=reshape(repmat(tr_temp,[1,no_fish,1]),no_frames,no_fish,no_fish,no_dim);
%     elseif iscell(tr_temp)
%         for k=1:numel(tr_temp)
%             
%         end
    end
end

% end


% end

%% Generate structure for output
movementdata.options=options;
if isnumeric(tr_temp)
    movementdata.trajectory=tr(:,:,:,1:no_dim_orig);
else
     movementdata.trajectory=tr_temp;
end
if strcmpi(smooth_method,'adaptive_moving') || strcmpi(smooth_method,'adaptive_moving_acc')
    movementdata.moving_average_width = mavgwidth;
end
if disp_options
    disp(['[' mfilename '] ' 'Trajectories are ready! (It took ' num2str(toc) ' seconds to prepare it)'])
end
