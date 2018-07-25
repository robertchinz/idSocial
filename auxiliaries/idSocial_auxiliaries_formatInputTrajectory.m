function [trajectory,vel,acc,no_frames,no_fish,no_dim] = idSocial_auxiliaries_formatInputTrajectory(trajectory,invertY,spline_derivs,addNbDim_flag,add3rdDim_flag)

if nargin < 2 || isempty(invertY)
    invertY=false;
end
invertY=false; % invertY is not longer used. The "y-axis-inversion-problem" is now handled in idSocial_auxiliaries_rotateToFocalSystem
% by (temporarily) inverting the y-axis to calculate relative positions. 

if nargin < 3 || isempty(spline_derivs)
    spline_derivs=false;
end

if nargin < 4 || isempty(addNbDim_flag)
    addNbDim_flag=true;
end
if nargin < 5 || isempty(add3rdDim_flag)
    add3rdDim_flag=false;
end
if iscell(trajectory) && isstruct(trajectory{1}) % from splines or tr structs
    [trajectory, vel, acc] = idSocial_auxiliaries_spline2arrays(trajectory,spline_derivs);
    no_fish=size(trajectory,2);
    no_frames=size(trajectory,1);
    
    if invertY
        if ndims(trajectory)==4
            trajectory(:,:,:,2)=-trajectory(:,:,:,2);
            vel(:,:,:,2) = -vel(:,:,:,2);
            acc(:,:,:,2) = -acc(:,:,:,2);
        elseif ndims(trajectory)==3
            trajectory(:,:,2)=-trajectory(:,:,2);
            vel(:,:,2) = -vel(:,:,2);
            acc(:,:,2) = -acc(:,:,2);
        end
    end
    
elseif  isstruct(trajectory) && isfield(trajectory,'Tr') && ...
        isfield(trajectory,'Vel') &&  isfield(trajectory,'Acc')
    vel = trajectory.Vel;
    acc = trajectory.Acc;
    trajectory = trajectory.Tr;
    no_fish=size(trajectory,2);
    no_frames=size(trajectory,1);
elseif  isstruct(trajectory) && isfield(trajectory,'Tr') && ...
        ~isfield(trajectory,'Vel') &&  ~isfield(trajectory,'Acc')
    trajectory = trajectory.Tr;
    no_fish=size(trajectory,2);
    no_frames=size(trajectory,1);
    vel=NaN(size(trajectory));
    acc=NaN(size(trajectory));
    vel(1:no_frames-1,:,:,:)=diff(trajectory,1,1);
    acc(1:no_frames-2,:,:,:)=diff(trajectory,2,1);
    
else
    no_fish=size(trajectory,2);
    no_frames=size(trajectory,1);
    vel=NaN(size(trajectory));
    acc=NaN(size(trajectory));
    
    if invertY
        if ndims(trajectory)==4
            trajectory(:,:,:,2)=-trajectory(:,:,:,2);
        elseif ndims(trajectory)==3
            trajectory(:,:,2)=-trajectory(:,:,2);
        end
    end
    
    vel(1:no_frames-1,:,:,:)=diff(trajectory,1,1);
    acc(1:no_frames-2,:,:,:)=diff(trajectory,2,1);
end

%% Add nb dimension
if addNbDim_flag
    if ndims(trajectory)==3
        no_dim=size(trajectory,3);
        tr2=NaN(no_frames,no_fish,no_fish,no_dim);
        for ff=1:no_fish
            tr2(:,ff,:,:)=repmat(trajectory(:,ff,:),[1 no_fish 1]);
        end
        trajectory=tr2;
        clear tr2
    else
        no_dim=size(trajectory,4);
    end
    if ndims(vel)==3
        vel2=NaN(no_frames,no_fish,no_fish,no_dim);
        acc2=NaN(no_frames,no_fish,no_fish,no_dim);
        for ff=1:no_fish
            vel2(:,ff,:,:)=repmat(vel(:,ff,:),[1 no_fish 1]);
            acc2(:,ff,:,:)=repmat(acc(:,ff,:),[1 no_fish 1]);
        end
        vel=vel2;
        acc=acc2;
        clear vel2 acc2
    end
else
    if ndims(trajectory)==4
        no_dim=size(trajectory,4);
        tr2=NaN(no_frames,no_fish,no_dim);
        for ff=1:no_fish
            tr2(:,ff,:)=squeeze(trajectory(:,ff,ff,:));
        end
        trajectory=tr2;
        clear tr2
    else
        no_dim=size(trajectory,3);
    end
    if ndims(vel)==4
        vel2=NaN(no_frames,no_fish,no_dim);
        acc2=NaN(no_frames,no_fish,no_dim);
        for ff=1:no_fish
            vel2(:,ff,:)=squeeze(vel(:,ff,ff,:));
            acc2(:,ff,:)=squeeze(acc(:,ff,ff,:));
        end
        vel=vel2;
        acc=acc2;
        clear vel2 acc2
    end
end
%% Add 3rd dimension
if add3rdDim_flag
    if size(trajectory,ndims(trajectory))==2
        if ndims(trajectory)==4
            trtemp=zeros(no_frames,no_fish,no_fish,3);
            trtemp(:,:,:,1:no_dim)=trajectory;
        else
            trtemp=zeros(no_frames,no_fish,3);
            trtemp(:,:,1:no_dim)=trajectory;
        end
        trajectory=trtemp;
    end
    if size(vel,ndims(vel))==2
        if ndims(trajectory)==4
            veltemp=zeros(no_frames,no_fish,no_fish,3);
            veltemp(:,:,:,1:no_dim)=vel;
        else
            veltemp=zeros(no_frames,no_fish,3);
            veltemp(:,:,1:no_dim)=vel;
        end
        vel=veltemp;
    end
    if size(acc,ndims(acc))==2
        if ndims(trajectory)==4
            acctemp=zeros(no_frames,no_fish,no_fish,3);
            acctemp(:,:,:,1:no_dim)=acc;
        else
            acctemp=zeros(no_frames,no_fish,3);
            acctemp(:,:,1:no_dim)=acc;
        end
        acc=acctemp;
    end
    no_dim = 3;
    
end