function [foc_to_nb_vec_rot,vel_rot,acc_rot,turnmatr]=idSocial_auxiliaries_rotateToFocalSystem(trajectory,vel,acc)


no_fish = size(trajectory,2);
no_frames = size(trajectory,1);

if ~iscell(trajectory) && ( nargin < 2 || isempty(vel))
    vel = NaN(size(trajectory));
    vel(1:end-1,:,:,:) = diff(trajectory,1,1);
end

if ~iscell(trajectory) && (nargin < 3 || isempty(acc))
    acc = NaN(size(trajectory));
    acc(1:end-2,:,:,:) = diff(trajectory,2,1);
end

if iscell(trajectory)
    [trajectory,vel,acc,no_frames,no_fish,no_dim] = idSocial_auxiliaries_formatInputTrajectory(trajectory);
end
rand_check = idSocial_auxiliaries_trRandCheck(trajectory);

xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation
rotate_z = false;

trIsNaN = isnan(trajectory);
trIsNaN = any(trIsNaN,4);

vel_magn=sqrt(nansum(vel.^2,4));
vel_magn(trIsNaN) = NaN;
vel_magn(cat(1,trIsNaN(2:end,:,:),true(1,no_fish,no_fish))) = NaN;
vel_norm=bsxfun(@rdivide,vel,vel_magn);
vel_norm(repmat(vel_magn==0,[1,1,1,size(trajectory,4)]))=0;

% vel_norm will be exclusively used to determine the
% movement direction as a proxy for body orientation. If 
% vel_norm = 0, we assume that orientation has not changed
% and interpolate:

% b=~isnan(Y);
% c=cumsum(b);
% d=Y(b);
% Y2 = d(c)

if any(vel_magn(:)==0)
    for fr = 2:no_frames
        for ff = 1:no_fish
            for nf = 1:no_fish
                if ~rand_check(ff,nf)
                    if vel_magn(fr,ff,nf) == 0
                        vel_norm(fr,ff,nf,:) = vel_norm(fr-1,ff,nf,:)/1000;
                    end
                end
            end
        end
    end
end

foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,3);
ang_foc_vel_norm_and_nb_vel_norm=NaN(no_frames,no_fish,no_fish);



for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,nf,ff,:)=squeeze(trajectory(:,nf,ff,:)-trajectory(:,ff,ff,:));
            ang_foc_vel_norm_and_nb_vel_norm(:,ff,nf)=...
                mod(atan2(vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,2)-vel_norm(:,ff,ff,1).*vel_norm(:,nf,ff,2),vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,1)+vel_norm(:,nf,ff,2).*vel_norm(:,ff,ff,2)),2*pi);
        end
    end
end
% trajectory(:,:,:,2)=-trajectory(:,:,2); % This is a dirty trick in order to make the orientation of the coordinate system coherent with the videos

ang_vel_norm_z_plane_and_y_axis=...
    mod(atan2(vel_norm(:,:,:,1)*yax(2)-yax(1)*vel_norm(:,:,:,2),vel_norm(:,:,:,1)*yax(1)+vel_norm(:,:,:,2)*yax(2)),2*pi);

clear vel_norm vel_magn

%% Measures which depend on relative position/orientation of foc and nb: Rotate vectors
a1=NaN(no_frames,no_fish);
for ff=1:no_fish
    a1(:,ff)=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,ff,ff,:));
end
a1=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,:,1,:));
if rotate_z
    a2=NaN(no_frames,no_fish);
    for ff=1:no_fish
        a2(:,ff)=squeeze(ang_vel_norm_and_xy_plane(:,ff,ff,:));
    end
else
    a2=zeros(size(a1));
end
clear ang_vel_norm_z_plane_and_y_axis

a3=zeros(size(a1));
turnmatr=NaN(3,3,size(a1,1),no_fish);

turnmatr(1,1,:,:)=cos(a3).*cos(a1);
turnmatr(1,2,:,:)=-cos(a2).*sin(a1)+cos(a1).*sin(a3).*sin(a2);
turnmatr(1,3,:,:)=sin(a2).*sin(a1)+cos(a2).*sin(a3).*cos(a1);

turnmatr(2,1,:,:)=cos(a3).*sin(a1);
turnmatr(2,2,:,:)=sin(a2).*sin(a3).*sin(a1)+cos(a2).*cos(a1);
turnmatr(2,3,:,:)=-sin(a2).*cos(a1)+cos(a2).*sin(a3).*sin(a1);

turnmatr(3,1,:,:)=-sin(a3);
turnmatr(3,2,:,:)=sin(a2).*cos(a3);
turnmatr(3,3,:,:)=cos(a2).*cos(a3);

foc_to_nb_vec_rot=NaN(no_frames,no_fish,no_fish,3);
acc_rot=NaN(no_frames,no_fish,no_fish,3);
vel_rot=NaN(no_frames,no_fish,no_fish,3);

% vel_trans=NaN(size(vel_rot));
% vel_trans(:,:,:,1)=vel_rot(:,:,:,1);

for ff=1:no_fish
    for nf=1:no_fish
        %         if ff~=nf
        % Idx (:,ff,nf,:) means: Corresponding value (vel,acc,..) for focal ff 'seen by' neighbor nf;
        % In case of foc_to_nb_vec, vel and acc, 'seen
        % by' refers to neighbor-depend filtering (e.g.,
        % value at (:,ff,nf,:) is NaN because neighbor moves slower than
        % the lower threshold speed).
        % For rotated values, (:,ff,nf,:) means: vector
        % for focal ff in the system where neighbor velocity
        % vector is parallel to the y-axis
        if ~rand_check(ff,nf)
            foc_to_nb_vec_rot(:,nf,ff,:)=...
                [sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
                sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
                sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
            vel_rot(:,ff,nf,:)=...
                [sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
                sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
                sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
            acc_rot(:,ff,nf,:)=...
                [sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
                sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
                sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
%             vel_trans(:,ff,nf,2)=vel_rot(:,ff,nf,2)-vel_rot(:,ff,ff,2);
            %         end
        end
    end
end

