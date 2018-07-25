function [tr_smooth, dist2o, delta2o]= idSocial_auxiliaries_smoothAcc2Vel2Trajectory(trajectories,sm_param,method)

if nargin < 2 || isempty(sm_param)
    sm_param = 1;
end
if nargin < 3 || isempty(method)
    method = 'moving';
end

no_fish = size(trajectories,2);
no_frames = size(trajectories,1);
no_dim = size(trajectories,3);

vel = vertcat(diff(trajectories,1,1),NaN(1,no_fish,size(trajectories,3)));
acc = vertcat(diff(trajectories,2,1),NaN(2,no_fish,size(trajectories,3)));

acc2 = NaN(size(trajectories));

for ff = 1:no_fish
    for dm = 1:no_dim
        if strcmpi(method,'moving')
            acc2(:,ff,dm) = smooth(acc(:,ff,dm),sm_param,'moving');
        elseif strcmpi(method,'spline')
            
            tol = no_frames*sm_param^2;
            
            accnan = isnan(acc(:,ff,dm));
            idx = 1:size(trajectories,1); idx = idx(~accnan);
%             [ppaccx, p] = csaps(idx,squeeze(acc(~accnan,ff,dm)),sm_param);
            [ppaccx, ~, p] = spaps(idx,squeeze(acc(~accnan,ff,dm)),tol);
            acc2(~accnan,ff,dm) = fnval(idx,ppaccx);
            
        end
    end
end
acc=acc2;

acc4vel = NaN(size(trajectories));
acc4vel(2:end,:,:) = acc(1:end-1,:,:); acc4vel(1,:,:) = zeros(1,no_fish,no_dim);
acc2vel = cumsum(acc4vel,1);
v0 = -nanmean(acc2vel-vel);
acc2vel = acc2vel + repmat(v0,[no_frames 1 1]);

figure; plot(abs(squeeze(acc2vel(:,1,1))-squeeze(vel(:,1,1))))

vel4tr = NaN(size(trajectories));
vel4tr(2:end,:,:) = acc2vel(1:end-1,:,:); vel4tr(1,:,:) = zeros(1,no_fish,no_dim);
vel2tr = cumsum(vel4tr,1);
tr0 = -nanmean(vel2tr-trajectories);
vel2tr = vel2tr + repmat(tr0,[no_frames 1 1]);

figure; plot(abs(squeeze(vel2tr(:,2,1))-squeeze(trajectories(:,2,1))))
figure; plot((squeeze(vel2tr(:,2,1)))); hold on; plot(squeeze(trajectories(:,2,1)),'r');


tr_smooth = vel2tr;

delta2o = sqrt(nansum((acc2vel-acc).^2,3));
dist2o = sqrt(nansum((tr_smooth-trajectories).^2,3));
% 
idces = 1:size(trajectories,1);
figure; plot(squeeze(vel(idces,2,1)),'--'); hold on; %plot(squeeze(vel(idces,ff,2)),':b');
plot(acc2vel(idces,2,1),'--g','LineWidth',2)
figure; plot(squeeze(trajectories(idces,1,2)),'--'); hold on; %plot(squeeze(vel(idces,ff,2)),':b');
plot(vel2tr(idces,1,2),'--g','LineWidth',2)

