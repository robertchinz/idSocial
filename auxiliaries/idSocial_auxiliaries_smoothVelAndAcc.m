function [trcell, vel, acc] = idSocial_auxiliaries_smoothVelAndAcc(trajectories,sm_param,method)

no_fish = size(trajectories,2);

if nargin<2 || isempty(sm_param)
    sm_tr = 1;
    sm_vel = 1;
    sm_acc = 1;
else
    sm_tr = sm_param(1);
    sm_vel = sm_param(2);
    sm_acc = sm_param(3);
end

if nargin<3 || isempty(method)
    method = 'moving';
end


    

vel = vertcat(diff(trajectories,1,1),NaN(1,no_fish,size(trajectories,3)));
acc = vertcat(diff(trajectories,2,1),NaN(2,no_fish,size(trajectories,3)));

% disp('Smoothing parameter is so random!')

trX = NaN(size(trajectories,1),1);
trY = NaN(size(trajectories,1),1);
velx = NaN(size(trajectories,1),1);
vely = NaN(size(trajectories,1),1);
accx = NaN(size(trajectories,1),1);
accy = NaN(size(trajectories,1),1);
trcell = cell(no_fish,1);
for ff = 1:no_fish
    trnanx = isnan(trajectories(:,ff,1));
    trnany = isnan(trajectories(:,ff,2));
    idxx = 1:size(trajectories,1); idxx = idxx(~trnanx);
    idxy = 1:size(trajectories,1); idxy = idxy(~trnany);
    
    if strcmpi(method,'spline')
        pptrx = csaps(idxx,squeeze(trajectories(~trnanx,ff,1)),sm_tr);
        pptry = csaps(idxy,squeeze(trajectories(~trnany,ff,2)),sm_tr);
        trX(idxx) = fnval(idxx,pptrx);
        trY(idxy) = fnval(idxy,pptry);
    else
        trX = smooth(trajectories(:,ff,1),sm_tr,'moving');
        trY = smooth(trajectories(:,ff,2),sm_tr,'moving');
    end
    
    velnanx = isnan(vel(:,ff,1));
    velnany = isnan(vel(:,ff,2));
    idxx = 1:size(trajectories,1); idxx = idxx(~velnanx);
    idxy = 1:size(trajectories,1); idxy = idxy(~velnany);
    
    if strcmpi(method,'spline')
        ppvelx = csaps(idxx,squeeze(vel(~velnanx,ff,1)),sm_vel);
        ppvely = csaps(idxy,squeeze(vel(~velnany,ff,2)),sm_vel);
        velx(idxx) = fnval(idxx,ppvelx);
        vely(idxy) = fnval(idxy,ppvely);
    else
        velx = smooth(vel(:,ff,1),sm_vel,'moving');
        vely = smooth(vel(:,ff,2),sm_vel,'moving');
    end
    
    
    accnanx = isnan(acc(:,ff,1));
    accnany = isnan(acc(:,ff,2));
    idxx = 1:size(trajectories,1); idxx = idxx(~accnanx);
    idxy = 1:size(trajectories,1); idxy = idxy(~accnany);
    
    if strcmpi(method,'spline')
        ppaccx = csaps(idxx,squeeze(acc(~accnanx,ff,1)),sm_acc);
        ppaccy = csaps(idxy,squeeze(acc(~accnany,ff,2)),sm_acc);
        accx(idxx) = fnval(idxx,ppaccx);
        accy(idxy) = fnval(idxy,ppaccy);
    else
        accx = smooth(acc(:,ff,1),sm_acc,'moving');
        accy = smooth(acc(:,ff,2),sm_acc,'moving');
    end
    
    %         ppaccx = csaps(1:size(trajectories,1),squeeze(acc(:,ff,1)),.5);
    %         ppaccy = csaps(1:size(trajectories,1),squeeze(acc(:,ff,2)),.5);
    %         accx = fnval(1:size(trajectories,1),ppaccx);
    %         accy = fnval(1:size(trajectories,1),ppaccy);
    
    vel(1:size(trajectories,1),ff,1) = velx;
    vel(1:size(trajectories,1),ff,2) = vely;
    acc(1:size(trajectories,1),ff,1) = accx;
    acc(1:size(trajectories,1),ff,2) = accy;
    
    trcell{ff}.X = trX;%squeeze(trajectories(:,ff,1));
    trcell{ff}.Y = trY;%squeeze(trajectories(:,ff,2));
    trcell{ff}.velx = velx;
    trcell{ff}.vely = vely;
    trcell{ff}.accx = accx;
    trcell{ff}.accy = accy;
end

%% Try to reconstruct trajectory from smooth acc
idces = 1:size(trajectories,1);
ff=2;
accx4velx = NaN(size(trajectories,1),1);
accx4velx(2:end) = accx(2:end); accx4velx(1) = 0;%3.55;%3.363;%6.36999;%nanmean(cumsum(accx));
accx2velx = cumsum(accx4velx);
vx0 = -nanmean(accx2velx-vel(idces,ff,1));
accx4velx(2:end) = accx(2:end); accx4velx(1) = vx0;
accx2velx = cumsum(accx4velx);
figure; plot(squeeze(vel(idces,ff,1)),'--'); hold on; %plot(squeeze(vel(idces,ff,2)),':b');
plot(accx2velx(idces),'--g','LineWidth',2)

velx4trx = NaN(size(trajectories,1),1);
velx4trx(2:end) = accx2velx(2:end); velx4trx(1) = 0;
velx2trx = cumsum(velx4trx);
trx0 = -nanmean(velx2trx-squeeze(trajectories(idces,ff,1)));

velx4trx(2:end) = accx2velx(2:end); velx4trx(1) = trx0;
velx2trx = cumsum(velx4trx);

figure; plot(squeeze(trajectories(idces,ff,1)),'--'); hold on;% plot(squeeze(trajectories(idces,1,2)),':b');
plot(velx2trx(idces),'--g','LineWidth',2)

%
%         accy4vely = NaN(size(trajectories,1),1);
%         accy4vely(2:end) = accy(2:end); accy4vely(1) = -3;%nanmean(cumsum(accx));
%         accy2vely = cumsum(accy4vely);
%          figure; plot(squeeze(vel(idces,1,2)),'--'); hold on; plot(squeeze(vel(idces,1,2)),':b');
%          plot(accy2vely(idces),'--g','LineWidth',2)
%          

%          
%          vely4try = NaN(size(trajectories,1),1);
%          vely4try(2:end) = accy2vely(2:end); vely4try(1) = trajectories(1,1,1);%nanmean(cumsum(accy));
%          vely2try = cumsum(vely4try);
%          figure; plot(squeeze(trajectories(idces,1,1)),'--'); hold on; plot(squeeze(trajectories(idces,1,2)),':b');
%          plot(vely2try(idces),'--g','LineWidth',2)