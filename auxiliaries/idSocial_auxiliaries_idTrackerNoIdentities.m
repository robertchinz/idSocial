% load('G:\trajectories.mat')
function [trajectories_new, good_frames] = idSocial_auxiliaries_idTrackerNoIdentities(trajectories,probtrajectories,prob4vel,prob4filter,prob4filterNeighbor)
% prob4vel = .9;
prob4acc = prob4vel;
% prob4filter = 0.9;
% prob4filterNeighbor = 0.1;

no_frames = size(probtrajectories,1);
no_fish = size(probtrajectories,2);
no_dims = size(trajectories,3);
ndims_orig = ndims(trajectories);

% if ndims_orig==4
%     trtemp = NaN(no_frames,no_fish,no_dims);
%     for ff = 1:no_fish
%         trtemp(:,ff,:) = trajectories(:,ff,:);
%     end
% end

vel = vertcat(sqrt(nansum(diff(trajectories,1,1).^2,3)),NaN(1,no_fish));
velGood = vel;
velGood(probtrajectories>prob4vel) = NaN;
disp([mfilename ': ' num2str(sum(probtrajectories(:)>prob4vel)) ' samples to calculate max. velocity.'])
velFilter = vel<prctile(velGood(:),95);

acc = vertcat(sqrt(nansum(diff(trajectories,2,1).^2,3)),NaN(2,no_fish));
accGood = acc;
accGood(probtrajectories>prob4acc) = NaN;
disp([mfilename ': ' num2str(sum(probtrajectories(:)>prob4acc)) ' samples to calculate max. acceleration.'])
accFilter = acc<prctile(accGood(:),95);
%
% figure; hist(velGood(:),200)
% figure; hist(accGood(accFilter(:)),200)
%
probFilter = false(size(probtrajectories));
probFilter(velFilter(:) & accFilter(:) & probtrajectories(:)>prob4filter) = true;

probFilterNeighbor = false(size(probtrajectories));
probFilterNeighbor(probtrajectories(:)>prob4filterNeighbor) = true;

 
trajectories_new = NaN(no_frames,no_fish,no_fish,no_dims);
good_frames = NaN(1,no_fish);
for ff=1:no_fish
    
    focalFilter = squeeze(probFilter(:,ff));
    neighborFilter = sum(squeeze(probFilterNeighbor(:,setxor(ff,1:no_fish)))>prob4filterNeighbor,2)==no_fish-1;
    trajectories_new(focalFilter & neighborFilter,:,ff,:) = trajectories(focalFilter & neighborFilter,:,:);
    good_frames(ff) = sum(focalFilter & neighborFilter);
end
good_frames_percentage = good_frames/no_frames*100;

if ndims_orig==3
    trtemp = NaN(no_frames,no_fish,no_dims);
    for ff = 1:no_fish
        trtemp(:,ff,:) = trajectories_new(:,ff,ff,:);
    end
    trajectories_new = trtemp;
end
% if ndims_orig==4
%     tr=zeros(no_frames,no_fish,no_fish,no_dim);
%     tr(:,:,:,1:no_dim)=reshape(repmat(trajectories_new,[1,no_fish]),no_frames,no_fish,no_fish,no_dim);
%     trajectories_new=tr;
% end
% [good_frames, all_members_frames,all_members_present_for_worst_focal ] = idSocial_auxiliaries_checkDataQuality(trajectories_new);
if min(good_frames_percentage)<70
    warning([mfilename ': % of frames for worst focal after reconstruction: ' num2str(min(good_frames_percentage))])
else
    disp([mfilename ': % of frames for worst focal after reconstruction: ' num2str(min(good_frames_percentage))])
end

%%
% fr=1:5000;
% ff=2;
% figure; 
% col = jet(no_fish);
% hold on
% for ff = 1:no_fish
%     plot(trajectories_new(fr,ff,ff,1),trajectories_new(fr,ff,ff,2),'o-','Color',col(ff,:))
% end
% axis equal