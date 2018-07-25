function [good_frames_worst_indiv, good_frames_all_present, all_members_present_for_worst_focal,all_members_present_for_best_focal] = idSocial_auxiliaries_checkDataQuality(tr)

no_frames = size(tr,1);

if iscell(tr) || isstruct(tr)
    [tr,~,~,no_frames] = idSocial_auxiliaries_formatInputTrajectory(tr);
end

if ndims(tr)==4
    rand_check = idSocial_auxiliaries_trRandCheck(tr);
    
    % Percentage of good frames for worst individual (all the others in the
    % video have a higher percentage)
    good_frames_worst_indiv=squeeze(sum(~isnan(tr(:,:,:,1)),1)./no_frames);
    good_frames_worst_indiv=round(min(good_frames_worst_indiv(logical(~rand_check)))*100);
    
    % Old good_frames_all_present: Checks if all FOCALS (ff,ff) are present simultaneously,
    % but misses if (ff,ff) is filtered out(in this case due to "focal in center filter").
    % % Percentage of frames in which all focals are present (before filtering!! Checks only "the diagonal")
    % trPerm=squeeze(permute(tr(:,:,:,1),[1,4,2,3]));
    % trPerm = trPerm(:,logical(eye(size(tr,2),size(tr,2))));
    % good_frames_all_present = sum(sum(~isnan(trPerm),2)==size(tr,2))./size(tr,1);
    % good_frames_all_present = good_frames_all_present*100;
    
    rand_check = reshape(rand_check,[1 size(rand_check,1) size(rand_check,1)]);
    rand_check = logical(repmat(rand_check,[no_frames 1 1]));
    
    
    % New good_frames_all_present: Checks if focal coordinate exists at all.
    trCollapse = squeeze(nanmean(tr(:,:,:,1),3));
    good_frames_all_present = sum(sum(~isnan(trCollapse)|all(rand_check,3),2)==size(trCollapse,2))./no_frames;
    good_frames_all_present = good_frames_all_present*100;
    
    % Percentage of frames in which all neighbors are present for worst focal (all the others in the
    % video have a higher percentage). Similar to 'good_frames_all_present', but can be affected by filtering.
    % all_members_present_for_worst_focal = min(sum(squeeze(sum(~isnan(tr(:,:,:,1)),2)==size(tr,2)),1))/no_frames*100;
    % all_members_present_for_best_focal = max(sum(squeeze(sum(~isnan(tr(:,:,:,1)),2)==size(tr,2)),1))/no_frames*100;
    all_members_present_for_worst_focal = min(sum(squeeze(sum(~isnan(tr(:,:,:,1))|rand_check,2)==size(tr,2)),1))/no_frames*100;
    all_members_present_for_best_focal = max(sum(squeeze(sum(~isnan(tr(:,:,:,1))|rand_check,2)==size(tr,2)),1))/no_frames*100;
    
elseif ndims(tr) == 3
    rand_check = idSocial_auxiliaries_trRandCheck(tr);

    
    
    % Percentage of good frames for worst individual (all the others in the
    % video have a higher percentage)
    good_frames_worst_indiv=squeeze(sum(~isnan(tr(:,:,:,1)),1)./no_frames);
    good_frames_worst_indiv=round(min(good_frames_worst_indiv(logical(~rand_check)))*100);
    
    % Old good_frames_all_present: Checks if all FOCALS (ff,ff) are present simultaneously,
    % but misses if (ff,ff) is filtered out(in this case due to "focal in center filter").
    % % Percentage of frames in which all focals are present (before filtering!! Checks only "the diagonal")
    % trPerm=squeeze(permute(tr(:,:,:,1),[1,4,2,3]));
    % trPerm = trPerm(:,logical(eye(size(tr,2),size(tr,2))));
    % good_frames_all_present = sum(sum(~isnan(trPerm),2)==size(tr,2))./size(tr,1);
    % good_frames_all_present = good_frames_all_present*100;
    
    
    rand_check = reshape(rand_check,[1 size(rand_check,2)]);
    rand_check = logical(repmat(rand_check,[no_frames 1 1]));
    
    
    % New good_frames_all_present: Checks if focal coordinate exists at all.
    trCollapse = squeeze(nanmean(tr(:,:,1),2));
    good_frames_all_present = sum(sum(~isnan(trCollapse)|all(rand_check,2),2)==size(trCollapse,2))./no_frames;
    good_frames_all_present = good_frames_all_present*100;
    
    % Percentage of frames in which all neighbors are present for worst focal (all the others in the
    % video have a higher percentage). Similar to 'good_frames_all_present', but can be affected by filtering.
    % all_members_present_for_worst_focal = min(sum(squeeze(sum(~isnan(tr(:,:,:,1)),2)==size(tr,2)),1))/no_frames*100;
    % all_members_present_for_best_focal = max(sum(squeeze(sum(~isnan(tr(:,:,:,1)),2)==size(tr,2)),1))/no_frames*100;
    all_members_present_for_worst_focal =  good_frames_all_present;% min(sum(squeeze(sum(~isnan(tr(:,:,:,1))|rand_check,2)==size(tr,2)),1))/no_frames*100;
    all_members_present_for_best_focal =  good_frames_all_present;
    
    
    
    
end

% Note: all_members_present_for_worst_focal
% can be GREATER than good_frames_all_present,
% because good_frames_all_present considers only
% frames where all focals (the diagonals) are
% present, which can be filtered e.g. by the
% focal velocity filter.
% all_members_present_for_worst_focal
% however looks if there is one focal with all
% its neighbors present, and individuals can
% appear as neighbors even though they do not
% appear as focals.