function [distancedistribution,distancedistributionCell,distance,distanceMax]= idSocial_distance2COM(tr,edges,bodylength,method)
% Calculates distance distributions

no_fish=size(tr,2);
no_frames=size(tr,1);
no_dim=size(tr,4);

if nargin<2 || isempty(edges)
    edges = [];
end

if nargin<3 || isempty(bodylength)
    bodylength=1;
end
if nargin<4 || isempty(method)
    method = @nanmean;
end


if ischar(method)
    switch method
        case 'mean'
            method = @nanmean;
        case 'median'
            method = @nanmedian;
    end
end
% bodylength = 1;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,[],[],false);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

cmass =nanmean(tr,2);

dist2com = sqrt((tr(:,:,1)-repmat(cmass(:,:,1),[1 no_fish 1])).^2+(tr(:,:,2)-repmat(cmass(:,:,2),[1 no_fish 1])).^2)/bodylength;

allpresent = all(~isnan(dist2com),2);
%Extend this to next index since Ddist is calculated from idx and idx+1:
Dallpresent = allpresent;% vertcat(allpresent(1:end-1) & allpresent(2:end),false);

dist2com(~Dallpresent,:) = NaN(sum(~Dallpresent),no_fish);
dist2comMax = max(dist2com,[],2);
% Ddist2com(~Dallpresent,:) = NaN(sum(~Dallpresent),no_fish);
disp([mfilename ': % of frames with all animals present: ' num2str(sum(allpresent)/no_frames*100)])
% You cannot really calculate the center of mass not knowwing the position
% of all animals.

if isempty(edges)
    edges = 0:.5:max(dist2com(:));
end

no_bins=size(edges,2)-1;
distancedistribution=NaN(no_fish,no_fish,no_bins);
distance = NaN(no_fish,no_fish,no_frames);
distanceMax = NaN(no_fish,no_fish,no_frames);
distanceMax(1,1,:)=dist2comMax;
bins = NaN(no_frames,no_fish);

for ffish = 1:no_fish
    distance(ffish,ffish,:)=squeeze(dist2com(:,ffish));
    [hi,bins(:,ffish)] = histc(squeeze(dist2com(:,ffish)),edges);
    distancedistribution(ffish,ffish,:) = hi(1:end-1);
end

distancedistributionCell=mat2cell(distancedistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));


