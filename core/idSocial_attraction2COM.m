function [distChangeVsDist, distChangeVsDistCell, distancedistribution,distance]= idSocial_attraction2COM(tr,edges,framerate,bodylength,method,speed_normalization,cut_percentile,individual_cm)
% Calculates distance distributions

no_fish=size(tr,2);
no_frames=size(tr,1);
no_dim=size(tr,4);

if nargin<2 || isempty(edges)
    edges = [];
end
if nargin<3 || isempty(framerate)
    framerate=1;
end
if nargin<4 || isempty(bodylength)
    bodylength=1;
end
if nargin<5 || isempty(method)
    method = @nanmean;
end

if nargin<5 || isempty(speed_normalization)
    speed_normalization = false;
end
if nargin<6 || isempty(cut_percentile)
    cut_percentile = 0;
end

if nargin<7 || isempty(individual_cm)
    individual_cm = false;
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

if ~individual_cm
    cmass =nanmean(tr,2);

    dist2com = sqrt((tr(:,:,1)-repmat(cmass(:,:,1),[1 no_fish 1])).^2+(tr(:,:,2)-repmat(cmass(:,:,2),[1 no_fish 1])).^2)/bodylength;
    Ddist2com = cat(1,diff(dist2com,1,1),NaN(1,no_fish))*framerate;

    dist2com2 = sqrt((tr(2:end,:,1)-repmat(cmass(1:end-1,:,1),[1 no_fish 1])).^2+(tr(2:end,:,2)-repmat(cmass(1:end-1,:,2),[1 no_fish 1])).^2)/bodylength;
    dist2com2 = cat(1,dist2com2,NaN(1,no_fish));
    Ddist2com = (dist2com2-dist2com)*framerate;
else
    cmass = NaN(no_frames,no_fish);
    dist2com = NaN(no_frames,no_fish);
    dist2com2 = NaN(no_frames,no_fish);
    Ddist2com = NaN(no_frames,no_fish);
    for ff = 1:no_fish
        cmass(:,ff) = nanmean(tr(:,setxor(1:no_fish,ff),:),2);
        dist2com(:,ff) = sqrt((tr(:,:,1)-repmat(cmass(:,:,1),[1 no_fish 1])).^2+(tr(:,:,2)-repmat(cmass(:,:,2),[1 no_fish 1])).^2)/bodylength;
        Ddist2com(:,ff) = cat(1,diff(dist2com(:,ff),1,1),NaN(1,1))*framerate;
        
        dist2com2(:,ff) = sqrt((tr(2:end,ff,1)-repmat(cmass(1:end-1,ff,1),[1 1 1])).^2+(tr(2:end,ff,2)-repmat(cmass(1:end-1,ff,2),[1 1 1])).^2)/bodylength;
        dist2com2(:,ff) = cat(1,dist2com2(:,ff),NaN(1,1));
        Ddist2com(:,ff) = (dist2com2(:,ff)-dist2com(:,ff))*framerate;
    end
    
end
% D2dist2com = cat(3,diff(dist2com,2,3),NaN(size(dist2com,1),size(dist2com,2),2,no_fish))*framerate;

% figure; hist(Ddist2com(:),200)
% figure; plot(Ddist2com(:,1))
% hold on; plot(smooth(Ddist2com(:,1),30,'moving'))
for ffish = 1:no_fish
    Ddist2com(:,ffish) = smooth(Ddist2com(:,ffish),30,'moving');
end
Ddist2com(Ddist2com(:)<prctile(Ddist2com(:),cut_percentile) | Ddist2com(:)>prctile(Ddist2com(:),100-cut_percentile))=NaN;

allpresent = all(~isnan(dist2com),2);
%Extend this to next index since Ddist is calculated from idx and idx+1:
Dallpresent = allpresent;% vertcat(allpresent(1:end-1) & allpresent(2:end),false);

dist2com(~Dallpresent,:) = NaN(sum(~Dallpresent),no_fish);
Ddist2com(~Dallpresent,:) = NaN(sum(~Dallpresent),no_fish);
disp([mfilename ': % of frames with all animals present: ' num2str(sum(allpresent)/no_frames*100)])
% You cannot really calculate the center of mass not knowwing the position
% of all animals.

if isempty(edges)
    edges = 0:.5:max(dist2com(:));
end

no_bins=size(edges,2)-1;
distancedistribution=NaN(no_fish,no_fish,no_bins);
distance = NaN(no_fish,no_fish,no_frames);
bins = NaN(no_frames,no_fish);

for ffish = 1:no_fish
    distance(ffish,ffish,:)=squeeze(dist2com(:,ffish));
    [hi,bins(:,ffish)] = histc(squeeze(dist2com(:,ffish)),edges);
    distancedistribution(ffish,ffish,:) = hi(1:end-1);
end

distChangeVsDist =  NaN(no_fish,numel(edges));

for ffish = 1:no_fish
    good_bins = bins(:,ffish)~=0;
    
    if speed_normalization
        vexp = diff(squeeze(tr(:,ffish,:)),1,1);
        vexpMagn = sqrt(sum(vexp.^2,2));
        speed = nanmean(vexpMagn)/bodylength*framerate;
        map_val = squeeze(Ddist2com(good_bins,ffish))/speed;
    else
        map_val = squeeze(Ddist2com(good_bins,ffish));
    end
    
    tt = accumarray(squeeze(bins(good_bins,ffish)),map_val,[numel(edges) 1],method);
    distChangeVsDist(ffish,:) = tt';
end

distChangeVsDistCell =  cell(no_fish,no_fish,numel(edges));

try
for ffish = 1:no_fish
    good_bins = bins(:,ffish)~=0;
    if speed_normalization
        vexp = diff(squeeze(tr(:,ffish,:)),1,1);
        vexpMagn = sqrt(sum(vexp.^2,2));
        speed = nanmean(vexpMagn)/bodylength*framerate;
        map_val = squeeze(Ddist2com(good_bins,ffish))/speed;
    else
        map_val = squeeze(Ddist2com(good_bins,ffish));
    end
    if any(good_bins)
        distChangeVsDistCell (ffish,ffish,:) = accumarray(squeeze(bins(good_bins,ffish)),map_val,[numel(edges) 1],@(x) {x});
    end
end
catch
    keyboard
end

 distancedistribution=mat2cell(distancedistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));


