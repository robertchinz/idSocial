function [distancedistribution, distance_median,distance_mode,inbins, distance,funcinfo]=idSocial_distanceFromCenterDistribution(tr,edges,bodylength,method,kds_bandwidth, ...
    distance_units,arena_center,arena_radius,normalization)
% Calculates mutual distances between group members

if nargin<3 || isempty(bodylength)
    bodylength=1;
end

if nargin<4 || isempty(method)
    method='hist';
end

if nargin<5 || isempty(kds_bandwidth)
    kds_bandwidth=.5;
end


if nargin<6 || isempty(distance_units)
    distance_units='pixels';
end
if nargin<7 || isempty(arena_center)
    arena_center=[0 0];
end
if nargin<8 || isempty(arena_radius)
    arena_radius = 1;
end
if nargin<9 || isempty(normalization)
    normalization='none';
end

[tr,~,~,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

trIsNaN = isnan(tr);
if no_fish==1
    no_fishOrig=1;
    no_fish=2;
    tr_temp = NaN(no_frames,2,2,no_dim);
    tr_temp(:,1,1,:) = tr;
    tr_temp(:,1,2,:) = tr;
    tr_temp(:,2,1,1:2) = repmat(arena_center',[no_frames,1, 1,1]);
    tr_temp(:,2,2,1:2) = repmat(arena_center',[no_frames,1, 1,1]);
    tr_temp(:,2,2,3) = 0;
     tr_temp(:,2,1,3) = 0;
    
    tr = tr_temp;
    clear tr_temp;
else
    no_fishOrig=no_fish;

    % Substitute neighbors by center
    for ff = 1:no_fish
        for nf = 1:no_fish
            if ff~=nf
                tr(:,nf,ff,1:2) = repmat(arena_center',[no_frames,1, 1,1]);
            end
        end
    end
end
tr(trIsNaN) = NaN; % To recover focal-neighbor filtering.

% ctr = arena_center;
% dist2ctr = NaN(no_frames,no_fish);
% if invertY
%     ctr(2)=ctr(2)*nanmean(sign(tr(:,1,1,2)));
% end
% for ff = 1:no_fish
%     dist2ctr(:,ff) = sqrt(sum((squeeze(tr(:,ff,ff,1:2))-repmat(ctr',[no_frames 1])).^2,2))/arena_radius;
%     
% end

foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,no_dim);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,ff,nf,:)=squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
        end
    end
end
distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,4));

if strcmpi(distance_units,'arena_radius')
    distance_focal_neighbour = distance_focal_neighbour/arena_radius;
elseif strcmpi(distance_units,'bodylength')
    distance_focal_neighbour = distance_focal_neighbour./bodylength;
end

distance_focal_neighbour=permute(distance_focal_neighbour,[2,3,1]);
distance_focal_neighbour(logical(repmat(eye(no_fish,no_fish),[1 1 no_frames])))=NaN;
%%
no_bins=size(edges,2)-1;
distancedistribution=NaN(no_fish,no_fish,no_bins);
distance_median=NaN(no_fish,no_fish);
distance=NaN(no_fish,no_fish,no_frames);
distance_mode=NaN(no_fish,no_fish);
inbins=cell(no_fish,no_fish,no_bins);

locs = [];
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            val=squeeze(distance_focal_neighbour(ff,nf,:));
            
            if strcmpi(method,'hist')
                
                [hitemp,bins]=histc(val,edges);
                
                if ~isempty(edges)
                    switch normalization
                        case 'density'
                            distancedistribution(ff,nf,:)=hitemp(1:end-1)/no_frames/(edges(2)-edges(1));
                        case 'no_frames'
                            distancedistribution(ff,nf,:)=hitemp(1:end-1)/no_frames;
                        case 'none'
                            distancedistribution(ff,nf,:)=hitemp(1:end-1);
                    end
                end
                if ~all(bins==0)
                    inbins(ff,nf,:)=accumarray(bins(~isnan(val) & bins>0),ones(1,sum(~isnan(val)& bins>0)),[no_bins 1],@(x) {x});
                    
%                     inbins(ff,nf,:)=accumarray(bins(~isnan(val) & bins>0),val(~isnan(val)& bins>0),[no_bins 1],@(x) {x});
                end
            elseif strcmpi(method,'ksdensity_gauss') || strcmpi(method,'ksdensity_epanechnikov')...
                    || strcmpi(method,'ksdensity_triangular')
%                 keyboard
                edges_kds=edges(1:end-1);%(1:end-1)+(edges(2)-edges(1))/2;                
                [hitemp,prob] = idSocial_auxiliaries_kerneldensity(val,edges_kds,kds_bandwidth,method);
%                 figure; plot(edges_kds,hitemp)
                distancedistribution(ff,nf,:)=hitemp;
                for k=1:no_bins
                    inbins{ff,nf,k}=prob(~isnan(prob(:,k)) & prob(:,k)>0,k);
                end
                [~,locs] = idSocial_auxiliaries_findpeaks(hitemp);

            end
            distance_median(ff,nf)=nanmedian(val);
            distance(ff,nf,:)=val;
          
%             [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
%            [~,locs] = findpeaks(hitemp);

            if ~isempty(locs)
                distance_mode(ff,nf)=edges(min(locs));
            end
            
                   
            
        end
    end
end

distancedistribution=mat2cell(distancedistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));


%%
fstack=dbstack(3);
funcinfo.Function = mfilename;
funcinfo.callerFunction = fstack(1).name;
funcinfo.no_fish = no_fish;
funcinfo.XTick = edges(1:end-1);
funcinfo.XTickLabel = cellfun(@(x) strtrim(x),cellstr((num2str(edges(1:end-1)')))','UniformOutput',false);


if no_fishOrig==1
    distancedistribution = distancedistribution(1,2,:);
    inbins = inbins(1,2,:);
%     funcinfo = funcinfo(1,2,:);
    distance_median = distance_median(1,2,:);
    distance_mode = distance_mode(1,2,:);
    distance = distance(1,2,:);
end