function [speeddistribution, speed_median,speed_mode,inbins,speed,funcinfo]=idSocial_turningDistribution(tr,edges,framerate,bodylength,normalization,method,kds_bandwidth)
% Calculates turning acceleration distributions

if nargin<3 || isempty(framerate)
    framerate=1;
end
if nargin<4 || isempty(bodylength)
    bodylength=1;
end
if nargin<2 || isempty(edges)
    edges=0:1:20;
end
if nargin<4 || isempty(bodylength)
    bodylength=1;
end
if nargin<5 || isempty(normalization)
    normalization='none';
end
if nargin<6 || isempty(method)
    method='hist';
end

if nargin<7 || isempty(kds_bandwidth)
    kds_bandwidth=1;
end

[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

% Rotate
[~,~,focal_a_rotated]=idSocial_auxiliaries_rotateToFocalSystem(tr,vel,acc);

no_bins=size(edges,2)-1;
speeddistribution=NaN(no_fish,no_fish,no_bins);
speed_median=NaN(no_fish,no_fish);
speed_mode=NaN(no_fish,no_fish);
inbins=cell(no_fish,no_fish,no_bins);
speed = NaN(no_fish,no_fish,no_frames);
if strcmpi(method,'histabs')
    focal_a_rotated = abs(focal_a_rotated);
end
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf
            val=squeeze(focal_a_rotated(:,ff,nf,1));
            speed(ff,nf,:)=val;
            if strcmpi(method,'hist') || strcmpi(method,'histabs')
                [hitemp,bins]=histc(val,edges);
                switch normalization
                    case 'density'
                        speeddistribution(ff,nf,:)=hitemp(1:end-1)/no_frames/(edges(2)-edges(1));
                    case 'no_frames'
                        speeddistribution(ff,nf,:)=hitemp(1:end-1)/no_frames;
                    case 'none'
                        speeddistribution(ff,nf,:)=hitemp(1:end-1);
                end
                if ~all(bins==0)
                    inbins(ff,nf,:)=accumarray(bins(~isnan(val) & bins>0),ones(sum(~isnan(val)& bins>0),1),[no_bins 1],@(x) {x});
                end
                
            
            elseif strcmpi(method,'ksdensity_gauss') || strcmpi(method,'ksdensity_epanechnikov')...
                    || strcmpi(method,'ksdensity_triangular')
                
                edges_kds=edges(1:end-1);%(1:end-1)+(edges(2)-edges(1))/2;
                [hitemp,prob] = idSocial_auxiliaries_kerneldensity(val,edges_kds,kds_bandwidth,method);
                speeddistribution(ff,nf,:)=hitemp;
                for k=1:no_bins
                    inbins{ff,nf,k}=prob(~isnan(prob(:,k)) & prob(:,k)>0,k);
                end
                
            end
            speed_median(ff,nf)=nanmedian(val);

            [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
%             [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
            if ~isempty(locs)
                speed_mode(ff,nf)=edges(locs);
            end                
            
            
        end
    end
end

speeddistributionCell=mat2cell(speeddistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));


%%
fstack=dbstack(3);
funcinfo.Function = mfilename;
if ~isempty(fstack); funcinfo.callerFunction = fstack(1).name; end
funcinfo.no_fish = no_fish;
funcinfo.XTick = edges(1:end-1);
funcinfo.XTickLabel = cellfun(@(x) strtrim(x),cellstr((num2str(edges(1:end-1)')))','UniformOutput',false);
