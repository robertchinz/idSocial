function [distancedistribution, distance_median,distance_mode,inbins,distance,nb_rank, edges]=idSocial_accelerationOrientationDistribution(trajectory,edges,normalization,method,kds_bandwidth,collapse,periodic)
% Calculates mutual distances between group members



if nargin<3 || isempty(normalization)
    normalization='density';
end
if nargin<4 || isempty(method)
    method='hist';
end

if nargin<5 || isempty(kds_bandwidth)
    kds_bandwidth=.5;
end
if nargin<6 || isempty(collapse)
    collapse = true;
end
if nargin<7 || isempty(periodic)
    periodic = false;
end

if periodic && collapse
    warning([mfilename ': Options ''periodic'' and ''collapse'' are not compatible. ''periodic'' has been set to ''false''.'])
end

invertY = true;
[trajectory,~,~,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(trajectory,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(trajectory);

%% Prepare

vel=NaN(size(trajectory));
vel(1:no_frames-1,:,:,:)=diff(trajectory,1,1);
vel_magn=sqrt(sum(vel.^2,4));
vel_norm=bsxfun(@rdivide,vel,vel_magn);

acc=NaN(size(trajectory));
acc(1:no_frames-2,:,:,:)=diff(trajectory,2,1);
acc_magn=sqrt(sum(acc.^2,4));
acc_norm=bsxfun(@rdivide,acc,acc_magn);

% Calculations
% Focal to neighbor direction:
foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,3);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,nf,ff,:)=squeeze(trajectory(:,nf,ff,:)-trajectory(:,ff,ff,:));
        end
    end
end
% distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,4))./bodylength;
foc_to_nb_dir=foc_to_nb_vec./repmat(sqrt(sum(foc_to_nb_vec.^2,4)),[1,1,1,3]);
% Angle between focal movement direction and vector pointing to neighbor position
ang_foc_vel_norm_z_plane_and_nb_pos=NaN(no_frames,no_fish,no_fish);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            if ~collapse
                 ang_foc_vel_norm_z_plane_and_nb_pos(:,ff,nf)=...
                    mod(atan2(foc_to_nb_dir(:,nf,ff,1).*acc_norm(:,ff,ff,2)-acc_norm(:,ff,ff,1).*foc_to_nb_dir(:,nf,ff,2),acc_norm(:,ff,ff,1).*foc_to_nb_dir(:,nf,ff,1)+acc_norm(:,ff,ff,2).*foc_to_nb_dir(:,nf,ff,2)),2*pi);
%   
            else
                ang_foc_vel_norm_z_plane_and_nb_pos(:,ff,nf)=...
                    acos(sum(foc_to_nb_dir(:,nf,ff,:).*vel_norm(:,ff,ff,:),4));
            end
        end
    end
end
angle_focal_movement_direction_and_neighbour_position=ang_foc_vel_norm_z_plane_and_nb_pos;

%%
no_bins=size(edges,2)-1;
distancedistribution=NaN(no_fish,no_fish,no_bins);
distance=NaN(no_fish,no_fish,no_frames);
distance_median=NaN(no_fish,no_fish);
inbins=cell(no_fish,no_fish,no_bins);
distance_mode=NaN(no_fish,no_fish);



for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            
            val=squeeze(angle_focal_movement_direction_and_neighbour_position(:,ff,nf));
            
            if periodic
                de = edges(2)-edges(1);
                valhi = vertcat(val, val(val>=edges(end)-de)-edges(end), val(val<edges(1)+de)+edges(end));
            end
            if strcmpi(method,'hist')
                %                 tic
                
                [hitemp,bins]=histc(valhi,edges);
                
                %                 toc
                
                switch normalization
                    case 'density'
                        distancedistribution(ff,nf,:)=hitemp(1:end-1)/no_frames/(edges(2)-edges(1));
                    case 'no_frames'
                        distancedistribution(ff,nf,:)=hitemp(1:end-1)/no_frames;
                    case 'none'
                        distancedistribution(ff,nf,:)=hitemp(1:end-1);
                end
                if ~all(bins==0)
                    inbins(ff,nf,:)=accumarray(bins(~isnan(valhi) & bins>0),ones(sum(~isnan(valhi)& bins>0 & bins <= no_bins),1),[no_bins 1],@(x) {x});
                end
            elseif strcmpi(method,'ksdensity_gauss') || strcmpi(method,'ksdensity_epanechnikov')...
                        || strcmpi(method,'ksdensity_triangular')
               
                edges_kds=edges(1:end-1);%(1:end-1)+(edges(2)-edges(1))/2;                
                [hitemp,prob] = idSocial_auxiliaries_kerneldensity(valhi,edges_kds,kds_bandwidth,method);
                distancedistribution(ff,nf,:)=hitemp;
               
                for k=1:no_bins
                    inbins{ff,nf,k}=prob(~isnan(prob(:,k)) & prob(:,k)>0,k);
                end
                
            end
            distance_median(ff,nf)=nanmedian(val);
            distance(ff,nf,:)=val;

%             [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
            [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
            if ~isempty(locs)
                distance_mode(ff,nf)=edges(locs);
            end
            
           
            
        end
    end
end
nb_rank = NaN(size(distance));
for ff=1:no_fish
    if ~ any(rand_check(ff,:),2)
        
        [~, sort_idces_temp] = sort(distance(ff,:,:));
        sort_idces_temp(sort_idces_temp==ff ) = NaN;
        sort_idces_temp(:,:,all(isnan(distance(ff,:,:)),2)) = NaN;
        nb_rank(ff,1:ff-1,:) = sort_idces_temp(1,1:ff-1,:);
        nb_rank(ff,ff+1:end,:) = sort_idces_temp(1,ff:end-1,:);
    end
end

distancedistribution=mat2cell(distancedistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));