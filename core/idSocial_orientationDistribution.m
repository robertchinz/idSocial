function [distancedistribution, distance_median,distance_mode,inbins,distance]=idSocial_orientationDistribution(trajectory,edges,normalization,method,statistical_operation,kds_bandwidth,collapse)
% Calculates mutual distances between group members



if nargin<3 || isempty(normalization)
    normalization='density';
end
if nargin<4 || isempty(method)
    method='hist';
end
if nargin<5 || isempty(statistical_operation)
    statistical_operation='median';
end
if nargin<6 || isempty(kds_bandwidth)
    kds_bandwidth=.5;
end
if nargin<7 || isempty(collapse)
    collapse = true;
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
                    mod(atan2(foc_to_nb_dir(:,nf,ff,1).*vel_norm(:,ff,ff,2)-vel_norm(:,ff,ff,1).*foc_to_nb_dir(:,nf,ff,2),vel_norm(:,ff,ff,1).*foc_to_nb_dir(:,nf,ff,1)+vel_norm(:,ff,ff,2).*foc_to_nb_dir(:,nf,ff,2)),2*pi);
%   
            else
                ang_foc_vel_norm_z_plane_and_nb_pos(:,ff,nf)=...
                    acos(sum(foc_to_nb_dir(:,nf,ff,:).*vel_norm(:,ff,ff,:),4));
            end
        end
    end
end
angle_focal_movement_direction_and_neighbour_position=ang_foc_vel_norm_z_plane_and_nb_pos;

%% Debug
% idx = 15420:35520;
% sidx = idx(round(numel(idx)/2));
% ff=1; nf=2;
% sc=200;
% figure; 
% hold on
% plot(squeeze(trajectory(idx,ff,ff,1)),squeeze(trajectory(idx,ff,ff,2)))
% plot(squeeze(trajectory(idx,nf,ff,1)),squeeze(trajectory(idx,nf,ff,2)))
% axis equal
% % figure
% plot([squeeze(trajectory(sidx,ff,ff,1)) squeeze(trajectory(sidx,ff,ff,1))+squeeze(foc_to_nb_dir(sidx,nf,ff,1))*sc], ...
%     [squeeze(trajectory(sidx,ff,ff,2)) squeeze(trajectory(sidx,ff,ff,2))+squeeze(foc_to_nb_dir(sidx,nf,ff,2))*sc],'g')
% 
% plot([squeeze(trajectory(sidx,ff,ff,1)) squeeze(trajectory(sidx,ff,ff,1))+squeeze(vel_norm(sidx,ff,ff,1))*sc], ...
%     [squeeze(trajectory(sidx,ff,ff,2)) squeeze(trajectory(sidx,ff,ff,2))+squeeze(vel_norm(sidx,ff,ff,2))*sc],'g')
% 
% plot(squeeze(trajectory(sidx,nf,ff,1)),squeeze(trajectory(sidx,nf,ff,2)),'ro')
% title(['Angle ' num2str(angle_focal_movement_direction_and_neighbour_position(sidx,ff,nf)/pi*180)])
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
                
            if strcmpi(method,'hist')
%                 tic
%  keyboard
                [hitemp,bins]=histc(val,edges);
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
                    inbins(ff,nf,:)=accumarray(bins(~isnan(val) & bins>0),ones(sum(~isnan(val)& bins>0 & bins <= no_bins),1),[no_bins 1],@(x) {x});
                end
            elseif strcmpi(method,'ksdensity_gauss') || strcmpi(method,'ksdensity_epanechnikov')...
                        || strcmpi(method,'ksdensity_triangular')
               
                edges_kds=edges(1:end-1);%(1:end-1)+(edges(2)-edges(1))/2;                
                [hitemp,prob] = idSocial_auxiliaries_kerneldensity(val,edges_kds,kds_bandwidth,method);
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

distancedistribution=mat2cell(distancedistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));