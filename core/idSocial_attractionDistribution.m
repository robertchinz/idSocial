function [attractiondistribution, attraction_median,attraction_mode,attraction_circvar, attraction_ratio,inbins, attraction]=idSocial_attractionDistribution(tr,edges,angle_range,bodylength,normalization,method,statistical_operation,kds_bandwidth,angles_2pi,rotate_z)
% Calculates mutual distances between group members
no_fish=size(tr,2);
no_frames=size(tr,1);
no_dim=size(tr,4);


xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation
if nargin < 3 || isempty(angle_range)
    angle_range = 'all';
end

if nargin<4 || isempty(bodylength)
    bodylength=1;
end
if nargin<5 || isempty(normalization)
    normalization='density';
end
if nargin<7 || isempty(method)
    method='hist';
end
if nargin<8 || isempty(statistical_operation)
    statistical_operation='median';
end
if nargin<9 || isempty(kds_bandwidth)
    kds_bandwidth=.5;
end
if nargin<10 || isempty(angles_2pi)
    angles_2pi = false;
end
if nargin<11 || isempty(rotate_z)
    rotate_z=false;
end

invertY = true;
[tr,~,~,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(tr);


%% Prepare

if no_dim==2
    trtemp=zeros(no_frames,no_fish,no_fish,3);
    trtemp(:,:,:,1:no_dim)=tr;
    tr=trtemp;
end
tr(:,:,:,2)=-tr(:,:,:,2); % This is a dirty trick in order to make the orientation of the coordinate system coherent with the videos



vel=NaN(size(tr));
acc=NaN(size(tr));


vel(1:no_frames-1,:,:,:)=diff(tr,1,1);
acc(1:no_frames-2,:,:,:)=diff(tr,2,1);
vel_magn=sqrt(sum(vel.^2,4));
vel_norm=bsxfun(@rdivide,vel,vel_magn);
acc_magn=sqrt(sum(acc.^2,4));
acc_norm=bsxfun(@rdivide,acc,acc_magn);

foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,no_dim);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,nf,ff,:)=squeeze(tr(:,nf,ff,1:no_dim)-tr(:,ff,ff,1:no_dim));
        end
    end
end
foc_to_nb_vec_magn=sqrt(sum(foc_to_nb_vec.^2,4));
foc_to_nb_vec_norm=bsxfun(@rdivide,foc_to_nb_vec,foc_to_nb_vec_magn);

distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,4))./bodylength;

distance_focal_neighbour=permute(distance_focal_neighbour,[2,3,1]);
distance_focal_neighbour(logical(repmat(eye(no_fish,no_fish),[1 1 no_frames])))=NaN;


for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,nf,ff,:)=squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
        end
    end
end

if ~strcmpi(angle_range,'all')
    
    % tr(:,:,:,2)=-tr(:,:,2); % This is a dirty trick in order to make the orientation of the coordinate system coherent with the videos
    ang_foc_vel_norm_and_nb_vel_norm=NaN(no_frames,no_fish,no_fish);
    for ff=1:no_fish
        for nf=1:no_fish
            if ff~=nf && ~rand_check(ff,nf)
                ang_foc_vel_norm_and_nb_vel_norm(:,ff,nf)=...
                    mod(atan2(vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,2)-vel_norm(:,ff,ff,1).*vel_norm(:,nf,ff,2),vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,1)+vel_norm(:,nf,ff,2).*vel_norm(:,ff,ff,2)),2*pi);
            end
        end
    end
    ang_vel_norm_z_plane_and_y_axis=...
        mod(atan2(vel_norm(:,:,:,1)*yax(2)-yax(1)*vel_norm(:,:,:,2),vel_norm(:,:,:,1)*yax(1)+vel_norm(:,:,:,2)*yax(2)),2*pi);
    
    
    ang_vel_norm_z_plane_and_z_axis=...
        acos(vel_norm(:,:,:,1)*zax(1)+vel_norm(:,:,:,2)*zax(2)+vel_norm(:,:,:,3)*zax(3));
    ang_vel_norm_and_xy_plane=NaN(size(ang_vel_norm_z_plane_and_z_axis));
    ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis>=0)=...
        pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis>=0);
    ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis<0)=...
        -pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis<0);
    
    % Angle between focal movement direction and vector pointing to neighbor position
    foc_to_nb_dir=foc_to_nb_vec./repmat(sqrt(sum(foc_to_nb_vec.^2,4)),[1,1,1,3]);
    ang_foc_vel_norm_z_plane_and_nb_pos=NaN(no_frames,no_fish,no_fish);
    for ff=1:no_fish
        for nf=1:no_fish
            if ff~=nf && ~rand_check(ff,nf)
                if ~angles_2pi 
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
    %% Measures which depend on relative position/orientation of foc and nb: Rotate vectors
    a1=NaN(no_frames,no_fish);
    for ff=1:no_fish
        a1(:,ff)=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,ff,ff,:));
    end
    a1=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,:,1,:));
    if rotate_z
        a2=NaN(no_frames,no_fish);
        for ff=1:no_fish
            a2(:,ff)=squeeze(ang_vel_norm_and_xy_plane(:,ff,ff,:));
        end
    else
        a2=zeros(size(a1));
    end
    
    a3=zeros(size(a1));
    turnmatr=NaN(3,3,size(a1,1),no_fish);
    
    turnmatr(1,1,:,:)=cos(a3).*cos(a1);
    turnmatr(1,2,:,:)=-cos(a2).*sin(a1)+cos(a1).*sin(a3).*sin(a2);
    turnmatr(1,3,:,:)=sin(a2).*sin(a1)+cos(a2).*sin(a3).*cos(a1);
    
    turnmatr(2,1,:,:)=cos(a3).*sin(a1);
    turnmatr(2,2,:,:)=sin(a2).*sin(a3).*sin(a1)+cos(a2).*cos(a1);
    turnmatr(2,3,:,:)=-sin(a2).*cos(a1)+cos(a2).*sin(a3).*sin(a1);
    
    turnmatr(3,1,:,:)=-sin(a3);
    turnmatr(3,2,:,:)=sin(a2).*cos(a3);
    turnmatr(3,3,:,:)=cos(a2).*cos(a3);
    
    foc_to_nb_vec_rot=NaN(no_frames,no_fish,no_fish,3);
    acc_rot=NaN(no_frames,no_fish,no_fish,3);
    vel_rot=NaN(no_frames,no_fish,no_fish,3);
    
    vel_trans=NaN(size(vel_rot));
    vel_trans(:,:,:,1)=vel_rot(:,:,:,1);
    
    for ff=1:no_fish
        for nf=1:no_fish
            %         if ff~=nf && ~rand_check(ff,nf)
            % Idx (:,ff,nf,:) means: Corresponding value (vel,acc,..) for focal ff 'seen by' neighbor nf;
            % In case of foc_to_nb_vec, vel and acc, 'seen
            % by' refers to neighbor-depend filtering (e.g.,
            % value at (:,ff,nf,:) is NaN because neighbor moves slower than
            % the lower threshold speed).
            % For rotated values, (:,ff,nf,:) means: vector
            % for focal ff in the system where neighbor velocity
            % vector is parallel to the y-axis
            foc_to_nb_vec_rot(:,nf,ff,:)=...
                [sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
                sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
                sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
            vel_rot(:,ff,nf,:)=...
                [sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
                sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
                sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
            acc_rot(:,ff,nf,:)=...
                [sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
                sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
                sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
            vel_trans(:,ff,nf,2)=vel_rot(:,ff,nf,2)-vel_rot(:,ff,ff,2);
            %         end
        end
    end
    focal_to_neighbour_vector_rotated=foc_to_nb_vec_rot;
end
%%
no_bins=size(edges,2)-1;
attractiondistribution=NaN(no_fish,no_fish,no_bins);
attraction_median=NaN(no_fish,no_fish);
attraction_circvar=NaN(no_fish,no_fish);
attraction_ratio=NaN(no_fish,no_fish);
attraction=NaN(no_fish,no_fish,no_frames);
attraction_mode=NaN(no_fish,no_fish);
inbins=cell(no_fish,no_fish,no_bins);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            distance_vector = NaN(no_frames,1);
            if strcmpi(angle_range,'all')
                distance_vector = ones(no_frames,1);
            elseif strcmpi(angle_range,'lateral')
                distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > pi/4 & ...
                    angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < 3/4*pi) = 1;
                
               
                
            elseif strcmpi(angle_range,'frontal')
                distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < pi/4 ) = 1;
                distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > 3/4*pi) = 1;
                
               
            end
            zone_filter = ~isnan(distance_vector);
            
            if ~angles_2pi
                val=acos(sum(squeeze(acc_norm(:,ff,ff,1:2))'.*squeeze(foc_to_nb_vec_norm(:,nf,ff,1:2))',1)');
                %             squeeze(distance_focal_neighbour(ff,nf,:));
            else
                a=squeeze(acc_norm(:,ff,ff,1:2));
                b=squeeze(foc_to_nb_vec_norm(:,nf,ff,1:2));
                val=...
                    atan2(a(:,1).*b(:,2)-b(:,1).*a(:,2),a(:,1).*b(:,1)+a(:,2).*b(:,2));
            end
            val(~zone_filter) = NaN;
            %             keyboard
            if strcmpi(method,'hist')
                
                [hitemp,bins]=histc(val,edges);
                
%                 keyboard
                switch normalization
                    case 'density'
                        attractiondistribution(ff,nf,:)=hitemp(1:end-1)/no_frames/(edges(2)-edges(1));
                    case 'no_frames'
                        attractiondistribution(ff,nf,:)=hitemp(1:end-1)/no_frames;
                    case 'none'
                        attractiondistribution(ff,nf,:)=hitemp(1:end-1);
                end
                if ~all(bins==0)
                    try
                        inbins(ff,nf,:)=accumarray(bins(~isnan(val) & bins>0 & bins <= no_bins),ones(sum(~isnan(val)& bins>0 & bins <= no_bins),1),[no_bins 1],@(x) {x});
                        %                     inbins(ff,nf,:)=accumarray(bins(~isnan(val) & bins>0 & bins <= no_bins),val(~isnan(val)& bins>0 & bins <= no_bins),[no_bins 1],@(x) {x});
                    catch
                        keyboard
                    end
                end
            elseif strcmpi(method,'ksdensity_gauss') || strcmpi(method,'ksdensity_epanechnikov')...
                        || strcmpi(method,'ksdensity_triangular')
                
                edges_kds=edges;%(1:end-1);%(1:end-1)+(edges(2)-edges(1))/2;                
                [hitemp,prob] = idSocial_auxiliaries_kerneldensity(val,edges_kds,kds_bandwidth,method);
                attractiondistribution(ff,nf,:)=hitemp;
                for k=1:no_bins
                    inbins{ff,nf,k}=prob(~isnan(prob(:,k)) & prob(:,k)>0,k);
                end
                
            end
            attraction_median(ff,nf)=nanmedian(val);
            try
                attraction(ff,nf,:)=val;
            catch
                keyboard
            end
            
            attraction_circvar(ff,nf)=idSocial_circ_var(val);
            attraction_ratio(ff,nf)=nansum(val<pi/2)/sum(~isnan(val));
            
            %             [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
            [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
            
            if ~isempty(locs)
                attraction_mode(ff,nf)=edges(locs);
            end
            
                   
            
        end
    end
end

 attractiondistribution=mat2cell(attractiondistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));