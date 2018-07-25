function [accdistribution, acc_median,acc_mode, acc_ratio,inbins, acc,funcinfo]=idSocial_accelerationReactionDistribution(tr,edges,angle_range,framerate,bodylength,normalization,method,kds_bandwidth,collapse_left_right,xIsDistance,arena_center)
% Calculates acceleration acceleration (component perpendicular to focal's
% direction of movement)

collapse = true;

if nargin < 3 || isempty(angle_range)
    angle_range = 'all';
end

if nargin < 2 || isempty(edges)
    edges = 0:1:20;
end

if nargin<4 || isempty(framerate)
    framerate=1;
end
if nargin<5 || isempty(bodylength)
    bodylength=1;
end
if nargin<6 || isempty(normalization)
    normalization='density';
end
if nargin<7 || isempty(method)
    method='hist';
end

if nargin<8 || isempty(kds_bandwidth)
    kds_bandwidth=1;
end

if nargin<9 || isempty(collapse_left_right)
    collapse_left_right=false;
end

if nargin<10 || isempty(xIsDistance)
    xIsDistance = 'acceleration';
end

invertY = true;
[tr,vel,acc2,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

%% Prepare

if strcmpi(xIsDistance,'distance_from_center') % Change neighbor coordinates to center of arena, also change acceleration
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
end

[focal_to_neighbour_vector_rotated,~,focal_a_rotated]=idSocial_auxiliaries_rotateToFocalSystem(tr,vel,acc2);
vel_magn=sqrt(sum(vel.^2,4));
vel_norm=bsxfun(@rdivide,vel,vel_magn);


if strcmpi(angle_range,'lateral') || strcmpi(angle_range,'frontal')
    %     Angle between focal movement direction and vector pointing to neighbor position
    foc_to_nb_dir=foc_to_nb_vec./repmat(sqrt(sum(foc_to_nb_vec.^2,4)),[1,1,1,3]);
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
end


%%
dim=2;
no_bins=size(edges,2)-1;
accdistribution=NaN(no_fish,no_fish,no_bins);
acc_median=NaN(no_fish,no_fish);
acc_circvar=NaN(no_fish,no_fish);
acc_ratio=NaN(no_fish,no_fish,no_frames);
acc=NaN(no_fish,no_fish,no_frames);
acc_mode=NaN(no_fish,no_fish);
inbins=cell(no_fish,no_fish,no_bins);

for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            val=focal_a_rotated(:,ff,ff,dim)*framerate^2/bodylength;
            
            if strcmpi(xIsDistance,'distance') || strcmpi(xIsDistance,'distance_to_center')
                distance_vector = focal_to_neighbour_vector_rotated(:,nf,ff,:)/bodylength;
                sign_distance = squeeze(distance_vector(:,1,1,dim))<0;
                distance_vector = sqrt(sum(distance_vector.^2,4));
                distance_vector(sign_distance) = distance_vector(sign_distance)*-1;
            elseif strcmpi(xIsDistance,'left_right')
                distance_vector = focal_to_neighbour_vector_rotated(:,nf,ff,dim)/bodylength;
            elseif strcmpi(xIsDistance,'acceleration')
                distance_vector = val;
            end
            
            if strcmpi(xIsDistance,'distance')
                if strcmpi(angle_range,'all')
                elseif strcmpi(angle_range,'lateral')
                    distance_vector(~(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > pi/4 & ...
                        angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < 3/4*pi & ...
                        focal_to_neighbour_vector_rotated(:,nf,ff,dim) > 0)) = NaN;
                    
                    distance_vector(~(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > pi/4 & ...
                        angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < 3/4*pi & ...
                        focal_to_neighbour_vector_rotated(:,nf,ff,dim) < 0)) = NaN;
                    
                elseif strcmpi(angle_range,'frontal')
                    distance_vector(~(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < pi/4 & ...
                        focal_to_neighbour_vector_rotated(:,nf,ff,dim) > 0)) = NaN;
                    distance_vector(~(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > 3/4*pi& ...
                        focal_to_neighbour_vector_rotated(:,nf,ff,dim) > 0)) = NaN;
                    
                    distance_vector(~(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < pi/4 & ...
                        focal_to_neighbour_vector_rotated(:,nf,ff,dim) < 0)) = NaN;
                    distance_vector(~(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > 3/4*pi& ...
                        focal_to_neighbour_vector_rotated(:,nf,ff,dim) < 0)) = NaN;
                end
                
                ttt = val;
                distance_vector = distance_vector(~isnan(ttt)&ttt~=0);
                ttt = ttt(~isnan(ttt)&ttt~=0);
                if collapse_left_right
                    ttt(distance_vector<0) = -ttt(distance_vector<0);
                    distance_vector = abs(distance_vector);
                    edges(edges<0) = [];
                    
                end
                acc(ff,nf,1:numel(ttt))= ttt;
            else 
                ttt = val.*sign(distance_vector);
                ttt = ttt(~isnan(ttt));
                acc(ff,nf,1:numel(ttt))= ttt;
            end
            
            
            
            
            if strcmpi(method,'hist')
                try
                    [hitemp,bins]=histc(ttt,edges);
                catch
                    keyboard
                end
                
                switch normalization
                    case 'density'
                        accdistribution(ff,nf,:)=hitemp(1:end-1)/no_frames/(edges(2)-edges(1));
                    case 'no_frames'
                        accdistribution(ff,nf,:)=hitemp(1:end-1)/no_frames;
                    case 'none'
                        accdistribution(ff,nf,:)=hitemp(1:end-1);
                end
                if ~all(bins==0)
                    try
                        inbins(ff,nf,:)=accumarray(bins(~isnan(ttt) & bins>0 & bins <= no_bins),ttt(~isnan(ttt)& bins>0 & bins <= no_bins),[no_bins 1],@(x) {x});
                    catch
                        keyboard
                    end
                end
            elseif strcmpi(method,'ksdensity_gauss') || strcmpi(method,'ksdensity_epanechnikov')...
                    || strcmpi(method,'ksdensity_triangular')
                
                edges_kds=edges(1:end-1);%(1:end-1)+(edges(2)-edges(1))/2;
                [hitemp,prob] = idSocial_auxiliaries_kerneldensity(ttt,edges_kds,kds_bandwidth,method);
                accdistribution(ff,nf,:)=hitemp;
                for k=1:no_bins
                    inbins{ff,nf,k}=prob(~isnan(prob(:,k)) & prob(:,k)>0,k);
                end
                
            end
            acc_median(ff,nf)=nanmedian(ttt);
            
            acc_circvar(ff,nf)=idSocial_circ_var(ttt);
            acc_ratio(ff,nf,1:numel(ttt))=ttt;
            
            [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
            if ~isempty(locs)
                acc_mode(ff,nf)=edges(locs);
            end
            
            edgesOut=edges(1:end-1);
        end
    end
end


accdistribution=mat2cell(accdistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));

%%
fstack=dbstack(3);
funcinfo.Function = mfilename;
if ~isempty(fstack)
    if ~isempty(fstack); funcinfo.callerFunction = fstack(1).name; end
end
funcinfo.no_fish = no_fish;
funcinfo.XTick = edgesOut;
funcinfo.XTickLabel = cellfun(@(x) strtrim(x),cellstr((num2str(edges(1:end-1)')))','UniformOutput',false);

if strcmpi(xIsDistance,'distance_from_center') && no_fishOrig==1
    accdistribution = accdistribution(1,2,:);
    inbins = inbins(1,2,:);
    acc_median = acc_median(1,2,:);
    acc_mode = acc_mode(1,2,:);
    acc_ratio = acc_ratio(1,2,:);
    acc = acc(1,2,:);
end