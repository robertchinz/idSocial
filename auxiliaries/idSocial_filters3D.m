function [tr_filter,allFilters] = idSocial_filters3D(tr,options,bodylength,framerate,roi,CustomFrames)

if nargin<3 || isempty(bodylength)
   bodylength = 1;
end
if nargin<4 || isempty(framerate)
   framerate = 1;
end

if nargin<5 || isempty(roi)
   roi = [];
   croi = [];
else
%     if size(roi,1)>1 % Rectangular arena
%     else % Circular Arena
        croi = roi;
        if size(croi,1)>1 && size(croi,2)==1
            croi = croi';
        end
%     end
end

if nargin < 6 
    CustomFrames = -1;
end


minx = -inf;
maxx = inf;
miny = -inf;
maxy = inf;
minr = -inf;
maxr = inf;
center = [0,0,0];
def_options.filter_distancelimits_bl = [-inf inf];
def_options.filter_neighbor_speedlimits_bl_per_s = [-inf inf];
def_options.filter_focal_speedlimits_bl_per_s = [-inf inf];
def_options.filter_neighbor_accelerationlimits_bl_per_s2 = [-inf inf];
def_options.filter_focal_accelerationlimits_bl_per_s2  = [-inf inf];
def_options.filter_focal_rectangularROI = [minx maxx miny maxy];
def_options.filter_neighbor_rectangularROI = [minx maxx miny maxy];
def_options.filter_focal_circularROI = [center minr maxr];
def_options.filter_neighbor_circularROI = [center minr maxr];
def_options.filter_focal_spatial_sectors = [];
def_options.filter_framesWithAllNeighborsOnly = false;
def_options.filter_focal_list = [];
def_options.filter_neighbor_list = [];
% def_options.filter_CustomFrames = [];


[~, options] = idSocial_readparams([],options,def_options);
optionfieldnames = fieldnames(options);
for fn = 1:size(optionfieldnames)
    if isempty(options.(optionfieldnames{fn}))
        options.(optionfieldnames{fn}) = def_options.(optionfieldnames{fn});
    end
end

filter_focal_list = options.filter_focal_list;

filter_neighbor_list = options.filter_neighbor_list;

filter_distancelimits_bl = options.filter_distancelimits_bl;

filter_neighbor_speedlimits_bl_per_s = options.filter_neighbor_speedlimits_bl_per_s;

filter_focal_speedlimits_bl_per_s = options.filter_focal_speedlimits_bl_per_s;


filter_neighbor_accelerationlimits_bl_per_s2 = options.filter_neighbor_accelerationlimits_bl_per_s2;

filter_focal_accelerationlimits_bl_per_s2 = options.filter_focal_accelerationlimits_bl_per_s2;

filter_focal_rectangularROI = options.filter_focal_rectangularROI;

filter_neighbor_rectangularROI = options.filter_neighbor_rectangularROI;


if ~isempty(croi) && ~any(isnan(croi)) && ~isempty(options.filter_focal_circularROI) && ...
        iscell(options.filter_focal_circularROI) % Only radius given, take center from croi and scale radius according to second entry
    
    if strcmpi(options.filter_focal_circularROI{3},'BL') % 
        filter_focal_circularROI = [croi(1:3) options.filter_focal_circularROI{1}*bodylength options.filter_focal_circularROI{2}*bodylength];
    elseif strcmpi(options.filter_focal_circularROI{3},'Arena')
        filter_focal_circularROI = [croi(1:3) croi(4)*options.filter_focal_circularROI{1} croi(4)*options.filter_focal_circularROI{2}];
    else %  Pixels
         filter_focal_circularROI = [croi(1:3) options.filter_focal_circularROI{1} options.filter_focal_circularROI{2}];
    end
elseif (isempty(croi) || any(isnan(croi)) ) && iscell(options.filter_focal_circularROI)
    filter_focal_circularROI = [center minr maxr];
else %if isnumeric(options.filter_focal_circularROI)
    filter_focal_circularROI = options.filter_focal_circularROI;
    
end


if ~isempty(croi) && ~any(isnan(croi)) && ~isempty(options.filter_neighbor_circularROI) && ...
        iscell(options.filter_neighbor_circularROI) % Only radius given, take center from croi and scale radius according to second entry
    
    if strcmpi(options.filter_neighbor_circularROI{3},'BL') % 
        filter_neighbor_circularROI = [croi(1:3) options.filter_neighbor_circularROI{1}*bodylength options.filter_neighbor_circularROI{2}*bodylength];
    elseif strcmpi(options.filter_neighbor_circularROI{3},'Arena')
        filter_neighbor_circularROI = [croi(1:3) croi(4)*options.filter_neighbor_circularROI{1} croi(4)*options.filter_neighbor_circularROI{2}];
    else %  Pixels
         filter_neighbor_circularROI = [croi(1:3) options.filter_neighbor_circularROI{1} options.filter_neighbor_circularROI{2}];
    end
elseif (isempty(croi) || any(isnan(croi)) ) && iscell(options.filter_neighbor_circularROI)
    filter_neighbor_circularROI = [center minr maxr];
else %if isnumeric(options.filter_neighbor_circularROI)
    filter_neighbor_circularROI = options.filter_neighbor_circularROI;
    
end

% filter_neighbor_circularROI = options.filter_neighbor_circularROI;
% filter_neighbor_circularROI(1:2) = filter_neighbor_circularROI([2 1]);

filter_focal_spatial_sectors = options.filter_focal_spatial_sectors;

% if ndims(tr)==3
%     no_dim = size(tr,3);
%     % Add neighbour dimension if necessary
%     tr2=zeros(no_frames,no_focals,no_neighbors,no_dim);
%     tr2(:,:,:,1:no_dim)=reshape(repmat(tr,[1,no_fish]),no_frames,no_focals,no_neighbors,no_dim);
%     tr=tr2;
%     clear tr2;
% end    
%  no_dim = size(tr,4);
celltr=false;
if iscell(tr)|| isstruct(tr); celltr = true; tr_orig=tr; end;
invertY = false;%true;
% warning([mfilename ': invertY switched off. Please check!']);
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);

no_frames = size(tr,1);
no_focals = size(tr,2);
no_neighbors = size(tr,3);

% Use "rand_check" to filter focal and neighbors also in not-rand case
if ~isempty(filter_focal_list)
    for ff=setxor(filter_focal_list,1:no_focals)
        for nf=setxor(filter_focal_list,1:no_focals)
            tr(:,ff,nf,1)=NaN(1,no_frames);
            tr(:,ff,nf,2)=Inf(1,no_frames);
            
        end
    end
end
if ~isempty(filter_neighbor_list)
    for nf=setxor(filter_neighbor_list,1:no_neighbors)
        for ff=1:no_focals
            tr(:,ff,nf,1)=NaN(1,no_frames);
            tr(:,ff,nf,2)=Inf(1,no_frames);
            
            tr(:,nf,ff,1)=NaN(1,no_frames);
            tr(:,nf,ff,2)=Inf(1,no_frames);
        end
    end
end
    
    
    rand_check = idSocial_auxiliaries_trRandCheck(tr);


% if isempty(filter_focal_list)
% filter_focal_list = 1:no_focals;
% else
%     no_focals = numel(filter_focal_list);
% end
% if isempty(filter_neighbor_list)
% filter_neighbor_list = 1:no_neighbors;
%     else
%     no_neighbors = numel(filter_neighbor_list);
% end

% if invertY && ~isempty(filter_focal_circularROI)
%     filter_focal_circularROI(2) = -filter_focal_circularROI(2);
% end
% if invertY && ~isempty(filter_neighbor_circularROI)
%     filter_neighbor_circularROI(2) = -filter_neighbor_circularROI(2);
% end
%% Spatial sectors
if ~isempty(filter_focal_spatial_sectors) && (strcmpi(filter_focal_spatial_sectors,'lateral') || strcmpi(filter_focal_spatial_sectors,'frontal') || strcmpi(filter_focal_spatial_sectors,'behind'))
    
    collapse = true;
    foc_to_nb_vec=NaN(no_frames,no_focals,no_neighbors,3);
    ang_foc_vel_norm_and_nb_vel_norm=NaN(no_frames,no_focals,no_neighbors);
  
    vel_magn=sqrt(sum(vel.^2,4));
    vel_norm=bsxfun(@rdivide,vel,vel_magn);
    for ff=1:no_focals
        for nf=1:no_neighbors
            if ff~=nf && ~rand_check(ff,nf)
                foc_to_nb_vec(:,nf,ff,:)=squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
                ang_foc_vel_norm_and_nb_vel_norm(:,ff,nf)=...
                    mod(atan2(vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,2)-vel_norm(:,ff,ff,1).*vel_norm(:,nf,ff,2),vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,1)+vel_norm(:,nf,ff,2).*vel_norm(:,ff,ff,2)),2*pi);
            end
        end
    end
    
    % Angle between focal movement direction and vector pointing to neighbor position
    foc_to_nb_dir=foc_to_nb_vec./repmat(sqrt(sum(foc_to_nb_vec.^2,4)),[1,1,1,3]);
    ang_foc_vel_norm_z_plane_and_nb_pos=NaN(no_frames,no_focals,no_neighbors);
    for ff=1:no_focals
        for nf=1:no_neighbors
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
%     focal_to_neighbour_vector_rotated=idSocial_auxiliaries_rotateToFocalSystem(tr,vel,acc);

    
   
    if strcmpi(filter_focal_spatial_sectors,'lateral')
        
        for ff=1:no_focals
            for nf=1:no_neighbors
                if ff~=nf && ~rand_check(ff,nf)
                    distance_vector = false(no_frames,1);
                    distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > pi/4 & ...
                        angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < 3/4*pi) = true;
                    tr(~distance_vector,nf,ff,:) = NaN;
                    
%                     focal_to_neighbour_vector_rotated=idSocial_auxiliaries_rotateToFocalSystem(tr,vel,acc);
%                     figure; plot(focal_to_neighbour_vector_rotated(:,nf,ff,1),focal_to_neighbour_vector_rotated(:,nf,ff,2),'.')
                end
            end
        end
    elseif strcmpi(filter_focal_spatial_sectors,'frontal')
         for ff=1:no_focals
            for nf=1:no_neighbors
                if ff~=nf && ~rand_check(ff,nf)
                    distance_vector = false(no_frames,1);
                    distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < pi/4 ) = true;
                    distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > 3/4*pi) = true;
                    tr(~distance_vector,nf,ff,:) = NaN;
                end
            end
        end
    elseif strcmpi(filter_focal_spatial_sectors,'behind')
         for ff=1:no_focals
            for nf=1:no_neighbors
                if ff~=nf && ~rand_check(ff,nf)
                    distance_vector = false(no_frames,1);
                    distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) > pi/2 ) = true;
%                     distance_vector(angle_focal_movement_direction_and_neighbour_position(:,ff,nf) < 3/2*pi) = true;
                    tr(~distance_vector,nf,ff,:) = NaN;
                end
            end
        end
        
    end
end
 
%% Calculate some velocity etc. necessary for filtering

if iscell(filter_focal_speedlimits_bl_per_s) && size(filter_focal_speedlimits_bl_per_s,1)==2
    filter_focal_speedlimits_bl_per_sTEMP = NaN(1,2);
    for k=1:2
        if strcmpi(filter_focal_speedlimits_bl_per_s{k,1},'prctile')
            filter_focal_speedlimits_bl_per_sTEMP(k) = prctile(vel_magn(:),filter_focal_speedlimits_bl_per_s{k,2});
        end
    end
    filter_focal_speedlimits_bl_per_s = filter_focal_speedlimits_bl_per_sTEMP;
elseif iscell(filter_focal_speedlimits_bl_per_s) && size(filter_focal_speedlimits_bl_per_s,1)~=2
    error([mfilename ': Unknown format for filter_focal_speedlimits_bl_per_s']);    
end

% vel_magn=sqrt(sum(vel.^2,4));
% acc_magn=sqrt(sum(acc.^2,4));

foc_to_nb_vec=NaN(no_frames,no_dim);
% foc_to_center_vec=NaN(no_frames,no_dim);
% nb_to_center_vec=NaN(no_frames,no_dim);
if celltr
          vel_out = vel;%(:,no_focals,no_neighbors,:);
          acc_out = acc;%(:,no_focals,no_neighbors,:);
end
tr_filter = tr;

% tr_filter = tr(:,no_focals,no_neighbors,:);


for ff=1:no_focals
    for nf=1:no_neighbors
        vel_magn=sqrt(sum(squeeze(vel(:,ff,nf,:)).^2,2));
        acc_magn=sqrt(sum(squeeze(acc(:,ff,nf,:)).^2,2));
        
        vel_magnNb=sqrt(sum(squeeze(vel(:,nf,ff,:)).^2,2));
        acc_magnNb=sqrt(sum(squeeze(acc(:,nf,ff,:)).^2,2));
        if ff~=nf
            foc_to_nb_vec =squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
            
        end
        foc_to_center_vec =repmat(filter_focal_circularROI(1:no_dim),[no_frames 1 1 1])-squeeze(tr(:,ff,ff,:));
        nb_to_center_vec =repmat(filter_neighbor_circularROI(1:no_dim),[no_frames 1 1 1])-squeeze(tr(:,nf,ff,:));
        
        distance_focal_center=sqrt(sum(foc_to_center_vec.^2,2));%./bodylength;
        distance_neighbor_center=sqrt(sum(nb_to_center_vec.^2,2));%./bodylength;
        distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,2))./bodylength;
        
        % Apply filters for focal - neighbor pair
        if ff==nf
            vel_filter = vel_magn < filter_focal_speedlimits_bl_per_s(1) | ...
                vel_magn > filter_focal_speedlimits_bl_per_s(2);
            
            
            acc_filter = acc_magn < filter_focal_accelerationlimits_bl_per_s2(1) | ...
                acc_magn > filter_focal_accelerationlimits_bl_per_s2(2);
            
            
            circRoi = (distance_focal_center < filter_focal_circularROI(4)) | (distance_focal_center > filter_focal_circularROI(5));
            rectRoi = (tr(:,ff,nf,1) < filter_focal_rectangularROI(1,1)) | ...
                (tr(:,ff,nf,2) < filter_focal_rectangularROI(1,3)) | ...
                (tr(:,ff,nf,1) > filter_focal_rectangularROI(1,2)) | ...
                (tr(:,ff,nf,2) > filter_focal_rectangularROI(1,4));
            %             else
            %                 circRoi = false(no_frames,1);
            %                 rectRoi = false(no_frames,1);
            %             end
            
            try
            allFilters = vel_filter | acc_filter | rectRoi |  circRoi ;
            catch 
                keyboard
            end
            
        elseif nf ~= ff
            dist_filter = distance_focal_neighbour < filter_distancelimits_bl(1) | ...
                distance_focal_neighbour > filter_distancelimits_bl(2);
            velnb_filter = vel_magnNb < filter_neighbor_speedlimits_bl_per_s(1) | ...
                vel_magnNb > filter_neighbor_speedlimits_bl_per_s(2);
            accnb_filter = acc_magnNb < filter_neighbor_accelerationlimits_bl_per_s2(1) | ...
                acc_magnNb > filter_neighbor_accelerationlimits_bl_per_s2(2);
            rectRoiNb = (tr(:,nf,ff,1) < filter_neighbor_rectangularROI(1,1)) | ...
                (tr(:,nf,ff,2) < filter_neighbor_rectangularROI(1,3)) | ...
                (tr(:,nf,ff,1) > filter_neighbor_rectangularROI(1,2)) | ...
                (tr(:,nf,ff,2) > filter_neighbor_rectangularROI(1,4));
            circRoiNb = (distance_neighbor_center < filter_neighbor_circularROI(4)) | (distance_neighbor_center > filter_neighbor_circularROI(5));
            %             allFilters = allFilters | velnb_filter | accnb_filter | dist_filter |rectRoiNb |  circRoiNb;
            allFilters =  velnb_filter | accnb_filter | dist_filter |rectRoiNb |  circRoiNb;

        end
        
        tr_filter(allFilters,ff,nf,:) = ...
            NaN(sum(allFilters),no_dim);
        
        if celltr
            vel_out(allFilters,ff,nf,:) = ...
                NaN(sum(allFilters),no_dim);
            acc_out(allFilters,ff,nf,:) = ...
                NaN(sum(allFilters),no_dim);
        end
        
        % End apply filter
        
        
    end
end
% distance_focal_center=sqrt(sum(foc_to_center_vec.^2,4));%./bodylength;
% distance_neighbor_center=sqrt(sum(nb_to_center_vec.^2,4));%./bodylength;
% distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,4))./bodylength;



%% Apply filters
% tr_filter = tr;
% if celltr
%           vel_out = vel;
%           acc_out = acc;
% end
% for ff = 1:no_focals
%     for nf = 1:no_neighbors

%         if ff==nf
%             vel_filter = vel_magn(:,ff,nf) < filter_focal_speedlimits_bl_per_s(1) | ...
%                 vel_magn(:,ff,nf) > filter_focal_speedlimits_bl_per_s(2);
%             
%             
%             acc_filter = acc_magn(:,ff,nf) < filter_focal_accelerationlimits_bl_per_s2(1) | ...
%                 acc_magn(:,ff,nf) > filter_focal_accelerationlimits_bl_per_s2(2);
%             
%             
%             circRoi = (distance_focal_center(:,ff,nf) < filter_focal_circularROI(4)) | (distance_focal_center(:,ff,nf) > filter_focal_circularROI(5));
%             rectRoi = (tr(:,ff,nf,1) < filter_focal_rectangularROI(1,1)) | ...
%                 (tr(:,ff,nf,2) < filter_focal_rectangularROI(1,3)) | ...
%                 (tr(:,ff,nf,1) > filter_focal_rectangularROI(1,2)) | ...
%                 (tr(:,ff,nf,2) > filter_focal_rectangularROI(1,4));
%             %             else
%             %                 circRoi = false(no_frames,1);
%             %                 rectRoi = false(no_frames,1);
%             %             end
%             
%             
%             allFilters = vel_filter | acc_filter | rectRoi |  circRoi ;
%             
%         elseif nf ~= ff
%             dist_filter = distance_focal_neighbour(:,nf,ff) < filter_distancelimits_bl(1) | ...
%                 distance_focal_neighbour(:,nf,ff) > filter_distancelimits_bl(2);
%             velnb_filter = vel_magn(:,nf,ff) < filter_neighbor_speedlimits_bl_per_s(1) | ...
%                 vel_magn(:,nf,ff) > filter_neighbor_speedlimits_bl_per_s(2);
%             accnb_filter = acc_magn(:,nf,ff) < filter_neighbor_accelerationlimits_bl_per_s2(1) | ...
%                 acc_magn(:,nf,ff) > filter_neighbor_accelerationlimits_bl_per_s2(2);
%             rectRoiNb = (tr(:,nf,ff,1) < filter_neighbor_rectangularROI(1,1)) | ...
%                 (tr(:,nf,ff,2) < filter_neighbor_rectangularROI(1,3)) | ...
%                 (tr(:,nf,ff,1) > filter_neighbor_rectangularROI(1,2)) | ...
%                 (tr(:,nf,ff,2) > filter_neighbor_rectangularROI(1,4));
%             circRoiNb = (distance_neighbor_center(:,ff,nf) < filter_neighbor_circularROI(4)) | (distance_neighbor_center(:,ff,nf) > filter_neighbor_circularROI(5));
%             %             allFilters = allFilters | velnb_filter | accnb_filter | dist_filter |rectRoiNb |  circRoiNb;
%             allFilters =  velnb_filter | accnb_filter | dist_filter |rectRoiNb |  circRoiNb;
% 
%         end
%         
%         tr_filter(allFilters,ff,nf,:) = ...
%             NaN(sum(allFilters),no_dim);
%         
%         if celltr
%             vel_out(allFilters,ff,nf,:) = ...
%                     NaN(sum(allFilters),no_dim);
%                 acc_out(allFilters,ff,nf,:) = ...
%                     NaN(sum(allFilters),no_dim);
%         end
%         
%     end
% end
% Apply CustomFrames filter
try
if ~isempty(CustomFrames) && ~(isnumeric(CustomFrames)&& numel(CustomFrames) == 1 && CustomFrames==-1)
    if (isnumeric(CustomFrames) || islogical(CustomFrames)) && ndims(CustomFrames)==3 && ...
            isequal(size(CustomFrames),[no_frames no_fish no_fish]) 
        tempFilter = repmat(CustomFrames,[1 1 1 no_dim]);
        tr_filter(~logical(tempFilter)) = NaN;
    elseif (isnumeric(CustomFrames) || islogical(CustomFrames)) && ndims(CustomFrames)==4 && ...
            isequal(size(CustomFrames),[no_frames no_fish no_fish no_dim]) 
        tr_filter(~logical(CustomFrames)) = NaN;
    end
elseif CustomFrames>-1
    tr_filter = NaN(no_frames, no_fish, no_fish, no_dim);
end
catch
    keyboard
end
if options.filter_framesWithAllNeighborsOnly
    for ff=1:no_fish
        allNbPresent = all(~any(isnan(tr_filter(:,:,ff,:)),4),2);
        tr_filter(~allNbPresent,:,ff,:) = ...
            NaN(sum(~allNbPresent),no_fish,1,no_dim);
        if celltr
            vel_out(~allNbPresent,nf,ff,:) = ...
                NaN(sum(~allNbPresent),no_fish,1,no_dim);
            acc_out(~allNbPresent,nf,ff,:) = ...
                NaN(sum(~allNbPresent),no_fish,1,no_dim);
        end
    end
end
if celltr
%     for k=1:numel(tr_orig)
          tr_out.Tr = tr_filter;
          tr_out.Vel = vel_out;
          tr_out.Acc = acc_out;
          tr_filter = tr_out;
%     end
end
%%
% figure; plot(tr_filter(:,1,2,1),tr_filter(:,1,2,2),'.'); axis equal; hold on; plot(tr_filter(:,1,1,1),tr_filter(:,1,1,2),'r.')
% hold on; plot(roi(1)+roi(4)*cos(0:.01:2*pi+.01),roi(2)+roi(4)*sin(0:.01:2*pi+.01),'k')