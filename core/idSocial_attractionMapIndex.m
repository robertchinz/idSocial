function varargout=idSocial_attractionMapIndex(trajectory,edges,spacing,framerate,bodylength,discard_percentage_outliers,rotate_z,method,angles_2pi)
%   Two-dimensional probability density of relative neighbor
%   positions in reference system of focal individual
%   (defined by focal velocity vector) and maps of the
%   average focal velocity and acceleration depending on
%   neighbor position.
%
%   The reference system is defined by the position and the
%   velocity vector of the respective focal individual (and
%   by the fixed external reference frame given by axis
%   x=[1 0 0], y=[0 1 0] and z=[0 0 1]):
%   The origin is defined by the coordinates of the focal
%   given by trajectory, the y-axis by the
%   direction of movement of the focal and the x-axis by the
%   direction perpendicular to the movement direction and
%   parallel to the x-y-plane of the fixed external
%   reference frame.
%   In case of a three-dimensional trajectory, the direction
%   of the z-axis is either fixed along the z-axis of the
%   external frame or calculated from the focal's movement
%   direction (see input parameter ROTATE_Z).
%
%   The probability density is given by the number of counts
%   which fall into a certain bin, divided by the total
%   length of the time series (the number of frames), and by
%   the bin width in x- and y-direction.
%
%   [rel_positions, ...
%       vertical_acc_map,...
%       turning_map,...
%       vertical_velocity_map,...
%       horizontal_velocity_map,...
%       ...
%       rel_positions_idx, ...
%       vertical_acc_values,...
%       turning_idx_values,...
%       vertical_velocity_values,...
%       horizontal_velocity_values] = ...
%
%               IDSOCIAL_DYNAMICMAPS(TR)
%
%   INPUTS:
%       TRAJECTORY: #Frames-by-#Individuals-by-#Dimensions
%           or  #Frames-by-#Individuals-by-#Individuals-by-
%           #Dimensions array containing the time series.
%           The second case enables non-symmetric filtering,
%           e.g., when values of the focal individual are
%           filtered depending on position, speed, etc. of
%           each neighbor individual
%           (see idSocial_prepareTrajectories3D).
%           #Dimensions can be 2 or 3.
%
%       EDGES: Following the convention of MATLAB(R)'s
%           hist3, edges is a two-element cell array of
%           numeric vectors with monotonically
%           non-decreasing values. Values will be sorted
%           into a 2-D grid of bins with edges at edges{1}
%           in the first dimension and at edges{2} in the
%           second.
%           Note that values in edges are in units of body
%           length and thus affected by input parameter
%           bodylength.
%           DEFAULT: EDGES={-10:2:10; -10:2:10}
%
%       FRAMERATE,BODYLENGTH: Optional scaling of the
%           results according to the given framerate
%           (normally the framerate at which the video has
%           been recorded) and the body length in pixels of
%           either each individual in a 1-by-#Individuals
%           vector or the same scalar value for all
%           individuals.
%           DEFAULT: FRAMERATE = 1, BODYLENGTH = 1
%
%       DISCARD_PERCENTAGE_OUTLIERS: If defined, absolute
%           values of velocity/acceleration greater than the
%           'discard_percentage_outliers'-th prctile will be
%           regarded as outlier and not used for creating
%           the maps.
%           DEFAULT: DISCARD_PERCENTAGE_OUTLIERS = 0
%
%       ROTATE_Z: (ONLY APPLICABLE FOR 3D-TRAJECTORIES)
%           NOTE: SO FAR THIS FUNCTION DOES NOT HAVE 3D
%           FUNCTIONALITY!
%           Defines if for the calculation of the focal's
%           reference system in 3 dimensions, the direction
%           of the z-axis is fixed along the vertic vertical
%           axis of the experimental setup (rotate_z=false),
%           or if it is calculated from the focal's velocity
%           vector and thus its direction is variable
%           (rotate_z=true).
%           DEFAULT: ROTATE_Z = FALSE
%
%   OUTPUTS:
%       REL_POSITIONS: #Indiviudals-by-#Individuals-by-
%           #BinsX-by-#BinsY, where #BinsX/Y equals the
%           number of values in the respective dimension of
%           edges, minus 1.
%           For each focal individual, rel_positions
%           contains the probability density of finding a
%           neighbor in a certain bin at a position defined
%           by the corresponding elements in edges.
%
%           EXAMPLE: rel_position(1,2,3,2) contains the
%           probability density of finding neighbor 2 at a
%           position falling between edges{1}(3) and
%           edges{1}(4) in the x-direction (left-right of
%           the focal)and between edges{2}(2) and
%           edges{2}(3) in the y-direction (back-front) in
%           the reference system of focal 1.

%       VERTICAL_ACC_MAP: #Bins-by-#Bins array
%       TURNING_MAP: #Bins-by-#Bins array
%       VERTICAL_VELOCITY_MAP: #Bins-by-#Bins array
%       HORIZONTAL_VELOCITY_MAP: #Bins-by-#Bins array
%       REL_POSITIONS_IDCES: #Bins-by-#Bins array
%       VERTICAL_ACC_VALUES: #Bins-by-#Bins array
%       TURNING_VALUES: #Bins-by-#Bins array
%       VERTICAL_VELOCITY_VALUES: #Bins-by-#Bins array
%       HORIZONTAL_VELOCITY_VALUES: #Bins-by-#Bins array
%
%
%   See also idSocial_prepareTrajectories3D,
%   idSocial_accelerationDistribution, idSocial_speed, prctile

%   2014 Robert C. Hinz, Gonzalo G. de Polavieja, Consejo Superior de Investigaciones Científicas

if nargin<1 || isempty(trajectory)
    error([mfilename ': Trajectory missing.'])
end

if nargin<2 || isempty(edges)
    edges={-10:2:10; -10:2:10};
end

if nargin<3 || isempty(spacing)
    spacing=2;
end

if nargin<4 || isempty(framerate)
    framerate=1;
end
if nargin<5 || isempty(bodylength)
    bodylength=1;
end
if nargin<6 || isempty(discard_percentage_outliers)
    discard_percentage_outliers=0;
end

if nargin<7 || isempty(rotate_z)
    rotate_z=false;
end

if nargin<8 || isempty(method)
    method='gauss';
end
if nargin<9 || isempty(angles_2pi)
    angles_2pi=0;
end

xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation

no_focals=size(trajectory,2);

no_frames=size(trajectory,1);

if strcmp(method,'gauss')
    heatmapsize=[numel(edges{1})-1,numel(edges{2})-1];
else
    heatmapsize=[(size(edges{2},2)-1)*spacing-spacing+1 ...
        (size(edges{1},2)-1)*spacing-spacing+1];
    
end

no_maps=2;

if ndims(trajectory)==3
    no_dim=size(trajectory,3);
    tr2=NaN(no_frames,no_focals,no_neighbors,no_dim);
    for ff=1:no_focals
        tr2(:,ff,:,:)=repmat(trajectory(:,ff,:),[1 no_focals 1]);
    end
    trajectory=tr2;
    clear tr2
else
    no_dim=size(trajectory,4);
end
no_neighbors=size(trajectory,3);
%% Prepare

if no_dim==2
    trtemp=zeros(no_frames,no_focals,no_neighbors,3);
    trtemp(:,:,:,1:no_dim)=trajectory;
    trajectory=trtemp;
end
trajectory(:,:,:,2)=-trajectory(:,:,:,2); % This is a dirty trick in order to make the orientation of the coordinate system coherent with the videos

vel=NaN(size(trajectory));
acc=NaN(size(trajectory));
foc_to_nb_vec=NaN(no_frames,no_focals,no_neighbors,3);
ang_foc_vel_norm_and_nb_vel_norm=NaN(no_frames,no_focals,no_neighbors);

vel(1:no_frames-1,:,:,:)=diff(trajectory,1,1);
acc(1:no_frames-2,:,:,:)=diff(trajectory,2,1);
vel_magn=sqrt(sum(vel.^2,4));
vel_norm=bsxfun(@rdivide,vel,vel_magn);


for ff=1:no_focals
    for nf=1:no_neighbors
        if ff~=nf
            foc_to_nb_vec(:,nf,ff,:)=squeeze(trajectory(:,nf,ff,:)-trajectory(:,ff,ff,:));
            ang_foc_vel_norm_and_nb_vel_norm(:,ff,nf)=...
                mod(atan2(vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,2)-vel_norm(:,ff,ff,1).*vel_norm(:,nf,ff,2),vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,1)+vel_norm(:,nf,ff,2).*vel_norm(:,ff,ff,2)),2*pi);
        end
    end
end
% trajectory(:,:,:,2)=-trajectory(:,:,2); % This is a dirty trick in order to make the orientation of the coordinate system coherent with the videos

ang_vel_norm_z_plane_and_y_axis=...
    mod(atan2(vel_norm(:,:,:,1)*yax(2)-yax(1)*vel_norm(:,:,:,2),vel_norm(:,:,:,1)*yax(1)+vel_norm(:,:,:,2)*yax(2)),2*pi);


ang_vel_norm_z_plane_and_z_axis=...
    acos(vel_norm(:,:,:,1)*zax(1)+vel_norm(:,:,:,2)*zax(2)+vel_norm(:,:,:,3)*zax(3));
ang_vel_norm_and_xy_plane=NaN(size(ang_vel_norm_z_plane_and_z_axis));
ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis>=0)=...
    pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis>=0);
ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis<0)=...
    -pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis<0);
%% Measures which depend on relative position/orientation of foc and nb: Rotate vectors
a1=NaN(no_frames,no_focals);
for ff=1:no_focals
    a1(:,ff)=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,ff,ff,:));
end
a1=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,:,1,:));
if rotate_z
    a2=NaN(no_frames,no_focals);
    for ff=1:no_focals
        a2(:,ff)=squeeze(ang_vel_norm_and_xy_plane(:,ff,ff,:));
    end
else
    a2=zeros(size(a1));
end

a3=zeros(size(a1));
turnmatr=NaN(3,3,size(a1,1),no_focals);

turnmatr(1,1,:,:)=cos(a3).*cos(a1);
turnmatr(1,2,:,:)=-cos(a2).*sin(a1)+cos(a1).*sin(a3).*sin(a2);
turnmatr(1,3,:,:)=sin(a2).*sin(a1)+cos(a2).*sin(a3).*cos(a1);

turnmatr(2,1,:,:)=cos(a3).*sin(a1);
turnmatr(2,2,:,:)=sin(a2).*sin(a3).*sin(a1)+cos(a2).*cos(a1);
turnmatr(2,3,:,:)=-sin(a2).*cos(a1)+cos(a2).*sin(a3).*sin(a1);

turnmatr(3,1,:,:)=-sin(a3);
turnmatr(3,2,:,:)=sin(a2).*cos(a3);
turnmatr(3,3,:,:)=cos(a2).*cos(a3);

foc_to_nb_vec_rot=NaN(no_frames,no_focals,no_neighbors,3);
acc_rot=NaN(no_frames,no_focals,no_neighbors,3);
vel_rot=NaN(no_frames,no_focals,no_neighbors,3);

vel_trans=NaN(size(vel_rot));
vel_trans(:,:,:,1)=vel_rot(:,:,:,1);

for ff=1:no_focals
    for nf=1:no_neighbors
        %         if ff~=nf
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
foc_to_nb_vec_rot_magn=sqrt(sum(foc_to_nb_vec_rot.^2,4));
foc_to_nb_vec_dir_rot=bsxfun(@rdivide,foc_to_nb_vec_rot,foc_to_nb_vec_rot_magn);
focal_a_rotated=acc_rot;
focal_v_rotated=vel_rot;
acc_rot_magn=sqrt(sum(acc_rot.^2,4));
acc_rot_norm=bsxfun(@rdivide,acc_rot,acc_rot_magn);
%%
heatmap=NaN(no_focals,no_neighbors,no_maps+1,heatmapsize(1),heatmapsize(2));
maps=NaN(no_focals,no_neighbors,no_maps+1,heatmapsize(1),heatmapsize(2));
stdmaps=NaN(no_focals,no_neighbors,no_maps+1,heatmapsize(1),heatmapsize(2));
indices_in_bin=cell(no_focals,no_neighbors);

if nargout==2*no_maps+2+1
    maps_allframes=NaN(no_focals,no_neighbors,1,no_frames,heatmapsize(1),heatmapsize(2));
end

if strcmp(method,'gauss')
    ctrsx=edges{1}(1:end-1)+(edges{1}(2)-edges{1}(1))/2;
    ctrsy=edges{2}(1:end-1)+(edges{2}(2)-edges{1}(1))/2;
    [grx, gry]=meshgrid(ctrsx,ctrsy);
    edgesArrayX=repmat(reshape(grx,[1 numel(ctrsx) numel(ctrsy)]),[no_frames 1]);
    edgesArrayY=repmat(reshape(gry,[1 numel(ctrsx) numel(ctrsy)]),[no_frames 1]);
    [ix, iy]=meshgrid(1:numel(ctrsx),1:numel(ctrsy));
end

map_vallist=cell(no_focals,no_neighbors,no_maps+1,heatmapsize(1),heatmapsize(2));
for focal=1:no_focals
    for neighbour=1:no_neighbors
        if focal~=neighbour
            
            % Calculate 'frame' coordinates which define the x- and y-axis of
            % the map
            
            coordinates=squeeze(focal_to_neighbour_vector_rotated(:,neighbour,focal,1:2))/bodylength;
            %# Debug:
            %             coordinates=[-1*rand(no_frames,1) rand(no_frames,1)]*500;
            
            mapval=NaN(no_maps,no_frames);
            mapval(1,:)=sum(squeeze(focal_a_rotated(:,focal,focal,1:2))'*framerate^2/bodylength.*squeeze(focal_to_neighbour_vector_rotated(:,neighbour,focal,1:2))',1);
            if ~angles_2pi
                mapval(2,:)=acos(sum(squeeze(acc_rot_norm(:,focal,focal,1:2))'.*squeeze(foc_to_nb_vec_dir_rot(:,neighbour,focal,1:2))',1));
            else
                a=squeeze(acc_rot_norm(:,focal,focal,1:2))';
                b=squeeze(foc_to_nb_vec_dir_rot(:,neighbour,focal,1:2))';
                mapval(2,:)=...
                    atan2(a(1,:).*b(2,:)-b(1,:).*a(2,:),a(1,:).*b(1,:)+a(2,:).*b(2,:));
            end
            %             mapval(1:2,:)=rand(2,no_frames)*20;
            %             mapval(3:4,:)=rand(2,no_frames)*20;
            %
            
            mapval_new=mapval;
            for mapno=1:no_maps
                mapval_new(mapno,mapval_new(mapno,:)<prctile(mapval_new(mapno,:),discard_percentage_outliers))=NaN;
                mapval_new(mapno,mapval_new(mapno,:)>prctile(mapval_new(mapno,:),100-discard_percentage_outliers))=NaN;
            end
            
            
            
            
            % Note: 'Switch' in idx-logic: focal, neighbor
            % will mean 'focal reacts to neighbor' instead
            % of 'focal in neighbor system'
            if strcmp(method,'gauss')
                
                
                cx=repmat(coordinates(:,1),[1 numel(ctrsx) numel(ctrsy)]);
                cy=repmat(coordinates(:,2),[1 numel(ctrsx) numel(ctrsy)]);
                tic
                prob=exp(...
                    -.5*((cx-edgesArrayX)/spacing).^2 + ...
                    -.5*((cy-edgesArrayY)/spacing).^2 ...
                    );
                
                maps_temp=cell(no_maps,1);
                for mp = 1:no_maps
                    maps_temp{mp} = prob.*repmat(mapval_new(mp,:)',[1 numel(ctrsx) numel(ctrsy)]);
                end
                
                
                %                 maps(focal,neighbour,1,:,:)=squeeze(nanmean(prob,1));
                %                 for mapno=2:no_maps+1
                %                     maps(focal,neighbour,mapno,:,:)=squeeze(nanmean(maps_temp{mapno-1},1));
                %                     stdmaps(focal,neighbour,mapno,:,:)=squeeze(nanstd(maps_temp{mapno-1},1));
                %                 end
                
               
                
                map_vallist(focal,neighbour,1,:,:)=num2cell(squeeze(nanmean(prob,1)));
                
                
                for mapno=2:no_maps+1
                    A=maps_temp{mapno-1};
%                     A(abs(A)<10^-8)=NaN;
                    
                    %                     map_vallist(focal,neighbour,mapno,:,:)=...
                    %                         arrayfun(@(k,m) A(~isnan(A(:,k,m)),k,m),iy,ix,'UniformOutput',false);
                    
                    for k=1:numel(ctrsx)
                        for m=1:numel(ctrsy)
                            map_vallist(focal,neighbour,mapno,k,m)={A(~isnan(A(:,k,m)),k,m)};
                        end
                    end
                    
                    
                end
                
                if nargout==2*no_maps+2+1
                    maps_allframes(focal,neighbour,1,:,:,:)=prob;
                    for mapno=2:2
                        A=maps_temp{mapno-1};
%                         A(abs(A)<10^-8)=NaN;
                        maps_allframes(focal,neighbour,mapno,:,:,:)=A;
                    end
                end
                toc
            else
                [heatmap(focal,neighbour,1,:,:),...
                    ~, ~, maps_temp,...
                    stdmaps_temp,...
                    indices_in_bin{focal,neighbour}]=...
                    idSocial_hist3indexmap_mosaic(coordinates,'map',mapval_new','edges',edges,'spacing',spacing);
                
                maps(focal,neighbour,1,:,:)=heatmap(focal,neighbour,1,:,:)/no_frames/...
                    ((edges{1}(2)-edges{1}(1))*(edges{2}(2)-edges{2}(1)));
                
                
                
                for mapno=2:no_maps+1
                    maps(focal,neighbour,mapno,:,:)=maps_temp{mapno-1};
                    stdmaps(focal,neighbour,mapno,:,:)=stdmaps_temp{mapno-1};
                end
            end
        end
    end
end

if ~strcmp(method,'gauss')
    %% Generate 'map_vallist', a no_fish x no_fish x no_maps x edges x edges cell array,
    % which contains the indices (map no.1) and corresponding mapped values (e.g. velocity at index)
    map_vallist=cell(no_focals,no_neighbors,no_maps+1,heatmapsize(1),heatmapsize(2));
    for ff=1:size(indices_in_bin,1)
        for nf=1:size(indices_in_bin,2)
            if ff~=nf
                map_vallist(ff,nf,1,:,:)=cellfun(@(x) length(x)/no_frames, indices_in_bin{ff,nf},'UniformOutput',false);
                
            end
        end
    end
    mapval=NaN(no_maps,no_frames);
    for mapno=2:no_maps+1
        for ff=1:size(indices_in_bin,1)
            for nf=1:size(indices_in_bin,2)
                if ff~=nf
                    act=indices_in_bin{ff,nf};
                    %%%
                    
                    mapval(1,:)=sum(squeeze(focal_a_rotated(:,ff,ff,1:2))'*framerate^2/bodylength.*squeeze(focal_to_neighbour_vector_rotated(:,nf,ff,1:2))',1);
                    if ~angles_2pi
                        mapval(2,:)=acos(sum(squeeze(acc_rot_norm(:,ff,ff,1:2))'.*squeeze(foc_to_nb_vec_dir_rot(:,nf,ff,1:2))',1));
                    else
                        a=squeeze(acc_rot_norm(:,ff,ff,1:2))';
                        b=squeeze(foc_to_nb_vec_dir_rot(:,nf,ff,1:2))';
                        mapval(2,:)=...
                            atan2(a(1,:).*b(2,:)-b(1,:).*a(2,:),a(1,:).*b(1,:)+a(2,:).*b(2,:));
                    end
                    %                                         keyboard
                    mapval_new=mapval;
                    mapval_new(mapno-1,mapval_new(mapno-1,:)<prctile(mapval_new(mapno-1,:),discard_percentage_outliers))=NaN;
                    mapval_new(mapno-1,mapval_new(mapno-1,:)>prctile(mapval_new(mapno-1,:),100-discard_percentage_outliers))=NaN;
                    %%%
                    
                    try
%                         keyboard
                        map_vallist(ff,nf,mapno,:,:)=cellfun(@(x) mapval_new(mapno-1,x(~isnan(x)))',act,'UniformOutput',false);
                        
                    catch
                        keyboard
                    end
                end
            end
        end
    end
end
%% Prepare output
mapcell=cell(no_maps,1);
varargout=cell(2*(no_maps+1),1);
stdcell=cell(no_maps,1);
for mp=1:no_maps+1
    mapcell{mp}=squeeze(maps(:,:,mp,:,:));
    stdcell{mp}=squeeze(stdmaps(:,:,mp,:,:));
    varargout{mp}=squeeze(maps(:,:,mp,:,:));
end
for mp=no_maps+2:2*no_maps+2
    varargout{mp}=squeeze(map_vallist(:,:,mp-no_maps-1,:,:));
end
varargout{2*no_maps+2+1}=(nanmean(maps(:,:,2,:),4));
varargout{2*no_maps+2+2}=(nanmean(maps(:,:,3,:),4));
if nargout==2*no_maps+2+3
    varargout{2*no_maps+2+3}=maps_allframes;
end
% varargout{mp+1}=indices_in_bin;
% varargout=mapcell;
%% For debugging
% ff=1; nf=2; fr=20093;
%
% vsc=150;
% asc=100;
%
% f2n=squeeze(foc_to_nb_vec(fr,nf,ff,:));
% v=squeeze(vel(fr,ff,ff,:));
% a=squeeze(acc(fr,ff,ff,:));
%
% f2nr=squeeze(focal_to_neighbour_vector_rotated(fr,nf,ff,:));
% vr=squeeze(vel_rot(fr,ff,ff,:));
% ar=squeeze(acc_rot(fr,ff,ff,:));
%
% figure;
% plot([0 f2n(1)],[0 f2n(2)],'r');
% hold on
% plot([0 v(1)*vsc],[0 v(2)*vsc],':b');
% plot([0 a(1)*asc],[0 a(2)*asc],':k');
%
% plot([0 f2nr(1)],[0 f2nr(2)],'g');
% plot([0 vr(1)*vsc],[0 vr(2)*vsc],'.-b');
% plot([0 ar(1)*asc],[0 ar(2)*asc],'.-k');
%
% a1(fr,ff)/pi*180
% turnmatr(:,:,fr,ff)
%
% axis equal