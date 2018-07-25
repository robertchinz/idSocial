function [movementdata, options]=idSocial_preparedata3D(trayectorias,options)
% Calculate dynamic data from trajectory.
% def_options.smooth_method='moving';
% def_options.smooth_degree=3;
% def_options.interpolation_mode='spline';
% def_options.interp_maxlen=inf;
% def_options.interpolate_polar_angle_and_radius=false;
% def_options.polar_center=[0 0 0];
% def_options.polar_axis=[xax; yax; zax];
% def_options.exclusion_radius=0;
% def_options.params_temporal_distance=[1 1 1];
% def_options.preparedata_centerofmass=false;
% def_options.arena_limits=[];

% Output: structure MOVEMENTDATA with fields:


if nargin<2 || isempty(options)
    options=[];
end
xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation
no_dim=3;
%%
act_method=mfilename;
act_method=act_method(10:end);

def_options.act_method=act_method;%(10:end);


% Def. options for dynamic data
def_options.smooth_method='moving';
def_options.smooth_degree=3;
def_options.interpolation_mode='spline';
def_options.interp_maxlen=inf;
def_options.interpolate_polar_angle_and_radius=false;
def_options.polar_center=[0 0 0];
def_options.polar_axis=[xax; yax; zax];
def_options.exclusion_radius=0;
def_options.params_temporal_distance=[1 1 1];
def_options.preparedata_centerofmass=false;
def_options.arena_limits=[];
def_options.rotate_z=true;
def_options.start_min=1;
def_options.end_min=inf;

defoptnames=fieldnames(def_options);
no_defoptions=size(defoptnames,1);
for optfield=1:no_defoptions
    if ~isfield(options,defoptnames{optfield}) %% || isempty(options.(defoptnames{optfield}))
        options.(defoptnames{optfield})=def_options.(defoptnames{optfield});
    end
end
disp(options)
tic;
%% Options for dynamics: 
tr=trayectorias;
polar_center=options.polar_center;
interpolate_polar_angle_and_radius=options.interpolate_polar_angle_and_radius;
polar_axis=options.polar_axis;
exclusion_radius=options.exclusion_radius;
params_temporal_distance=options.params_temporal_distance;
preparedata_centerofmass=options.preparedata_centerofmass;
arena_limits=options.arena_limits;
rotate_z=options.rotate_z;
idces=[options.start_min min(options.end_min,size(trayectorias,1))];
%%
no_frames=size(trayectorias,1);
no_fish=size(trayectorias,2);
no_dim_orig=size(trayectorias,4);
no_vertices=8;

if no_dim_orig==2
    trtemp=zeros(no_frames,no_fish,no_fish,3);
    trtemp(:,:,:,1:no_dim_orig)=trayectorias;
    tr=trtemp;
else
    tr=trayectorias;
end   

if isempty(arena_limits)
    minx=floor(min(min(min(trayectorias(:,:,1))))*1);
    miny=floor(min(min(min(trayectorias(:,:,2))))*1);
    minz=floor(min(min(min(trayectorias(:,:,3))))*1);
    maxx=floor(max(max(max(trayectorias(:,:,1))))*1);
    maxy=floor(max(max(max(trayectorias(:,:,2))))*1);
    maxz=floor(max(max(max(trayectorias(:,:,3))))*1);
    arena_limits=[minx miny minz; ...
        minx maxy minz; ...
        maxx maxy minz; ...
        maxx miny minz; ...
        minx miny maxz; ...
        minx maxy maxz; ...
        maxx maxy maxz; ...
        maxx miny maxz];
    options.arena_limits=arena_limits;
    
    
end
%% Interpolation and smoothing either from polar or cartesian coordinates
if interpolate_polar_angle_and_radius % Transform to polar/cylindrical coordinates, then interpolate (and smooth, but smoothing is disabled right now...disculpen las molestias) radius and angle.
     
    foc_to_polarcenter=cat(4,polar_center(1)-tr(:,:,:,1), polar_center(2)-tr(:,:,:,2), polar_center(3)-tr(:,:,:,3));
    foc_to_polarcenter_dir=bsxfun(@rdivide,foc_to_polarcenter,sqrt(sum(foc_to_polarcenter.^2,4)));
    foc_polar_radius=sqrt(sum(foc_to_polarcenter.^2,4));
    
    foc_polar_angle_xz=...
            atan2(polar_axis(1,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(1,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))),...
            (polar_axis(1,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(1,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))))...
            );
    foc_polar_angle_z=...
        atan2(-foc_to_polarcenter_dir(:,:,:,1)*zax(2)-zax(1)*(-foc_to_polarcenter_dir(:,:,:,2)),(-foc_to_polarcenter_dir(:,:,:,1))*zax(1)+(-foc_to_polarcenter_dir(:,:,:,2))*zax(2));
    
    foc_polar_angle_xy=NaN(size(foc_polar_angle_z));
    foc_polar_angle_xy(foc_polar_angle_z>=0)=...
        pi/2-foc_polar_angle_z(foc_polar_angle_z>=0);
    foc_polar_angle_xy(foc_polar_angle_z<0)=...
        -pi/2-foc_polar_angle_z(foc_polar_angle_z<0);
    % The following would calculate the angle in "real" 3D
    % foc_polar_angle_xz=...
    %     atan2(sqrt(...
    %         (polar_axis(2).*(-squeeze(foc_to_polarcenter_dir(:,ff,3))) - polar_axis(3).*(-squeeze(foc_to_polarcenter_dir(:,ff,2)))).^2+...
    %         (polar_axis(3).*(-squeeze(foc_to_polarcenter_dir(:,ff,1))) - polar_axis(1).*(-squeeze(foc_to_polarcenter_dir(:,ff,3)))).^2+...
    %         (polar_axis(1).*(-squeeze(foc_to_polarcenter_dir(:,ff,2))) - polar_axis(2).*(-squeeze(foc_to_polarcenter_dir(:,ff,1)))).^2)',...
    %         polar_axis*(-squeeze(foc_to_polarcenter_dir(:,ff,:))')...
    %     );
    
    if ~isempty(interp_mode)
        disp([act_method ': Interpolation using polar coordinates.'])
        for ff=1:no_fish
            foc_polar_angle_xz(:,ff,:)=repmat(interp_angles(foc_polar_angle_xz(:,ff,ff)),[1,1,no_fish]);
            foc_polar_angle_xy(:,ff,:)=repmat(interp_angles(foc_polar_angle_xy(:,ff,ff)),[1,1,no_fish]);
        end
        
        foc_polar_radius=idSocial_interpolateTrajectories(foc_polar_radius,interp_mode,interp_maxlen);
    end
    disp('¡¡¡ATTENTION!!! No moving averages in polar coordinates. Function mediamovil4angles is faulty. Fix it!')
    
    foc_to_polarcenter=tr;
  
else 
    foc_to_polarcenter=cat(4,polar_center(1)-tr(:,:,:,1), polar_center(2)-tr(:,:,:,2), polar_center(3)-tr(:,:,:,3));
    foc_to_polarcenter_dir=bsxfun(@rdivide,foc_to_polarcenter,sqrt(sum(foc_to_polarcenter.^2,4)));
    foc_polar_radius=sqrt(sum(foc_to_polarcenter.^2,4));
    foc_polar_angle_xz=...
        mod(atan2(polar_axis(1,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(1,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))),...
        (polar_axis(1,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(1,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))))),2*pi...
        );
    foc_polar_angle_yz=...
        mod(atan2(polar_axis(2,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(2,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))),...
        (polar_axis(2,1).*(-squeeze(foc_to_polarcenter_dir(:,:,:,1))) - polar_axis(2,2).*(-squeeze(foc_to_polarcenter_dir(:,:,:,2))))),2*pi...
        );
    foc_polar_angle_z=...
        atan2(-foc_to_polarcenter_dir(:,:,:,1)*zax(2)-zax(1)*(-foc_to_polarcenter_dir(:,:,:,2)),(-foc_to_polarcenter_dir(:,:,:,1))*zax(1)+(-foc_to_polarcenter_dir(:,:,:,2))*zax(2));
    foc_polar_angle_xy=NaN(size(foc_polar_angle_z));
    foc_polar_angle_xy(foc_polar_angle_z>=0)=...
        pi/2-foc_polar_angle_z(foc_polar_angle_z>=0);
    foc_polar_angle_xy(foc_polar_angle_z<0)=...
        -pi/2-foc_polar_angle_z(foc_polar_angle_z<0);
%     foc_polar_angle_xz=mediamovil4angles(foc_polar_angle_xz,interp_maxlen);
%     foc_polar_radius= moving_average(interp_maxlen,foc_polar_radius,[],2);
end
%% Some basic measures (velocity, acceleration, angles to axes) which depend only on the foc itself. 
vel=NaN(size(tr));
acc=NaN(size(tr));
foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,no_dim);
% foc_to_vertex_vec=NaN(no_frames,no_fish,no_vertices,no_dim);
ang_foc_vel_norm_and_nb_vel_norm=NaN(no_frames,no_fish,no_fish);

vel(1:no_frames-1,:,:,:)=diff(tr,1,1);
acc(1:no_frames-2,:,:,:)=diff(tr,2,1);
vel_magn=sqrt(sum(vel.^2,4));
vel_norm=bsxfun(@rdivide,vel,vel_magn);
acc_magn=sqrt(sum(acc.^2,4));
acc_norm=bsxfun(@rdivide,acc,acc_magn);

% for ff=1:no_fish
%     for vt=1:no_vertices
%         foc_to_vertex_vec(:,ff,vt,:)=repmat(arena_limits(vt,:),[no_frames,1])-squeeze(tr(:,ff,1,:));  
%     end
% end

for ff=1:no_fish
    for nf=1:no_fish
        foc_to_nb_vec(:,ff,nf,:)=squeeze(tr(:,nf,1,:)-tr(:,ff,1,:));  
        ang_foc_vel_norm_and_nb_vel_norm(:,ff,nf)=...
            mod(atan2(vel_norm(:,nf,1,1).*vel_norm(:,ff,1,2)-vel_norm(:,ff,1,1).*vel_norm(:,nf,1,2),vel_norm(:,nf,1,1).*vel_norm(:,ff,1,1)+vel_norm(:,nf,1,2).*vel_norm(:,ff,1,2)),2*pi);
    end
end
distance=sqrt(sum(foc_to_nb_vec.^2,4));

foc_dist_to_origin=sqrt(sum(tr.^2,4));
foc_pos_dir=bsxfun(@rdivide,tr,foc_dist_to_origin);
foc_to_nb_dir=foc_to_nb_vec./repmat(sqrt(sum(foc_to_nb_vec.^2,4)),[1,1,1,3]);
% foc_to_vertex_dir=foc_to_vertex_vec./repmat(sqrt(sum(foc_to_vertex_vec.^2,4)),[1,1,1,3]);
% distance_vertex=sqrt(sum(foc_to_vertex_vec.^2,4));

ang_vel_norm_z_plane_and_y_axis=...
    mod(atan2(vel_norm(:,:,:,1)*yax(2)-yax(1)*vel_norm(:,:,:,2),vel_norm(:,:,:,1)*yax(1)+vel_norm(:,:,:,2)*yax(2)),2*pi);
ang_vel_norm_z_plane_and_z_axis=...
    acos(vel_norm(:,:,:,1)*zax(1)+vel_norm(:,:,:,2)*zax(2)+vel_norm(:,:,:,3)*zax(3));
 


ang_vel_norm_and_xy_plane=NaN(size(ang_vel_norm_z_plane_and_z_axis));
ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis>=0)=...
    pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis>=0);
ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis<0)=...
    -pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis<0);

ang_el_foc_vel_norm_and_nb_vel_norm=...
    strucdyn_angle_difference(ang_vel_norm_and_xy_plane,permute(ang_vel_norm_and_xy_plane,[1 3 2]));
%% Misc. angels
foc_exclusion_distance=sqrt(foc_polar_radius.^2-exclusion_radius^2);
foc_exclusion_angle=asin(bsxfun(@rdivide,exclusion_radius,foc_polar_radius));
foc_polar_angle_xz_vel=vertcat(strucdyn_angle_difference(foc_polar_angle_xz(1:end-1,:,:),foc_polar_angle_xz(2:end,:,:)), NaN(1,no_fish,no_fish));
foc_polar_angle_xz_acceleration=vertcat(diff(foc_polar_angle_xz_vel,1,1), NaN(1,no_fish,no_fish));
foc_polar_radius_vel=vertcat(diff(foc_polar_radius,1,1), NaN(1,no_fish,no_fish));

angle_foc_to_polar_center_and_foc_nb_dir=NaN(no_frames,no_fish,no_fish);
polar_angle_distance_foc_to_nb=NaN(no_frames,no_fish,no_fish);
ang_foc_vel_norm_z_plane_and_nb_vel_norm=NaN(no_frames,no_fish,no_fish);
ang_foc_vel_norm_z_plane_and_nb_pos=NaN(no_frames,no_fish,no_fish);
ang_foc_vel_norm_z_plane_and_vertex_pos=NaN(no_frames,no_fish,no_vertices);
ang_foc_acc_norm_z_plane_and_nb_pos=NaN(no_frames,no_fish,no_fish);
ang_foc_acc_norm_z_plane_and_vel_norm=NaN(no_frames,no_fish,no_fish);
ang_foc_vel_norm_z_plane_t1_and_t2=NaN(no_frames,no_fish,no_fish);
for ff=1:no_fish
    
    ang_foc_acc_norm_z_plane_and_vel_norm(:,ff,:)=...
        repmat(...
        atan2(vel_norm(:,ff,1,1).*acc_norm(:,ff,1,2)-acc_norm(:,ff,1,1).*vel_norm(:,ff,1,2),acc_norm(:,ff,1,1).*vel_norm(:,ff,1,1)+acc_norm(:,ff,1,2).*vel_norm(:,ff,1,2)),...
        1,no_fish);
    
    ang_foc_vel_norm_z_plane_t1_and_t2(:,ff,:)=...
        repmat(...
        vertcat(atan2(vel_norm(2:end,ff,1,2).*vel_norm(1:end-1,ff,1,1)-vel_norm(1:end-1,ff,2).*vel_norm(2:end,ff,1,1),vel_norm(2:end,ff,1,1).*vel_norm(1:end-1,ff,1,1)+vel_norm(2:end,ff,1,2).*vel_norm(1:end-1,ff,1,2)), NaN),...
        1,no_fish);
    
    
    for nf=1:no_fish
        if ff~=nf
            
            angle_foc_to_polar_center_and_foc_nb_dir(:,ff,nf)=...
                atan2(foc_to_polarcenter_dir(:,ff,1,1).*foc_to_nb_dir(:,ff,nf,2)-foc_to_nb_dir(:,ff,nf,1).*foc_to_polarcenter_dir(:,ff,1,2),...
                foc_to_polarcenter_dir(:,ff,1,1).*foc_to_nb_dir(:,ff,nf,1)+foc_to_polarcenter_dir(:,ff,1,2).*foc_to_nb_dir(:,ff,nf,2));
            
            polar_angle_distance_foc_to_nb(:,ff,nf)=strucdyn_angle_difference(foc_polar_angle_xz(:,ff,1),foc_polar_angle_xz(:,nf,1));
            
            ang_foc_vel_norm_z_plane_and_nb_vel_norm(:,ff,nf)=...
                mod(atan2(vel_norm(:,nf,1,1).*vel_norm(:,ff,1,2)-vel_norm(:,ff,1,1).*vel_norm(:,nf,1,2),vel_norm(:,ff,1,1).*vel_norm(:,nf,1,1)+vel_norm(:,ff,1,2).*vel_norm(:,nf,1,2)),2*pi);
        
            ang_foc_vel_norm_z_plane_and_nb_pos(:,ff,nf)=...
                mod(atan2(foc_to_nb_dir(:,ff,nf,1).*vel_norm(:,ff,1,2)-vel_norm(:,ff,1,1).*foc_to_nb_dir(:,ff,nf,2),vel_norm(:,ff,1,1).*foc_to_nb_dir(:,ff,nf,1)+vel_norm(:,ff,1,2).*foc_to_nb_dir(:,ff,nf,2)),2*pi);

            ang_foc_acc_norm_z_plane_and_nb_pos(:,ff,nf)=...
                mod(atan2(foc_to_nb_dir(:,ff,nf,1).*acc_norm(:,ff,1,2)-acc_norm(:,ff,1,1).*foc_to_nb_dir(:,ff,nf,2),acc_norm(:,ff,1,1).*foc_to_nb_dir(:,ff,nf,1)+acc_norm(:,ff,1,2).*foc_to_nb_dir(:,ff,nf,2)),2*pi);

            
        end
    end
%     for vt=1:no_vertices
%         ang_foc_vel_norm_z_plane_and_vertex_pos(:,ff,vt)=...
%                 mod(atan2(foc_to_vertex_dir(:,ff,vt,1).*vel_norm(:,ff,1,2)-vel_norm(:,ff,1,1).*foc_to_vertex_dir(:,ff,vt,2),vel_norm(:,ff,1,1).*foc_to_vertex_dir(:,ff,vt,1)+vel_norm(:,ff,1,2).*foc_to_vertex_dir(:,ff,vt,2)),2*pi);
% 
%     end
    
end

ang_foc_vel_norm_z_plane_t1_and_t2_magn=abs(ang_foc_vel_norm_z_plane_t1_and_t2);
ang_acc_foc_vel_norm_z_plane_t1_and_t2=vertcat(diff(ang_foc_vel_norm_z_plane_t1_and_t2,1,1), NaN(1,no_fish,no_fish));

%% Temporal distance
% temporal_distance=distance/params_temporal_distance(1) + ...
%                     min(ang_foc_vel_norm_z_plane_and_nb_pos,2*pi-ang_foc_vel_norm_z_plane_and_nb_pos)*180/pi/params_temporal_distance(2)+ ...
%                     params_temporal_distance(3);          
%% Measures which depend on relative position/orientation of foc and nb: Rotate vectors
a1=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,:,1,:));
if rotate_z
    a2=squeeze(ang_vel_norm_and_xy_plane(:,:,1,:));
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

foc_to_nb_vec_rot=NaN(no_frames,no_fish,no_fish,no_dim);
foc_to_vertex_vec_rot=NaN(no_frames,no_fish,no_vertices,no_dim);
acc_rot=NaN(no_frames,no_fish,no_fish,no_dim);
vel_rot=NaN(no_frames,no_fish,no_fish,no_dim);
foc_to_nb_dir_rot=bsxfun(@rdivide,foc_to_nb_vec_rot,sqrt(sum(foc_to_nb_vec_rot.^2,4)));

vel_trans=NaN(size(vel_rot));
vel_trans(:,:,:,1)=vel_rot(:,:,:,1);

for ff=1:no_fish
    for nf=1:no_fish
        foc_to_nb_vec_rot(:,ff,nf,:)=...
            [sum(squeeze(foc_to_nb_vec(:,ff,nf,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
            sum(squeeze(foc_to_nb_vec(:,ff,nf,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
            sum(squeeze(foc_to_nb_vec(:,ff,nf,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
        vel_rot(:,ff,nf,:)=...
            [sum(squeeze(vel(:,ff,1,:)).*squeeze(turnmatr(:,1,:,nf))',2) ...
            sum(squeeze(vel(:,ff,1,:)).*squeeze(turnmatr(:,2,:,nf))',2) ...
            sum(squeeze(vel(:,ff,1,:)).*squeeze(turnmatr(:,3,:,nf))',2)];
        acc_rot(:,ff,nf,:)=...
            [sum(squeeze(acc(:,ff,1,:)).*squeeze(turnmatr(:,1,:,nf))',2) ...
            sum(squeeze(acc(:,ff,1,:)).*squeeze(turnmatr(:,2,:,nf))',2) ...
            sum(squeeze(acc(:,ff,1,:)).*squeeze(turnmatr(:,3,:,nf))',2)];
        vel_trans(:,ff,nf,2)=vel_rot(:,ff,nf,2)-vel_rot(:,ff,ff,2);
    end
%     for vt=1:no_vertices
%         foc_to_vertex_vec_rot(:,ff,vt,:)=...
%             [sum(squeeze(foc_to_vertex_vec(:,ff,vt,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
%             sum(squeeze(foc_to_vertex_vec(:,ff,vt,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
%             sum(squeeze(foc_to_vertex_vec(:,ff,vt,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
%     end
end
vel_norm_rot=bsxfun(@rdivide,vel_rot,sqrt(sum(vel_rot.^2,4)));
acc_norm_rot=bsxfun(@rdivide,acc_rot,sqrt(sum(acc_rot.^2,4)));
acc_trans=vertcat(diff(vel_trans,1,1), NaN(1,no_fish,no_fish,3));
%% ----------------Integrating over
% %angle_focal_movement_direction_t1_and_t2 until direction
% %changes---------------------------------------------------
% ang_foc_vel_norm_z_plane_t1_and_t2_integrated=NaN(size(ang_foc_vel_norm_z_plane_t1_and_t2));
% for ff=1:no_fish
%     da=ang_foc_vel_norm_z_plane_t1_and_t2(end:-1:1,ff,1);
%     da(isnan(ang_foc_vel_norm_z_plane_t1_and_t2(end:-1:1,ff,1)))=0;
%     da(da<0)=0;
%     
%     cx1=cumsum(da,1);
%     
%     cx2=zeros(size(cx1));
%     i=strfind([0 da'>0.'],[0 1]);
%     cx2(i)=diff([0 cx1(i)'-da(i)']);
%     x3=cx1-cumsum(cx2);
%     % x3(~x) = NaN;
%     x3(da<=0)=0;
%     
%     da=ang_foc_vel_norm_z_plane_t1_and_t2(end:-1:1,ff);
%     da(isnan(ang_foc_vel_norm_z_plane_t1_and_t2(end:-1:1,ff)))=0;
%     da(da>0)=0;
%     cx1=cumsum(da);
%     
%     cx3=zeros(size(cx1));
%     i=strfind([0 da'<0.'],[0 1]);
%     cx3(i)=diff([0 cx1(i)'-da(i)']);
%     x4=cx1-cumsum(cx3);
%     % x4(~x(:,1)) = NaN;
%     
%     x4(da>=0)=0;
%     ang_foc_vel_norm_z_plane_t1_and_t2_integrated(:,ff,:)=repmat(x3,1,no_fish);
%     ang_foc_vel_norm_z_plane_t1_and_t2_integrated(x3==0,ff,:)=repmat(x4(x3==0),1,no_fish);
%     ang_foc_vel_norm_z_plane_t1_and_t2_integrated(:,ff,:)=repmat(ang_foc_vel_norm_z_plane_t1_and_t2_integrated(end:-1:1,ff),1,no_fish);
% end
%% Generate structure for output
movementdata.options=options;
movementdata.angle_focal_movement_direction_t1_and_t2=ang_foc_vel_norm_z_plane_t1_and_t2;
% movementdata.angle_focal_movement_direction_t1_and_t2_magnitude=ang_foc_vel_norm_z_plane_t1_and_t2_magn;
movementdata.angle_acceleration_focal_movement_direction_t1_and_t2=ang_acc_foc_vel_norm_z_plane_t1_and_t2;
% movementdata.angle_focal_movement_direction_t1_and_t2_integrated=ang_foc_vel_norm_z_plane_t1_and_t2_integrated;
% movementdata.trajectory_raw=tr_raw(:,:,:,1:no_dim_orig);
% movementdata.trajectory=tr(:,:,:,1:no_dim_orig);
movementdata.velocity=vel(:,:,:,1:no_dim_orig);
movementdata.velocity_magnitude=vel_magn;
movementdata.acc_magnitude=acc_magn;
movementdata.acc=acc(:,:,:,1:no_dim_orig);
% movementdata.angle_focal_movement_direction_and_y_axis=ang_vel_norm_z_plane_and_y_axis;
movementdata.angle_az_focal_movdir_and_nb_movdir=ang_foc_vel_norm_z_plane_and_nb_vel_norm;
movementdata.angle_el_focal_movdir_and_nb_movdir=ang_el_foc_vel_norm_and_nb_vel_norm;
movementdata.angle_focal_movement_direction_and_neighbour_position=ang_foc_vel_norm_z_plane_and_nb_pos;
% movementdata.angle_focal_movement_direction_and_vertex_position=ang_foc_vel_norm_z_plane_and_vertex_pos;
movementdata.focal_to_neighbour_vector_rotated=foc_to_nb_vec_rot(:,:,:,1:no_dim_orig);
% movementdata.focal_to_vertex_vector_rotated=foc_to_vertex_vec_rot(:,:,:,1:no_dim_orig);
movementdata.focal_to_neighbour_direction_rotated=foc_to_nb_dir_rot(:,:,:,1:no_dim_orig);
movementdata.velocity_rotated=vel_rot(:,:,:,1:no_dim_orig);
% movementdata.velocity_transformed=vel_trans(:,:,:,1:no_dim_orig);
movementdata.acc_rotated=acc_rot(:,:,:,1:no_dim_orig);
% movementdata.acc_transformed=acc_trans(:,:,:,1:no_dim_orig);
movementdata.distance_focal_neighbour=distance;
% movementdata.distance_focal_vertex=distance_vertex;
% movementdata.temporal_distance=temporal_distance;
movementdata.movement_direction=vel_norm(:,:,:,1:no_dim_orig);
movementdata.movement_direction_rotated=vel_norm_rot(:,:,:,1:no_dim_orig);
movementdata.acceleration_direction_rotated=acc_norm_rot(:,:,:,1:no_dim_orig);
movementdata.focal_to_neighbour_direction=foc_to_nb_dir(:,:,:,1:no_dim_orig);
movementdata.focal_to_neighbour_vector=foc_to_nb_vec(:,:,:,1:no_dim_orig);
% movementdata.focal_to_vertex_vector=foc_to_vertex_vec(:,:,:,1:no_dim_orig);
movementdata.angle_focal_acceleration_and_neighbour_position=ang_foc_acc_norm_z_plane_and_nb_pos;
movementdata.focal_acceleration_direction=acc_norm(:,:,:,1:no_dim_orig);
movementdata.angle_focal_acceleration_and_movement_direction=ang_foc_acc_norm_z_plane_and_vel_norm;
% movementdata.angle_focal_position_and_x_axis=foc_polar_angle_xz;
% movementdata.angle_focal_position_and_xy_plane=foc_polar_angle_xy;
% movementdata.angle_focal_position_and_yz_plane=foc_polar_angle_yz;
% movementdata.focal_position_direction=foc_pos_dir(:,:,:,1:no_dim_orig);
% movementdata.focal_distance_to_origin=foc_dist_to_origin;
% movementdata.focal_polar_radius=foc_polar_radius;
% movementdata.focal_polar_angle=foc_polar_angle_xz;
% movementdata.focal_polar_angle_velocity=foc_polar_angle_xz_vel;
% movementdata.focal_polar_angle_acceleration=foc_polar_angle_xz_acceleration;
% movementdata.focal_polar_radius_velocity=foc_polar_radius_vel;
% movementdata.focal_x_to_polarcenter_direction=foc_to_polarcenter_dir(:,:,:,1:no_dim_orig);
% movementdata.focal_x_to_polarcenter=foc_to_polarcenter(:,:,:,1:no_dim_orig);
% movementdata.focal_exclusion_distance=foc_exclusion_distance;
% movementdata.focal_exclusion_angle=foc_exclusion_angle;
% movementdata.angle_focal_to_polar_center_and_focal_neighbour_direction=angle_foc_to_polar_center_and_foc_nb_dir;
% movementdata.polar_angle_distance_focal_to_neighbour=polar_angle_distance_foc_to_nb;


%% And finish
disp(['[' mfilename '] ' 'Movementdata is ready! (It took ' num2str(toc) ' seconds to prepare it)'])
