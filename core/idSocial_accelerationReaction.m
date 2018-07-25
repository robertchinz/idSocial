function [turningReact, turningReactInBins, turningReactCell,turning]=idSocial_accelerationReaction(tr,spacing,edges,framerate,bodylength,normalize_by_individual_max,rotate_z,method)

no_fish=size(tr,2);
no_frames=size(tr,1);
no_dim=size(tr,4);

xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation

if nargin<4 || isempty(framerate)
    framerate=1;
end
if nargin<5 || isempty(bodylength)
    bodylength=1;
end

if nargin<6 || isempty(normalize_by_individual_max)
    normalize_by_individual_max=false;
end
if nargin<7 || isempty(rotate_z)
    rotate_z=false;
end
if nargin<8 || isempty(method)
    method='gauss';
end
if iscell(edges)
    edges=edges{1};
end
% %% Prepare
% 
% if no_dim==2
%     trtemp=zeros(no_frames,no_fish,no_fish,3);
%     trtemp(:,:,:,1:no_dim)=tr;
%     tr=trtemp;
% end
% tr(:,:,:,2)=-tr(:,:,:,2); % This is a dirty trick in order to make the orientation of the coordinate system coherent with the videos
% 
% vel=NaN(size(tr));
% acc=NaN(size(tr));
% foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,3);
% ang_foc_vel_norm_and_nb_vel_norm=NaN(no_frames,no_fish,no_fish);
% 
% vel(1:no_frames-1,:,:,:)=diff(tr,1,1);
% acc(1:no_frames-2,:,:,:)=diff(tr,2,1);
% vel_magn=sqrt(sum(vel.^2,4));
% vel_norm=bsxfun(@rdivide,vel,vel_magn);
% 
% 
% for ff=1:no_fish
%     for nf=1:no_fish
%         if ff~=nf
%             foc_to_nb_vec(:,nf,ff,:)=squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
%             ang_foc_vel_norm_and_nb_vel_norm(:,ff,nf)=...
%                 mod(atan2(vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,2)-vel_norm(:,ff,ff,1).*vel_norm(:,nf,ff,2),vel_norm(:,nf,ff,1).*vel_norm(:,ff,ff,1)+vel_norm(:,nf,ff,2).*vel_norm(:,ff,ff,2)),2*pi);
%         end
%     end
% end
% % tr(:,:,:,2)=-tr(:,:,2); % This is a dirty trick in order to make the orientation of the coordinate system coherent with the videos
% 
% ang_vel_norm_z_plane_and_y_axis=...
%     mod(atan2(vel_norm(:,:,:,1)*yax(2)-yax(1)*vel_norm(:,:,:,2),vel_norm(:,:,:,1)*yax(1)+vel_norm(:,:,:,2)*yax(2)),2*pi);
% 
% 
% ang_vel_norm_z_plane_and_z_axis=...
%     acos(vel_norm(:,:,:,1)*zax(1)+vel_norm(:,:,:,2)*zax(2)+vel_norm(:,:,:,3)*zax(3));
% ang_vel_norm_and_xy_plane=NaN(size(ang_vel_norm_z_plane_and_z_axis));
% ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis>=0)=...
%     pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis>=0);
% ang_vel_norm_and_xy_plane(ang_vel_norm_z_plane_and_z_axis<0)=...
%     -pi/2-ang_vel_norm_z_plane_and_z_axis(ang_vel_norm_z_plane_and_z_axis<0);
% %% Measures which depend on relative position/orientation of foc and nb: Rotate vectors
% a1=NaN(no_frames,no_fish);
% for ff=1:no_fish
%     a1(:,ff)=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,ff,ff,:));
% end
% a1=-squeeze(ang_vel_norm_z_plane_and_y_axis(:,:,1,:));
% if rotate_z
%     a2=NaN(no_frames,no_fish);
%     for ff=1:no_fish
%         a2(:,ff)=squeeze(ang_vel_norm_and_xy_plane(:,ff,ff,:));
%     end
% else
%     a2=zeros(size(a1));
% end
% 
% a3=zeros(size(a1));
% turnmatr=NaN(3,3,size(a1,1),no_fish);
% 
% turnmatr(1,1,:,:)=cos(a3).*cos(a1);
% turnmatr(1,2,:,:)=-cos(a2).*sin(a1)+cos(a1).*sin(a3).*sin(a2);
% turnmatr(1,3,:,:)=sin(a2).*sin(a1)+cos(a2).*sin(a3).*cos(a1);
% 
% turnmatr(2,1,:,:)=cos(a3).*sin(a1);
% turnmatr(2,2,:,:)=sin(a2).*sin(a3).*sin(a1)+cos(a2).*cos(a1);
% turnmatr(2,3,:,:)=-sin(a2).*cos(a1)+cos(a2).*sin(a3).*sin(a1);
% 
% turnmatr(3,1,:,:)=-sin(a3);
% turnmatr(3,2,:,:)=sin(a2).*cos(a3);
% turnmatr(3,3,:,:)=cos(a2).*cos(a3);
% 
% foc_to_nb_vec_rot=NaN(no_frames,no_fish,no_fish,3);
% acc_rot=NaN(no_frames,no_fish,no_fish,3);
% vel_rot=NaN(no_frames,no_fish,no_fish,3);
% 
% vel_trans=NaN(size(vel_rot));
% vel_trans(:,:,:,1)=vel_rot(:,:,:,1);
% 
% for ff=1:no_fish
%     for nf=1:no_fish
%         %         if ff~=nf
%         % Idx (:,ff,nf,:) means: Corresponding value (vel,acc,..) for focal ff 'seen by' neighbor nf;
%         % In case of foc_to_nb_vec, vel and acc, 'seen
%         % by' refers to neighbor-depend filtering (e.g.,
%         % value at (:,ff,nf,:) is NaN because neighbor moves slower than
%         % the lower threshold speed).
%         % For rotated values, (:,ff,nf,:) means: vector
%         % for focal ff in the system where neighbor velocity
%         % vector is parallel to the y-axis
%         foc_to_nb_vec_rot(:,nf,ff,:)=...
%             [sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
%             sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
%             sum(squeeze(foc_to_nb_vec(:,nf,ff,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
%         vel_rot(:,ff,nf,:)=...
%             [sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
%             sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
%             sum(squeeze(vel(:,ff,nf,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
%         acc_rot(:,ff,nf,:)=...
%             [sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,1,:,ff))',2) ...
%             sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,2,:,ff))',2) ...
%             sum(squeeze(acc(:,ff,nf,:)).*squeeze(turnmatr(:,3,:,ff))',2)];
%         vel_trans(:,ff,nf,2)=vel_rot(:,ff,nf,2)-vel_rot(:,ff,ff,2);
%         %         end
%     end
% end
% focal_to_neighbour_vector_rotated=foc_to_nb_vec_rot;
% focal_a_rotated=acc_rot;
% % focal_v_rotated=vel_rot;


invertY = true;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(tr);
[focal_to_neighbour_vector_rotated,~,focal_a_rotated]=idSocial_auxiliaries_rotateToFocalSystem(tr,vel,acc);

%%
dim=2;
no_bins=size(edges,2)-1;
no_fish=size(focal_to_neighbour_vector_rotated,2);

% turningReactOneBin=NaN(no_fish,no_fish,no_frames);
turningReact=NaN(no_fish,no_fish,no_bins);
turningReactInBins=cell(no_fish,no_fish,no_bins);
turning=NaN(no_fish,no_fish,no_frames);

for f1=1:no_fish
    for f2=1:no_fish
        if f1~=f2 && ~rand_check(f1,f2)
            
            distance_vector=focal_to_neighbour_vector_rotated(:,f2,f1,dim)'/bodylength;
            map_vector=focal_a_rotated(:,f1,f1,dim)*framerate^2/bodylength;
            %             figure; hist(map_vector(:),200)
            %             map_vector(abs(map_vector)<50)=NaN;
            %             disp('del this!')
            
            
          
            if normalize_by_individual_max
                map_vector=map_vector/max(abs(map_vector));
            end
%             turning(f1,f2,:)= map_vector;%>0;
            
            %%%%
            %             spacing=1;
            if strcmp(method,'gauss')
                ctrs=edges(1:end-1)+(edges(2)-edges(1))/2;
                prob=exp(-.5*((repmat(distance_vector',[1 numel(ctrs)])-repmat(ctrs,[no_frames 1]))/spacing).^2);
                A = prob.*repmat(map_vector,[1 numel(ctrs)]);
                
                turningReact=nanmean(A,1);
                
                turningReactInBins(f1,f2,:)=arrayfun(@(k) A(~isnan(A(:,k)),k),1:numel(ctrs),'UniformOutput',false);
                
                %             figure; plot(nansum(A,1))
                %%%%
            else
                try
                    [~, bin]=histc(distance_vector,edges);
                catch
                    keyboard
                end
                
                good_idx=bin~=0 & bin<size(edges,2);
                
                lag=spacing;
                
                A=accumarray(bin(good_idx)',map_vector(good_idx),[size(edges,2)-1 1],@nanmean);
                
                
                A=smooth(A,'moving',lag);
                
                turningReact(f1,f2,ceil(lag/2):end-floor(lag/2))=A(ceil(lag/2):end-floor(lag/2));
                
                B=accumarray(bin(good_idx)',map_vector(good_idx),[size(edges,2)-1 1] ,@(x) {x},{NaN});
                
                
                turningReactInBins(f1,f2,ceil(lag/2):end-floor(lag/2))=B(ceil(lag/2):end-floor(lag/2));
            end
            %             keyboard
            %%%%%
%             ttt = map_vector.*distance_vector';
%             ttt = ttt(~isnan(ttt));
            turning(f1,f2,:)= map_vector;
        end
    end
end

turningReactCell=num2cell(turningReact);