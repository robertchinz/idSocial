function [turnVsOrient, turnVsOrientCell]=idSocial_turningVsOrientationAtBorder(tr,edges,framerate,bodylength,arena_center,arena_radius,border_radius_normalized)
% Calculates distance distributions

if nargin<2 || isempty(edges)
    edges = 0:.05:1;
end

if nargin<3 || isempty(framerate)
    framerate=1;
end
if nargin<4 || isempty(bodylength)
    bodylength=1;
end
if nargin<5 || isempty(arena_center)
    arena_center=[0 0];
end
if nargin<6 || isempty(arena_radius)
    arena_radius=1;
end
if nargin<7 || isempty(border_radius_normalized)
    border_radius_normalized=.7;
end

dish_R_pxl=arena_radius;

% Format
invertY = true;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(tr);


vexp = vertcat(diff(tr,1,1),NaN(1,no_fish,no_fish,no_dim));
vexpMagn = sqrt(nansum(vexp.^2,4));
vexpMagn(all(isnan(tr),4))=NaN;
vexpNorm = vexp./repmat(vexpMagn,[1,1,1,size(vexp,4)]);
yax = [0 1];
heading = -atan2(vexpNorm(1:end-1,2).*yax(1)-vexpNorm(1:end-1,1).*yax(2), ...
    vexpNorm(1:end-1,1).*yax(1)+vexpNorm(1:end-1,2).*yax(2));

heading_vector = NaN(no_frames,no_fish,no_fish,no_dim);
heading_vector(1:no_frames-1,:,:,1:2) =cat(4,cos(heading),sin(heading));

vec_to_com =  repmat(dish_R_pxl,[no_frames,no_fish,no_fish,no_dim])-tr;
vec_to_comNorm = vec_to_com./repmat(sqrt(sum(vec_to_com.^2,4)),[1 1 1 size(vec_to_com,4)]);

heading_to_com = acos(nansum(vec_to_comNorm.*heading_vector,4));

turn = NaN(no_frames,no_fish,no_fish);
turn(1:no_frames-1,:,:) = abs(atan2(vexpNorm(1:end-1,:,:,2).*vexpNorm(2:end,:,:,1)-vexpNorm(1:end-1,:,:,1).*vexpNorm(2:end,:,:,2), ...
                vexpNorm(1:end-1,:,:,1).*vexpNorm(2:end,:,:,1)+vexpNorm(1:end-1,:,:,2).*vexpNorm(2:end,:,:,2)));

turnVsOrient = NaN(no_fish,no_fish,numel(edges)-1);
turnVsOrientCell = cell(no_fish,no_fish,numel(edges)-1);


for ff = 1:no_fish
    for nf = 1:no_fish
            val = turn(:,ff,nf);
            [~,bins]=histc(squeeze(heading_to_com(:,ff,nf)),edges);
            turnVsOrientCell(ff,nf,:)=accumarray(bins(bins>0),val(bins>0),[numel(edges)-1 1],@(x) {x},{NaN});
            turnVsOrient = cellfun(@(x) nanmean(x),turnVsOrientCell(ff,nf,:));
%         end
    end
end



