function [frontbackhist, funcinfo] = idSocial_scoreFrontBack(trajectory,maxdist,framerate,bodylength)


if nargin<1 || isempty(trajectory)
    error([mfilename ': Trajectory missing.'])
end

if nargin<2 || isempty(maxdist)
    maxdist = 6;
end

invertY = true;
[trajectory,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(trajectory,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(trajectory);


foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,no_dim);
vel_magn=sqrt(sum(vel.^2,4));
vel_norm=bsxfun(@rdivide,vel,vel_magn);

for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,nf,ff,:)=squeeze(trajectory(:,nf,ff,:)-trajectory(:,ff,ff,:));
        end
    end
end

foc_to_nb_vec_rot = ...
    idSocial_auxiliaries_rotateToFocalSystem(trajectory,vel,acc);


foc_to_nb_distance = sqrt(sum(foc_to_nb_vec.^2,4));
dist_filter = foc_to_nb_distance/bodylength > maxdist;
foc_to_nb_vec_rot = foc_to_nb_vec_rot(:,:,:,2);
%%
foc_to_nb_vec_rot(dist_filter) = NaN;
foc_to_nb_vec_rot = permute(foc_to_nb_vec_rot,[2 3 1]);

frontbackhist=NaN(no_fish,no_fish,no_frames);
frontbackhist(foc_to_nb_vec_rot>0) = 1;
frontbackhist(foc_to_nb_vec_rot<0) = -1;


frontbackhist = nansum(frontbackhist,3);
%%
fstack=dbstack(3);
if ~isempty(fstack)
    funcinfo.Function = mfilename;
    funcinfo.callerFunction = fstack(1).name;
    funcinfo.no_fish = no_fish;
end
% funcinfo.XTick = edges(1:end-1);
% funcinfo.XTickLabel = cellfun(@(x) strtrim(x),cellstr((num2str(edges(1:end-1)')))','UniformOutput',false);