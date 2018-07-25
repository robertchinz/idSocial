function [speedVsTurning, speedVsTurningCell]=idSocial_speedVsTurning(tr,edges,framerate,bodylength)
% Calculates distance distributions

if nargin<2 || isempty(edges)
    edges = 0:.05:3;
end

if nargin<3 || isempty(framerate)
    framerate=1;
end
if nargin<4 || isempty(bodylength)
    bodylength=1;
end


% Format
invertY = true;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(tr);






speedVsTurning = NaN(no_fish,no_fish,numel(edges)-1);
speedVsTurningCell= cell(no_fish,no_fish,numel(edges)-1);

a = NaN(no_frames,1);
for ff = 1:no_fish
    for nf = 1:no_fish
%         if ff ~=nf
            % Turning distribution
            vexp = vertcat(diff(squeeze(tr(:,ff,ff,:)),1,1),NaN(1,3));
            vexpMagn = sqrt(nansum(vexp.^2,2));
            vexpMagn(all(isnan(squeeze(tr(:,ff,ff,:))),2))=NaN;
            vexpNorm = vexp./repmat(vexpMagn,[1,size(vexp,2)]);
            a(1:no_frames-1,:) = abs(atan2(vexpNorm(1:end-1,2).*vexpNorm(2:end,1)-vexpNorm(1:end-1,1).*vexpNorm(2:end,2), ...
                vexpNorm(1:end-1,1).*vexpNorm(2:end,1)+vexpNorm(1:end-1,2).*vexpNorm(2:end,2)));
            [~,bins]=histc(a,edges);
            speedVsTurningCell(ff,nf,:)=accumarray(bins(bins>0),vexpMagn(bins>0),[numel(edges)-1 1],@(x) {x},{NaN});
            speedVsTurning = cellfun(@(x) nanmean(x),speedVsTurningCell(ff,nf,:));
%         end
    end
end




