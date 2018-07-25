function [pos_ratio, meanTurn2Ctr, meanTurn, turningCell, turningAngleCell, turningAngleTowardsCell, turningAngleAwayCell, speedCell]=idSocial_turningToArenaCenter(tr,edges,framerate,bodylength,arena_center,arena_radius)
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


% Format
invertY = true;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

if no_fish==1
    no_fishOrig=1;
    no_fish=2;
    tr_temp = NaN(no_frames,2,2,no_dim);
    tr_temp(:,1,1,:) = tr;
    tr_temp(:,1,2,:) = tr;
    tr_temp(:,2,1,1:2) = repmat(arena_center'.*[1 -1],[no_frames,1, 1,1]);
    tr_temp(:,2,2,1:2) = repmat(arena_center'.*[1 -1],[no_frames,1, 1,1]);
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
                tr(:,nf,ff,1:2) = repmat(arena_center'.*[1 -1],[no_frames,1, 1,1]);
            end
        end
    end
end
% Re-calculate vel and acc
invertYReCalc = false;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertYReCalc);


% Rotate
[focal_to_neighbour_vector_rotated,~,focal_a_rotated]=idSocial_auxiliaries_rotateToFocalSystem(tr,vel,acc);

dim = 1; %Turning = 1;
turning = cell(no_fish,no_fish);
trnangle = cell(no_fish,no_fish);
trnangleTowards = cell(no_fish,no_fish);
trnangleAway = cell(no_fish,no_fish);
speed = cell(no_fish,no_fish);
turningCell = cell(no_fish,no_fish,numel(edges)-1);
speedCell= cell(no_fish,no_fish,numel(edges)-1);
turningAngleCell = cell(no_fish,no_fish,numel(edges)-1);
turningAngleTowardsCell = cell(no_fish,no_fish,numel(edges)-1);
turningAngleAwayCell = cell(no_fish,no_fish,numel(edges)-1);
pos_ratio = NaN(no_fish,no_fish,numel(edges)-1);
meanTurn2Ctr = NaN(no_fish,no_fish,numel(edges)-1);
meanTurn= NaN(no_fish,no_fish,numel(edges)-1);
a = NaN(no_frames,1);
for ff = 1:no_fish
    for nf = 1:no_fish
        if ff ~=nf
            % Turning distribution
            distance_vector = sign(focal_to_neighbour_vector_rotated(:,nf,ff,dim));
            val=focal_a_rotated(:,ff,ff,dim);
            ttt = val.*sign(distance_vector);
            turning{ff,nf}= ttt;
            
            vexp = vertcat(diff(squeeze(tr(:,ff,ff,:)),1,1),NaN(1,3));
            vexpMagn = sqrt(nansum(vexp.^2,2));
            vexpNorm = vexp./repmat(vexpMagn,[1,size(vexp,2)]);
            a(1:no_frames-1,:) = atan2(vexpNorm(1:end-1,2).*vexpNorm(2:end,1)-vexpNorm(1:end-1,1).*vexpNorm(2:end,2), ...
                vexpNorm(1:end-1,1).*vexpNorm(2:end,1)+vexpNorm(1:end-1,2).*vexpNorm(2:end,2));
            trnangle{ff,nf} = a.*sign(distance_vector);
            trnangleTowards{ff,nf} = trnangle{ff,nf};
            trnangleTowards{ff,nf}(trnangle{ff,nf}<0) = NaN;
            trnangleAway{ff,nf} = trnangle{ff,nf};
            trnangleAway{ff,nf}(trnangle{ff,nf}>0) = NaN;
            speed{ff,nf} = vexpMagn;
        end
    end
end



% Distance from center:
ctr = arena_center;
dist2ctr = NaN(no_frames,no_fish);

if invertY
    ctr(2)=ctr(2)*nanmean(sign(tr(:,1,1,2)));
end
for ff = 1:no_fish
    dist2ctr(:,ff) = sqrt(sum((squeeze(tr(:,ff,ff,1:2))-repmat(ctr',[no_frames 1])).^2,2))/arena_radius;
    
end

for ff = 1:no_fish
    for nf = 1:no_fish
        if ff ~=nf
            [hi,bins]=histc(dist2ctr(:,ff),edges);
            turningCell(ff,nf,:)=accumarray(bins(bins>0),turning{ff,nf}(bins>0),[numel(edges)-1 1],@(x) {x},{NaN});
            turningAngleCell(ff,nf,:)=accumarray(bins(bins>0),trnangle{ff,nf}(bins>0),[numel(edges)-1 1],@(x) {x},{NaN});
            turningAngleTowardsCell(ff,nf,:)=accumarray(bins(bins>0),trnangleTowards{ff,nf}(bins>0),[numel(edges)-1 1],@(x) {x},{NaN});
            turningAngleAwayCell(ff,nf,:)=accumarray(bins(bins>0),trnangleAway{ff,nf}(bins>0),[numel(edges)-1 1],@(x) {x},{NaN});
            speedCell(ff,nf,:)=accumarray(bins(bins>0),speed{ff,nf}(bins>0),[numel(edges)-1 1],@(x) {x},{NaN});
            
            for bn = 1:numel(edges)-1
                
                pos_ratio(ff,nf,bn) = idSocial_positive_ratio(turningCell{ff,nf,bn});
                meanTurn2Ctr(ff,nf,bn) = nanmean(turningCell{ff,nf,bn});
                meanTurn(ff,nf,bn) = nanmean(turningAngleCell{ff,nf,bn});
                
            end
        end % if ff~=nf
    end %nf
end %ff
if no_fishOrig==1
    turningCell = turningCell(1,2,:);
    turningAngleCell = turningAngleCell(1,2,:);
    turningAngleTowardsCell = turningAngleTowardsCell(1,2,:);
    turningAngleAwayCell = turningAngleAwayCell(1,2,:);
    speedCell = speedCell(1,2,:);
    pos_ratio = pos_ratio(1,2,:);
    meanTurn2Ctr = meanTurn2Ctr(1,2,:);
    meanTurn = meanTurn(1,2,:);
end


