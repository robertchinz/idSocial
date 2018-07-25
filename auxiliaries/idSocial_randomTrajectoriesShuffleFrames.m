function outtr = idSocial_randomTrajectoriesShuffleFrames(trayectorias,no_RandomNeighbors,mindist,frame_range,structOut)
% frame_range: if no_frames > no_framesOut, it is not clear which
% no_framesOut-Frames to take for focal,focal (=original trajectories in
% the diagonal). This can be set by frame_range.
if nargin<5 || isempty(structOut)
    structOut = false;
end

struct_flag = false;
if isstruct(trayectorias) && isfield(trayectorias,'Tr') && ...
        isfield(trayectorias,'Vel') &&  isfield(trayectorias,'Acc')
    vel = trayectorias.Vel;
    acc = trayectorias.Acc;
    trayectorias = trayectorias.Tr;
    struct_flag = true;
end
no_fish = size(trayectorias,2);
no_frames = size(trayectorias,1);

dims = ndims(trayectorias);
if nargin <2 || isempty(no_RandomNeighbors)
    no_RandomNeighbors = no_fish;
end
if nargin <4 || isempty(frame_range)
    frame_range = 1:no_frames;
    no_framesOut = no_frames;
else
    no_framesOut = numel(frame_range);
end

no_RandomFocals = max(no_RandomNeighbors,no_fish);
no_dims = size(trayectorias,dims);

if nargin < 3 || isempty(mindist)
    mindist = min(12000,round(size(trayectorias,1)/8));
end

rand_trayectorias=NaN(no_framesOut,no_RandomFocals,no_RandomNeighbors,no_dims);

if dims == 4
    randidx=repmat((1:no_frames)',[1 no_RandomFocals*no_RandomNeighbors])+randi([mindist,no_frames-mindist+1],[no_frames,no_RandomFocals*no_RandomNeighbors]);
    randidx(randidx>no_frames)=randidx(randidx>no_frames)-no_frames;
    randidx(randidx<0)=no_frames-randidx(randidx<0)+1;
    for focal=1:no_RandomFocals
        for neighbor=1:no_RandomNeighbors
            if focal>no_fish && neighbor>no_fish
                rand_trayectorias(:,focal,neighbor,1)=NaN(1,no_frames);
                rand_trayectorias(:,focal,neighbor,2)=Inf(1,no_frames);
            else
                if focal~=neighbor
                    
                    act_focal = mod(focal-1,no_fish)+1;
                    pos_nbs = setxor(1:no_fish,act_focal);
                    act_neighbor = pos_nbs(mod(neighbor-1,numel(pos_nbs))+1);
                    rand_trayectorias(:,focal,neighbor,:)=trayectorias(randidx(:,no_RandomNeighbors*(focal-1)+neighbor),act_focal,act_neighbor,:);
                else
                    rand_trayectorias(:,focal,focal,:)=trayectorias(:,mod(focal-1,no_fish)+1,mod(focal-1,no_fish)+1,:);
                    
                end
            end
        end
        
    end
elseif dims == 3
    if struct_flag
        rand_vel=NaN(no_framesOut,no_RandomFocals,no_RandomNeighbors,no_dims);
        rand_acc=NaN(no_framesOut,no_RandomFocals,no_RandomNeighbors,no_dims);
    end

    randidx=repmat((1:no_framesOut)',[1 no_RandomFocals*no_RandomNeighbors])+randi([mindist,no_frames-mindist+1],[no_framesOut,no_RandomFocals*no_RandomNeighbors]);
    randidx(randidx>no_frames)=randidx(randidx>no_frames)-no_frames;
    randidx(randidx<0)=no_frames-randidx(randidx<0)+1;
    for focal=1:no_RandomFocals
        for neighbor=1:no_RandomNeighbors
            if focal>no_fish && neighbor>no_fish
                rand_trayectorias(:,focal,neighbor,1)=NaN(1,no_framesOut);
                rand_trayectorias(:,focal,neighbor,2)=Inf(1,no_framesOut);
            else
                if focal~=neighbor
                    
                    act_focal = mod(focal-1,no_fish)+1;
                    rand_trayectorias(:,focal,neighbor,:)=squeeze(trayectorias(randidx(:,no_RandomNeighbors*(focal-1)+neighbor),act_focal,:));
                    if struct_flag || structOut
                        rand_vel(:,focal,neighbor,:)=squeeze(vel(randidx(:,no_RandomNeighbors*(focal-1)+neighbor),act_focal,:));
                        rand_acc(:,focal,neighbor,:)=squeeze(acc(randidx(:,no_RandomNeighbors*(focal-1)+neighbor),act_focal,:));
                    end

                else
                    rand_trayectorias(:,focal,focal,:)=squeeze(trayectorias(frame_range,mod(focal-1,no_fish)+1,:));
                    if struct_flag || structOut
                        rand_vel(:,focal,focal,:)=squeeze(vel(frame_range,mod(focal-1,no_fish)+1,:));
                        rand_acc(:,focal,focal,:)=squeeze(vel(frame_range,mod(focal-1,no_fish)+1,:));
                    end
                    
                end
            end
        end
        
    end
end
if struct_flag || structOut
    outtr.Tr = rand_trayectorias;
    clear rand_trayectorias
    outtr.Vel = rand_vel;
    clear rand_vel;
    outtr.Acc = rand_acc;
    clear rand_acc;
else
    outtr = rand_trayectorias;
end