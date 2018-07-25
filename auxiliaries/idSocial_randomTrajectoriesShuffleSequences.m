function outtr = idSocial_randomTrajectoriesShuffleSequences(trayectorias,no_RandomNeighbors,mindist,seq_length,frame_range,structOut)
% frame_range: if no_frames > no_framesOut, it is not clear which
% no_framesOut-Frames to take for focal,focal (=original trajectories in
% the diagonal). This can be set by frame_range.
if nargin<6 || isempty(structOut)
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

if nargin <4 || isempty(seq_length)
    seq_length = 60;
end

if nargin <5 || isempty(frame_range)
    frame_range = 1:no_frames;
    no_framesOut = no_frames;
else
    no_framesOut = numel(frame_range);
end

no_dims = size(trayectorias,dims);

if nargin < 3 || isempty(mindist)
    mindist = min(12000,round(size(trayectorias,1)/8));
end

no_seqs=floor(no_frames/seq_length);


rand_trayectorias=NaN(size(trayectorias,1),no_RandomNeighbors,no_RandomNeighbors,no_dims);

if struct_flag
    rand_vel=NaN(no_framesOut,no_RandomNeighbors,no_RandomNeighbors,no_dims);
    rand_acc=NaN(no_framesOut,no_RandomNeighbors,no_RandomNeighbors,no_dims);
end

if dims == 4
    
    for focal=1:no_RandomNeighbors
        
        for nb=1:no_RandomNeighbors
            if focal>no_fish && nb>no_fish
                rand_trayectorias(:,focal,nb,1)=NaN(1,no_frames);
                rand_trayectorias(:,focal,nb,2)=Inf(1,no_frames);
            else
                
                if nb~=focal
                    
                    for seq=1:no_seqs
                        seq_idx=   (1:no_frames-seq_length)+randi([mindist,no_frames-mindist+1]);%     round(rand*(no_frames-seq_length));
                        seq_idx(seq_idx>no_frames)=seq_idx(seq_idx>no_frames)-no_frames;
                        %                             rand_nb=posible_nbs(randi(no_fish-1));
                        rand_trayectorias((seq-1)*seq_length+1,focal,nb,:)=NaN;
                        rand_trayectorias((seq-1)*seq_length+2:seq*seq_length-1,focal,nb,:)=trayectorias(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,mod(nb-1,no_fish)+1,:);
                        rand_trayectorias(seq*seq_length,focal,nb,:)=NaN;
                        if struct_flag || structOut
                            rand_vel((seq-1)*seq_length+2:min(seq*seq_length-1,(seq-1)*seq_length+1+numel(seqIdces)),focal,nb,:)=squeeze(vel(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,mod(nb-1,no_fish)+1,:));
                            rand_acc((seq-1)*seq_length+2:min(seq*seq_length-1,(seq-1)*seq_length+1+numel(seqIdces)),focal,nb,:)=squeeze(acc(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,mod(nb-1,no_fish)+1,:));
                        end
                    end
                else
                    rand_trayectorias(:,focal,focal,:)=trayectorias(:,mod(focal-1,no_fish)+1,mod(focal-1,no_fish)+1,:);
                    if struct_flag || structOut
                        rand_vel(:,focal,focal,:)=squeeze(vel(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,mod(focal-1,no_fish)+1,:));
                        rand_acc(:,focal,focal,:)=squeeze(acc(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,mod(focal-1,no_fish)+1,:));
                    end
                    
                end
            end
        end
    end
elseif dims == 3
    for focal=1:no_RandomNeighbors
        
        for nb=1:no_RandomNeighbors
            if focal>no_fish && nb>no_fish
                rand_trayectorias(:,focal,nb,1)=NaN(1,no_frames);
                rand_trayectorias(:,focal,nb,2)=Inf(1,no_frames);
            else
                
                if nb~=focal
                    
                    for seq=1:no_seqs
                        seq_idx=   (1:no_framesOut-seq_length)+randi([mindist,no_frames-mindist+1]);%     round(rand*(no_frames-seq_length));
                        seq_idx(seq_idx>no_frames)=seq_idx(seq_idx>no_frames)-no_frames;
                        %                             rand_nb=posible_nbs(randi(no_fish-1));
                        rand_trayectorias((seq-1)*seq_length+1,focal,nb,:)=NaN;
                        try
                            seqIdces = seq_idx+1:min(seq_idx+seq_length-2,no_frames);
                        rand_trayectorias((seq-1)*seq_length+2:min(seq*seq_length-1,(seq-1)*seq_length+1+numel(seqIdces)) ,focal,nb,:)=squeeze(trayectorias(seqIdces,mod(focal-1,no_fish)+1,:));
                        if struct_flag || structOut 
                            rand_vel((seq-1)*seq_length+2:min(seq*seq_length-1,(seq-1)*seq_length+1+numel(seqIdces)),focal,nb,:)=squeeze(vel(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,:));
                            rand_acc((seq-1)*seq_length+2:min(seq*seq_length-1,(seq-1)*seq_length+1+numel(seqIdces)),focal,nb,:)=squeeze(acc(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,:));
                        end
                        catch
                            keyboard
                        end
                        rand_trayectorias(seq*seq_length,focal,nb,:)=NaN;
                    end
                else
                    rand_trayectorias(:,focal,focal,:)=squeeze(trayectorias(frame_range,mod(focal-1,no_fish)+1,:));
                    if struct_flag || structOut
                        rand_vel(:,focal,focal,:)=squeeze(vel(frame_range,mod(focal-1,no_fish)+1,:));
                        rand_acc(:,focal,focal,:)=squeeze(acc(frame_range,mod(focal-1,no_fish)+1,:));
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