function [rand_trayectorias,options]=idSocial_randomPrepareTrajectories3D(trayectorias,options,frame_range,structOut)

if nargin < 2 || isempty(options)
    options = [];
end
if nargin < 3 || isempty(frame_range)
    frame_range = [];
end
if nargin < 4 || isempty(structOut)
    structOut = false;
end

%% Default options

act_method=mfilename;
def_options.act_method=act_method(10:end);
def_options.project_path='';
def_options.interpolat_trajectories=false;
def_options.random_data_seq_length=60;
def_options.random_data_method='shuffle_frames';
def_options.smooth_method='moving';
def_options.smooth_degree=1;
def_options.interpolation_mode='spline';
def_options.speedlim_bl_per_s=0;
def_options.maxspeed_bl_per_s=Inf;
def_options.movementdata_newformat=1;
def_options.blpxl=[];
def_options.interpolate_trajectories=false;
def_options.min_idx_dist=min(12000,round(size(trayectorias,1)/8));
def_options.no_RandomNeighbors = [];

input_data(1,1).options=options;
[~, options]=idSocial_readparams(input_data,options,def_options,def_options.act_method);

%%
rand_mode=  options.random_data_method;
seq_length= options.random_data_seq_length;
mindist=options.min_idx_dist;
%%


iscellTr = false;
if iscell(trayectorias) 
    iscellTr = true;
    tr_orig = trayectorias;
	trayectorias = idSocial_auxiliaries_formatInputTrajectory(trayectorias);
end
if isstruct(trayectorias)
    no_fish=size(trayectorias.Tr,2);
    no_frames=size(trayectorias.Tr,1);
else
    no_fish=size(trayectorias,2);
    no_frames=size(trayectorias,1);
end


if ~isempty(options.no_RandomNeighbors) && isnumeric(options.no_RandomNeighbors) && options.no_RandomNeighbors > 1
    no_RandomNeighbors = options.no_RandomNeighbors;
else
    no_RandomNeighbors = no_fish;
end

% if ~isempty(options.filter_minProbIdentityAssignment) && filter_probs
%     trayectorias(...
%         repmat(probtrayectorias,...
%         [1 1 size(tr,3)])<...
%         options.filter_minProbIdentityAssignment) = NaN;
% end
%%
% disp([act_method ': ' rand_mode])

switch rand_mode
    case 'shuffle_sequences'
        options.interpolate_trajectories=false;
        rand_trayectorias = ...
            idSocial_randomTrajectoriesShuffleSequences(trayectorias,no_RandomNeighbors,mindist,seq_length,frame_range,structOut);
%         no_seqs=floor(no_frames/seq_length);
%         rand_trayectorias=NaN(size(trayectorias,1),no_RandomNeighbors,no_RandomNeighbors,size(trayectorias,4));
% 
%         %         rand_movementdata=cell(no_fish,no_fish);
%         for focal=1:no_RandomNeighbors
% %             posible_nbs=setxor(focal,1:no_fish);
%             
%             for nb=1:no_RandomNeighbors
%                 if focal>no_fish && nb>no_fish
%                     rand_trayectorias(:,focal,nb,1)=NaN(1,no_frames);
%                     rand_trayectorias(:,focal,nb,2)=Inf(1,no_frames);
%                 else
%                     
%                     if nb~=focal
%                         
%                         for seq=1:no_seqs
%                             seq_idx=   (1:no_frames-seq_length)+randi([mindist,no_frames-mindist+1]);%     round(rand*(no_frames-seq_length));
%                             seq_idx(seq_idx>no_frames)=seq_idx(seq_idx>no_frames)-no_frames;
% %                             rand_nb=posible_nbs(randi(no_fish-1));
%                             rand_trayectorias((seq-1)*seq_length+1,focal,nb,:)=NaN;
%                             rand_trayectorias((seq-1)*seq_length+2:seq*seq_length-1,focal,nb,:)=trayectorias(seq_idx+1:min(seq_idx+seq_length-2,no_frames),mod(focal-1,no_fish)+1,mod(nb-1,no_fish)+1,:);
%                             rand_trayectorias(seq*seq_length,focal,nb,:)=NaN;
%                         end
%                     else
%                         rand_trayectorias(:,focal,focal,:)=trayectorias(:,mod(focal-1,no_fish)+1,mod(focal-1,no_fish)+1,:);
%                         
%                     end
%                 end
%             end
%         end
        
    case 'shuffle_frames'
        
        rand_trayectorias = idSocial_randomTrajectoriesShuffleFrames(trayectorias,no_RandomNeighbors,mindist,frame_range,structOut);
        
        
%     case 'fixed_point_center'
%         options.interpolate_trajectories=false;
%         options.interpolation_mode=[];
%         disp('ATTENTION: Filter set to "off" for Rand mode fixed_point_center!')
%         options.filters_on=false;
%         if nargin>2 && ~isempty(roi)
%             center=[(roi(2,1)-roi(1,1))/2 (roi(2,2)-roi(1,2))/2];
%             rand_trayectorias=NaN(no_frames,no_fish,no_fish,2);
%             rand_movementdata=cell(no_fish,no_fish);
%             for focal=1:no_fish
%                 for nb=1:no_RandomNeighbors
%                     rand_trayectorias(:,focal,focal,:)=trayectorias(:,focal,:);
%                     if nb~=focal
%                         
%                         rand_trayectorias(:,focal,nb,1)=ones(no_frames,1)*center(1);
%                         rand_trayectorias(:,focal,nb,2)=ones(no_frames,1)*center(2);
%                         temp_tray=NaN(no_frames,2,2);
%                         temp_tray(:,1,:)=squeeze(rand_trayectorias(:,focal,focal,:));
%                         temp_tray(:,2,:)=squeeze(rand_trayectorias(:,focal,nb,:));
%                         %                         temp_md=...
%                         %                             idSocial_preparedata3D(temp_tray,options);
%                         
%                         %                         rand_movementdata{focal,nb}=temp_md{1,2};
%                         %                     rand_movementdata{nb,focal}=temp_md{2,1};
%                     end
%                 end
%             end
%         else
%             disp('Rand mode fixed_point_center input info is needed.')
%         end
        
end

if iscellTr
    tr_orig.Tr = rand_trayectorias;
    rand_trayectorias = tr_orig;
end

