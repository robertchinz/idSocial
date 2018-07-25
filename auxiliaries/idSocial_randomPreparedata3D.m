function [rand_movementdata rand_trayectorias]=idSocial_randomPreparedata3D(trayectorias,options, roi)


%% Default options

act_method=mfilename;
def_options.act_method=act_method(10:end);%(10:end);
def_options.project_path='';
def_options.interpolat_trajectories=false;
def_options.seq_length=90;
def_options.rand_mode='shuffle_frames';%'shuffle_sequences'; %'shuffle_frames'
def_options.smooth_method='moving';
def_options.smooth_degree=1;
def_options.interpolation_mode='spline';
def_options.speedlim_bl_per_s=0;
def_options.maxspeed_bl_per_s=Inf;
def_options.movementdata_newformat=1;
def_options.blpxl=[];
def_options.interpolate_trajectories=false;

input_data(1,1).options=options;
[~, options]=idSocial_readparams(input_data,def_options,def_options.act_method);

%%
rand_mode=  options.rand_mode;
seq_length= options.seq_length;

%%
no_fish=size(trayectorias,2);
no_frames=size(trayectorias,1);

%%
disp([act_method ': ' rand_mode])
switch rand_mode
    case 'shuffle_sequences'
        options.interpolate_trajectories=false;
        no_seqs=floor(no_frames/seq_length);
        rand_trayectorias=NaN(no_frames,no_fish,2);
        rand_movementdata=cell(no_fish,no_fish);
        for focal=1:no_fish
            posible_nbs=setxor(focal,1:no_fish);
            
            for nb=1:no_fish
                rand_trayectorias(:,focal,focal,:)=trayectorias(:,focal,:);
                if nb~=focal
                    
                    for seq=1:no_seqs
                        seq_idx=round(rand*(no_frames-seq_length));
                        rand_nb=posible_nbs(randi(no_fish-1));
                        rand_trayectorias((seq-1)*seq_length+1,focal,nb,:)=NaN;
                        rand_trayectorias((seq-1)*seq_length+2:seq*seq_length-1,focal,nb,:)=trayectorias(seq_idx+1:seq_idx+seq_length-2,rand_nb,:);
                        rand_trayectorias(seq*seq_length,focal,nb,:)=NaN;
                    end
                    
                    temp_tray=NaN(no_frames,2,2);
                    temp_tray(:,1,:)=squeeze(rand_trayectorias(:,focal,focal,:));
                    temp_tray(:,2,:)=squeeze(rand_trayectorias(:,focal,nb,:));
                    
                    rand_movementdata{nb,focal}=temp_md{2,1};
                    
                    
                end
            end
        end
%         temp_md=...
%             idSocial_preparedata3D(rand_trayectorias,options);
%         
%         rand_movementdata{focal,nb}=temp_md{1,2};
%         
        
    case 'shuffle_frames'
        options.interpolate_trajectories=false;
        options.interpolation_mode=[];
        rand_trayectorias=NaN(no_frames,no_fish,2);
%         rand_movementdata=cell(no_fish,no_fish);
        for focal=1:no_fish
            
%                 rand_trayectorias(:,focal,:)=trayectorias(:,focal,:);
                rand_trayectorias(:,focal,:)=trayectorias(randperm(no_frames),focal,:);
               
            
        end
        rand_movementdata=idSocial_preparedata3D(rand_trayectorias,options);
        
    case 'fixed_point_center'
        options.interpolate_trajectories=false;
        options.interpolation_mode=[];
        disp('ATTENTION: Filter set to "off" for Rand mode fixed_point_center!')
        options.filters_on=false;
        if nargin>2 && ~isempty(roi) 
            center=[(roi(2,1)-roi(1,1))/2 (roi(2,2)-roi(1,2))/2];
            rand_trayectorias=NaN(no_frames,no_fish,no_fish,2);
            rand_movementdata=cell(no_fish,no_fish);
            for focal=1:no_fish
                for nb=1:no_fish
                    rand_trayectorias(:,focal,focal,:)=trayectorias(:,focal,:);
                    if nb~=focal
                        
                        rand_trayectorias(:,focal,nb,1)=ones(no_frames,1)*center(1);
                        rand_trayectorias(:,focal,nb,2)=ones(no_frames,1)*center(2);
                        temp_tray=NaN(no_frames,2,2);
                        temp_tray(:,1,:)=squeeze(rand_trayectorias(:,focal,focal,:));
                        temp_tray(:,2,:)=squeeze(rand_trayectorias(:,focal,nb,:));
                        temp_md=...
                            idSocial_preparedata3D(temp_tray,options);
                        
                        rand_movementdata{focal,nb}=temp_md{1,2};
                        %                     rand_movementdata{nb,focal}=temp_md{2,1};
                    end
                end
            end
        else
            disp('Rand mode fixed_point_center input info is needed.')
        end
        
end
%
% idx=2000:3000;
% figure; plot(temp_tray(idx,1,1),temp_tray(idx,1,2),'g-');
% hold on
% plot(temp_tray(idx,2,1),temp_tray(idx,2,2),'r.');
% figure; plot(temp_md{1,2}.focal_a_rotated(1,idx),'.-')
% hold on
% % plot(temp_md{1,2}.focal_a_rotated(1,idx),'.-')

