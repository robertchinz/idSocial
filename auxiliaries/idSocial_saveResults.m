% function input_data=idSocial_saveResults(input_data,act_method,raw)
function input_data=idSocial_saveResults(input_data,act_method, output, output_rand)

%% Calculate some means
if ~isempty(output)
    
    sz=size(output);
    ndm=ndims(output);
    mean_groups=nanmean(output,1);
    mean_groups_std=nanstd(output,[],1);
    mean_trials=nanmean(output,2);
    mean_trials_std=nanstd(output,[],2);
    mean_parts=nanmean(output,3);
    mean_parts_std=nanstd(output,[],3);
    mean_foc=nanmean(output,4);
    mean_foc_std=nanstd(output,[],4);
    mean_nbs=nanmean(output,5);
    mean_nbs_std=nanstd(output,[],5);
    
    mean_nbs_prts=nanmean(mean_nbs(:,:,2:end-1,:,:,:,:,:,:),3);
    mean_nbs_prts_std=nanstd(mean_nbs(:,:,2:end-1,:,:,:,:,:,:),[],3);
    mean_nbs_trials=nanmean(mean_nbs,2);
    mean_nbs_trials_std=nanstd(mean_nbs,[],2);
    mean_parts_trials=nanmean(mean_parts,2);
    mean_parts_trials_std=nanstd(mean_parts,[],2);
    mean_nbs_foc=nanmean(mean_nbs,4);
    mean_nbs_foc_std=nanstd(mean_nbs,[],4);
    mean_trials_nbs=nanmean(mean_trials,5);
    mean_trials_nbs_std=nanstd(mean_trials,[],5);
    
%     clps_nbs_trials=reshape(permute(output,[setxor(1:ndm,[2 5]) [2 5]]),[sz(setxor(1:ndm,[2 5])) sz(2)*sz(5)]);
    
    mean_nbs_prts_trials=nanmean(mean_nbs_prts,2);
    mean_nbs_prts_trials_std=nanstd(mean_nbs_prts,[],2);
    mean_nbs_foc_prts=nanmean(mean_nbs_foc,3);
    mean_nbs_foc_prts_std=nanstd(mean_nbs_foc,[],3);
    mean_nbs_foc_trials=nanmean(mean_nbs_foc,2);
    mean_nbs_foc_trials_std=nanstd(mean_nbs_foc,[],2);
    mean_trials_nbs_foc=nanmean(mean_trials_nbs,4);
    mean_trials_nbs_foc_std=nanstd(mean_trials_nbs,[],4);
    
    mean_nbs_foc_prts_trials=nanmean(mean_nbs_foc_prts,2);
    mean_nbs_foc_prts_trials_std=nanstd(mean_nbs_foc_prts,[],2);
    %% Put everything into a structure
    input_data(1,1).(act_method).output=output;
    input_data(1,1).(act_method).output_std=zeros(size(output));
    input_data(1,1).(act_method).mean_groups=mean_groups;
    input_data(1,1).(act_method).mean_groups_std=mean_groups_std;
    input_data(1,1).(act_method).mean_trials=mean_trials;
    input_data(1,1).(act_method).mean_trials_std=mean_trials_std;
    input_data(1,1).(act_method).mean_prts=mean_parts;
    input_data(1,1).(act_method).mean_prts_std=mean_parts_std;
    input_data(1,1).(act_method).mean_foc=mean_foc;
    input_data(1,1).(act_method).mean_foc_std=mean_foc_std;
    input_data(1,1).(act_method).mean_nbs=mean_nbs;
    input_data(1,1).(act_method).mean_nbs_std=mean_nbs_std;
    
    input_data(1,1).(act_method).mean_nbs_prts=mean_nbs_prts;
    input_data(1,1).(act_method).mean_nbs_prts_std=mean_nbs_prts_std;
    input_data(1,1).(act_method).mean_prts_trials=mean_parts_trials;
    input_data(1,1).(act_method).mean_prts_trials_std=mean_parts_trials_std;
    input_data(1,1).(act_method).mean_nbs_trials=mean_nbs_trials;
    input_data(1,1).(act_method).mean_nbs_trials_std=mean_nbs_trials_std;
    input_data(1,1).(act_method).mean_nbs_foc=mean_nbs_foc;
    input_data(1,1).(act_method).mean_nbs_foc_std=mean_nbs_foc_std;
    input_data(1,1).(act_method).mean_trials_nbs=mean_trials_nbs;
    input_data(1,1).(act_method).mean_trials_nbs_std=mean_trials_nbs_std;
    input_data(1,1).(act_method).mean_trials_nbs_foc=mean_trials_nbs_foc;
    input_data(1,1).(act_method).mean_trials_nbs_foc_std=mean_trials_nbs_foc_std;
    
    input_data(1,1).(act_method).mean_nbs_prts_trials=mean_nbs_prts_trials;
    input_data(1,1).(act_method).mean_nbs_prts_trials_std=mean_nbs_prts_trials_std;
    input_data(1,1).(act_method).mean_nbs_foc_trials=mean_nbs_foc_trials;
    input_data(1,1).(act_method).mean_nbs_foc_trials_std=mean_nbs_foc_trials_std;
    
    input_data(1,1).(act_method).mean_nbs_foc_prts=mean_nbs_foc_prts;
    input_data(1,1).(act_method).mean_nbs_foc_prts_std=mean_nbs_foc_prts_std;
    input_data(1,1).(act_method).mean_nbs_foc_prts_trials=mean_nbs_foc_prts_trials;
    input_data(1,1).(act_method).mean_nbs_foc_prts_trials_std=mean_nbs_foc_prts_trials_std;
end

if ~isempty(output_rand)
   rand_mean_groups=nanmean(output_rand,1);
    rand_mean_groups_std=nanstd(output_rand,[],1);
    rand_mean_trials=nanmean(output_rand,2);
    rand_mean_trials_std=nanstd(output_rand,[],2);
    rand_mean_parts=nanmean(output_rand,3);
    rand_mean_parts_std=nanstd(output_rand,[],3);
    rand_mean_foc=nanmean(output_rand,4);
    rand_mean_foc_std=nanstd(output_rand,[],4);
    rand_mean_nbs=nanmean(output_rand,5);
    rand_mean_nbs_std=nanstd(output_rand,[],5);
    
    rand_mean_nbs_prts=nanmean(rand_mean_nbs(:,:,2:end-1,:,:,:,:,:,:),3);
    rand_mean_nbs_prts_std=nanstd(rand_mean_nbs(:,:,2:end-1,:,:,:,:,:,:),[],3);
    rand_mean_nbs_trials=nanmean(rand_mean_nbs,2);
    rand_mean_nbs_trials_std=nanstd(rand_mean_nbs,[],2);
    rand_mean_nbs_foc=nanmean(rand_mean_nbs,4);
    rand_mean_nbs_foc_std=nanstd(rand_mean_nbs,[],4);
    rand_mean_trials_nbs=nanmean(rand_mean_trials,5);
    rand_mean_trials_nbs_std=nanstd(rand_mean_trials,[],5);
    
    rand_mean_nbs_prts_trials=nanmean(rand_mean_nbs_prts,2);
    rand_mean_nbs_prts_trials_std=nanstd(rand_mean_nbs_prts,[],2);
    rand_mean_nbs_foc_prts=nanmean(rand_mean_nbs_foc,3);
    rand_mean_nbs_foc_prts_std=nanstd(rand_mean_nbs_foc,[],3);
    rand_mean_nbs_foc_trials=nanmean(rand_mean_nbs_foc,2);
    rand_mean_nbs_foc_trials_std=nanstd(rand_mean_nbs_foc,[],2);
    
    rand_mean_nbs_foc_prts_trials=nanmean(rand_mean_nbs_foc_prts,2);
    rand_mean_nbs_foc_prts_trials_std=nanstd(rand_mean_nbs_foc_prts,[],2);
    %% Put everything into a structure
    input_data(1,1).(act_method).output_rand=output_rand;
    input_data(1,1).(act_method).rand_mean_groups=rand_mean_groups;
    input_data(1,1).(act_method).rand_mean_groups_std=rand_mean_groups_std;
    input_data(1,1).(act_method).rand_mean_trials=rand_mean_trials;
    input_data(1,1).(act_method).rand_mean_trials_std=rand_mean_trials_std;
    input_data(1,1).(act_method).rand_mean_parts=rand_mean_parts;
    input_data(1,1).(act_method).rand_mean_parts_std=rand_mean_parts_std;
    input_data(1,1).(act_method).rand_mean_foc=rand_mean_foc;
    input_data(1,1).(act_method).rand_mean_foc_std=rand_mean_foc_std;
    input_data(1,1).(act_method).rand_mean_nbs=rand_mean_nbs;
    input_data(1,1).(act_method).rand_mean_nbs_std=rand_mean_nbs_std;
    
    input_data(1,1).(act_method).rand_mean_nbs_prts=rand_mean_nbs_prts;
    input_data(1,1).(act_method).rand_mean_nbs_prts_std=rand_mean_nbs_prts_std;
    input_data(1,1).(act_method).rand_mean_nbs_trials=rand_mean_nbs_trials;
    input_data(1,1).(act_method).rand_mean_nbs_trials_std=rand_mean_nbs_trials_std;
    input_data(1,1).(act_method).rand_mean_nbs_foc=rand_mean_nbs_foc;
    input_data(1,1).(act_method).rand_mean_nbs_foc_std=rand_mean_nbs_foc_std;
    input_data(1,1).(act_method).rand_mean_trials_nbs=rand_mean_trials_nbs;
    input_data(1,1).(act_method).rand_mean_trials_nbs_std=rand_mean_trials_nbs_std;
    
    input_data(1,1).(act_method).rand_mean_nbs_prts_trials=rand_mean_nbs_prts_trials;
    input_data(1,1).(act_method).rand_mean_nbs_prts_trials_std=rand_mean_nbs_prts_trials_std;
    input_data(1,1).(act_method).rand_mean_nbs_foc_trials=rand_mean_nbs_foc_trials;
    input_data(1,1).(act_method).rand_mean_nbs_foc_trials_std=rand_mean_nbs_foc_trials_std;
    
    input_data(1,1).(act_method).rand_mean_nbs_foc_prts=rand_mean_nbs_foc_prts;
    input_data(1,1).(act_method).rand_mean_nbs_foc_prts_std=rand_mean_nbs_foc_prts_std;
    input_data(1,1).(act_method).rand_mean_nbs_foc_prts_trials=rand_mean_nbs_foc_prts_trials;
    input_data(1,1).(act_method).rand_mean_nbs_foc_prts_trials_std=rand_mean_nbs_foc_prts_trials_std;
end
