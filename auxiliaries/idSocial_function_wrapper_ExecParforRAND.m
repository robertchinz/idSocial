function output_parfor = idSocial_function_wrapper_ExecParforRAND(paramspath)
if ischar(paramspath)
    load(paramspath)
elseif isstruct(paramspath)
    local_exec = true;
end

slave_id=paramstruct.slave_id;
idx_combs_temp=paramstruct.idx_combs_temp;
no_frames_temp=paramstruct.no_frames_temp; 
info_no_frames_part_array_temp=paramstruct.info_no_frames_part_array_temp; 
rand_argin_scaled_temp=paramstruct.rand_argin_scaled_temp; 
temp_save=paramstruct.temp_save; 
function_handle=paramstruct.function_handle; 
temp_output_folder=paramstruct.slave_output_folder{slave_id}; 
temp_savepath_root=paramstruct.master_output_folder; 
no_temp_saves=paramstruct.no_temp_saves;
temp_idx = paramstruct.temp_idx;
info_framerate_temp = paramstruct.info_framerate_temp;
info_bodylength_temp = paramstruct.info_bodylength_temp;
trajargin = paramstruct.trajargin;
applyfilter = paramstruct.applyfilter;
order_nb_function_temp= paramstruct.order_nb_function_temp;
no_funcs= paramstruct.no_funcs;
no_outfuncs=paramstruct.no_outfuncs;
good_idx= paramstruct.good_idx;
no_combs_per_temp_save=paramstruct.no_combs_per_temp_save;

output_parfor=cell(no_combs_per_temp_save,sum(good_idx));


psize=4;
myCluster = parcluster('local');
if psize>myCluster.NumWorkers
    psize=myCluster.NumWorkers;
end

try
    no_cores_open=matlabpool('size');
    if no_cores_open~=psize && (psize~=Inf || no_cores_open~=feature('numCores'))
        if no_cores_open~=0
            matlabpool close
        end
        if psize==Inf
            matlabpool open % Use default configuration
        else
            matlabpool('open','local',psize)
        end
    end
catch
    disp([mfilename ': Could not open matlabpool. Continue without parallel processing.'])
end

if exist(temp_output_folder,'dir')~=7
    mkdir(temp_output_folder)
end


parfor idx_count=  1:numel(temp_idx)
    
    idces_act=idx_combs_temp(idx_count,:);
    
    group=idces_act(1);
    subset=idces_act(2);
    trial=idces_act(3);
    part=idces_act(4);
    no_frames=no_frames_temp(idx_count,:);
    
    
    idx1=max((part-2)*floor(info_no_frames_part_array_temp(idx_count,:)/2)+1,1);
    idx2=min(floor(info_no_frames_part_array_temp(idx_count,:)/2)*part+1,no_frames);
    
    
    
    if idx2>idx1
        fprintf('Random, Group %i, subset %i, trial %i, time %.1f-%.1f min\n',group,subset,trial,idx1/info_framerate_temp(idx_count)/60,idx2/info_framerate_temp(idx_count)/60);
        %             fprintf('.')
       
        rand_argin_act=rand_argin_scaled_temp{idx_count};
        %%
        rand_tr=load(rand_argin_scaled_temp{idx_count}{trajargin});
        try
            rand_tr=rand_tr.rand_tr(limit_frames_act(1):min(size(rand_tr.rand_tr,1),limit_frames_act(2)),:,:,:);
        catch
            keyboard
        end
        rand_tr=rand_tr(idx1:idx2,:,:,:);
        if applyfilter
            rand_tr = idSocial_filters3D(rand_tr,...
                options,...
                info_bodylength_temp(idx_count),...
                info_framerate_temp(idx_count));
        end
        
        if order_nb_function_temp(idx_count)
            % I cannot load the already
            % existing ordered trajectory if I want
            % to apply filters before!
            rand_tr=idSocial_auxiliaries_nearestneighbours(rand_tr);
        end
        
        rand_argin_act{trajargin}=rand_tr;
        
        rand_alloutput=cell(no_funcs,1);
        [rand_alloutput{:}]=feval(function_handle,rand_argin_act{:});
        rand_output_parfor(idx_count,:)=rand_alloutput(good_idx);
        %%
    end
    
end

if no_temp_saves>1 
    ofbufferTic = tic;
    for of = 1:no_outfuncs
        disp(['Buffering ' 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat' ' ...'])
        tempOutput =  rand_output_parfor(:,of);
        savefast([temp_output_folder 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'tempOutput');
        dummy = 'dummy';
        try
        save([temp_output_folder 'REMOTE_' num2str(slave_id) '_READY.mat'],'dummy');
        catch
            keyboard
        end
    end
    no_cores_open=matlabpool('size');
    if no_cores_open~=0;
        matlabpool close;
    end
    disp(['Done (' num2str(toc(ofbufferTic),'%.1f') 's)'])
end
if ischar(paramspath)
    delete(paramspath)
end
