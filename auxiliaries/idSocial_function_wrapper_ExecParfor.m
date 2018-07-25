function output_parfor = idSocial_function_wrapper_ExecParfor(paramspath)
if ischar(paramspath)
    load(paramspath)
elseif isstruct(paramspath)
    local_exec = true;
end

slave_id=paramstruct.slave_id;
idx_combs_temp=paramstruct.idx_combs_temp;
no_frames_temp=paramstruct.no_frames_temp; 
info_no_frames_part_array_temp=paramstruct.info_no_frames_part_array_temp; 
argin_scaled_temp=paramstruct.argin_scaled_temp; 
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
        fprintf('Group %i, subset %i, trial %i, time %.1f-%.1f min\n',group,subset,trial,idx1/info_framerate_temp(idx_count)/60,idx2/info_framerate_temp(idx_count)/60);
        %             fprintf('.')
        argin_act=argin_scaled_temp{idx_count};
        if ischar(argin_scaled_temp{idx_count}{trajargin}) % && exist(argin_scaled{idx_count}{trajargin},'file')==2
            
            % The path may not be the same on various
            % computers:
            argin_scaled_temp{idx_count}{trajargin}(1) = temp_output_folder(1);
            
            
            tr=load(argin_scaled_temp{idx_count}{trajargin});
            opt=tr.tr.options;
            tr=tr.tr.trajectory(opt.start_frame:min(size(tr.tr.trajectory,1),opt.end_frame),:,:,:);
            tr=tr(idx1:idx2,:,:,:);
            if applyfilter
                tr = idSocial_filters3D(tr,...
                    options,...
                    info_bodylength_temp(idx_count),...
                    info_framerate_temp(idx_count));
            end
            
            if order_nb_function_temp(idx_count)
                % I cannot load the already
                % existing ordered trajectory if I want
                % to apply filters before!
                % Skip if ordering has been applied before.
                tr=idSocial_auxiliaries_nearestneighbours(tr);
            end
            argin_act{trajargin}=tr;
        end
        %             if subset==21 && trial==1 && part==2; keyboard; end
        alloutput=cell(no_funcs,1);
        [alloutput{:}]=feval(function_handle,argin_act{:});
        output_parfor(idx_count,:)=alloutput(good_idx);
        
    end
    
end

if no_temp_saves>1 
    ofbufferTic = tic;
    for of = 1:no_outfuncs
        disp(['Buffering ' 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat' ' ...'])
        tempOutput =  output_parfor(:,of);
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
