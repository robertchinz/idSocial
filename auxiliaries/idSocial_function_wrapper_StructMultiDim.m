function input_data=idSocial_function_wrapper_StructMultiDim(input_data, options, def_options_function, plot_mode, functionInfo, external_func_name)

mb_per_temp_save=2000;
skip_end_parts = true;
% Number of workers for parallel computations:
psize=4;
manual_actualize=false;
manual_actualize_plot_mode = true;
use_parallel=true;
use_remote = false;
min_number_loop=50; % Min. number of loop iterations for which parallel computation will be used
min_data_lim4parallel = 0; % Min. number of trajectories for which parallel computation will be used

myCluster = parcluster('local');
if psize>myCluster.NumWorkers
    psize=myCluster.NumWorkers;
end

if nargin < 6 || isempty(external_func_name)
    external_func_name = '';
end

%% Default options
minx = -inf;
maxx = inf;
miny = -inf;
maxy = inf;
maxr = inf;
center = [0,0,0];

% def_options=def_options_function;
def_options.statistical_test_type='ranksum';%'bootstrap'; %
def_options.statistical_test_significance=0.05;
def_options.plot_titlefontsize=16;
def_options.plot_axislabelsize=16;
def_options.plot_labelfontsize=14;
def_options.plot_legendfontsize=14;
def_options.plot_gr_tr_fontsize=12;
def_options.plot_tickfontsize=14;
def_options.project_path='';
% n = size(get(gcf,'Colormap'),1);
n = 128;
c = get(0,'defaultAxesColorOrder');
cmap = c(rem(0:n-1,size(c,1))+1,:);
def_options.plot_colororder=cmap;%[1 0 0; 0 1 0; 0 0 1; .7 .5 .8; 0 1 1; 1 1 0; .5 0 1; 1 .5 0; 0 1 .5; 0 .4 1; 1 0 1; 1 .5 .5; 0 1 1; 1 1 0;  .25 .5 .5; 1 .5 0; 0 .6 0; .5 .2 1; .5 1 1; 1 1 0];
def_options.plot_linestyle={'-'};% '.-'; '--'; '-.'; '+-'; '*-'; 'o-'; 'x-'; 's-'; '^-'; 'd-'};
def_options.plot_linewidth=2.5;
def_options.significance_between_groups=true;
def_options.significance_level=0.05;
def_options.plot_mode.newfigure={[]};
def_options.plot_mode.subplots={[]};
def_options.plot_mode.legend={''};
def_options.plot_mode.xaxis={'Time'};
def_options.plot_mode.data={[]};
def_options.plot_mode.display_mode='plot2d';
def_options.group_names={};
def_options.plot_xticklabel=[];
def_options.plot_axis_scaling=[1 1];
def_options.edges=-4:.5:4;
def_options.spacing=1;
def_options.random_data=false;
def_options.filter_focal_speedlimits_bl_per_s=[];
def_options.filter_neighbor_speedlimits_bl_per_s=[];
def_options.filter_neighbor_accelerationlimits_bl_per_s2 = [];
def_options.filter_focal_accelerationlimits_bl_per_s2  = [];
def_options.filter_distancelimits_bl=[];
def_options.filter_focal_rectangularROI = [];
def_options.filter_neighbor_rectangularROI = [];
def_options.filter_focal_circularROI = [];
def_options.filter_neighbor_circularROI = []; def_options.filter_AllMembersPresent = -inf; def_options.filter_AllMembersPresentWorstFocal = -inf; def_options.filter_WorstIndividual = -inf;
def_options.filter_focal_spatial_sectors = []; %'lateral' or 'frontal'
% Options for maps:
def_options.rowlabelfontsize=12;
% Options for territories:
def_options.discard_bottom_percent=10;
def_options.plot_bootstrap_repetitions=10000;
def_options.plot_bootstrap_statfunc='Mean';
def_options.plot_deviation=true;
def_options.temp_savepath='D:\~idSocialTemp\';
def_options.temp_savepath_slave='';
def_options.order_neighbors = false;
def_options.internal_random_controls = [];

% Get save path
if ~isfield(options,'temp_savepath_slave')
    options.temp_savepath_slave='';
end

[~, temp_savepath]= idSocial_recursiveGetOptionsFromOptionsCell(input_data(1,1,1).options,'temp_savepath');
if isfield(options,'temp_savepath') && ~isempty(options.temp_savepath)
    temp_savepath = options.temp_savepath;
end
if ~strcmpi(temp_savepath(end-3:end),'.mat') % temp_savepath is directory
    if ~strcmp(temp_savepath(end),'\')
        temp_savepath = [temp_savepath '\'];
    end
    temp_savepath_root = temp_savepath;
    temp_savepath = [temp_savepath 'input_data.mat'];
else % temp_savepath is file
    delidx=strfind(temp_savepath,'\');
    temp_savepath_root = temp_savepath(1:max(delidx));
end
if exist(temp_savepath_root,'dir')~=7
    mkdir(temp_savepath_root)
end

% Slave/remote
if isfield(options,'temp_savepath_slave') && ~isempty(options.temp_savepath_slave)
    temp_savepath_slave_full = options.temp_savepath_slave;
else
    temp_savepath_slave_full = [];
end
no_remotes = 1;
if ~isempty(temp_savepath_slave_full)
    
    if ~iscell(temp_savepath_slave_full)
        temp_savepath_slave_full = {temp_savepath_slave_full};
    end
    
    temp_savepath_slave_full = vertcat(temp_savepath_root,temp_savepath_slave_full);
    
    remote_computers = cell(numel(temp_savepath_slave_full),1);
    temp_savepath_slave = cell(numel(temp_savepath_slave_full),1);
    temp_savepath_slave_root = cell(numel(temp_savepath_slave_full),1);
    no_remotes = numel(temp_savepath_slave_full);
    
    temp_savepath_slave{1} = temp_savepath_root;
    for rm = 2:numel(temp_savepath_slave_full)
        remote_computers{rm} = temp_savepath_slave_full{rm}(1: strfind(temp_savepath_slave_full{rm},':')-3);
        temp_savepath_slave{rm} = temp_savepath_slave_full{rm}(strfind(temp_savepath_slave_full{rm},':')-1:end);
    end
    
    savefast([temp_savepath_root 'REMOTE_1_READY.mat'],'dummy');
    for rm = 1:numel(temp_savepath_slave_full)
        if ~strcmpi(temp_savepath_slave{rm}(end-3:end),'.mat') % temp_savepath is directory
            if ~strcmp(temp_savepath_slave{rm}(end),'\')
                temp_savepath_slave{rm} = [temp_savepath_slave{rm} '\'];
            end
            temp_savepath_slave_root{rm} = temp_savepath_slave{rm};
            temp_savepath_slave{rm} = [temp_savepath_slave{rm} 'input_data.mat'];
        else % temp_savepath is file
            delidx=strfind(temp_savepath_slave{rm},'\');
            temp_savepath_slave_root{rm} = temp_savepath_slave{rm}(1:max(delidx));
        end
        
        if exist(temp_savepath_root,'dir')~=7
            mkdir(temp_savepath_root)
        end
        savefast([temp_savepath_root 'REMOTE_' num2str(rm) '_READY.mat'],'dummy');
    end
else
    use_remote = false;
end

diary([temp_savepath_root 'log_' datestr(now, 30) '.txt'])

def_options_fct_fldnames=fieldnames(def_options_function);
for k=1:size(def_options_fct_fldnames,1)
    def_options.(def_options_fct_fldnames{k})=def_options_function.(def_options_fct_fldnames{k});
end

act_method=def_options.act_method;
core_act_method=strrep(func2str(functionInfo.handle),'idSocial_','');

already_exists =  isfield(input_data(1,1),act_method);

[input_data, options]=idSocial_readparams(input_data,options,def_options,act_method);

order_nb_function_single = options.order_neighbors;
function_option_fields=intersect(fieldnames(input_data(1,1).(act_method).options),fieldnames(def_options_function));



plotfields=strfind(function_option_fields,'plot');
plotfields=cellfun(@(x)~isempty(x),plotfields);
savepathfield = strcmpi(function_option_fields,'temp_savepath');
function_option_fields(ismember(function_option_fields,'act_method') | plotfields | savepathfield  )=[];
options_vals=cellfun(@(x) input_data(1,1).(act_method).options.(x),function_option_fields,'UniformOutput',false);
% keyboard
if isfield(input_data(1,1),core_act_method) && isfield(input_data(1,1).(core_act_method),'options') && ...
        already_exists
    %
    %
    %     % Check if base data has changed:
    %     new_data = false;
    %     if isfield(input_data(1,1).(core_act_method),'trajectory_ver') && ...
    %     end
    %     %%%%%%%
    %
    %
    
    func_opts_vals=cellfun(@(x) input_data(1,1).(core_act_method).options.(x),function_option_fields,'UniformOutput',false);
    actualize=~isequal(options_vals,func_opts_vals);
    if actualize
        disp('Previous data found. The following options have changed:')
        
        changed_opts = cellfun(@(x,y) ~isequal(x,y),options_vals,func_opts_vals);
        %         no_changed = sum(changed_opts);
        for k=find(changed_opts)'
            disp([upper(function_option_fields{k}) ': '])
            disp('Old value: '); display(func_opts_vals{k});
            disp('New value: '); display(options_vals{k});
        end
    end
else
    if ~isfield(input_data(1,1),core_act_method)
        disp('No previous data found. Start calculations...')
    end
    actualize=1;
end
%     if actualize==1; keyboard; end


if isfield(input_data(1,1),act_method) && isfield(input_data(1,1).(act_method),'plot_mode') && ~(actualize || manual_actualize)
    old_plot_mode = fieldnames(input_data(1,1).(act_method).plot_mode);
    new_plot_mode = fieldnames(plot_mode);
    old_plot_mode_vals=cellfun(@(x) input_data(1,1).(act_method).plot_mode.(x),old_plot_mode,'UniformOutput',false);
    new_plot_mode_vals=cellfun(@(x) plot_mode.(x),new_plot_mode,'UniformOutput',false);
    
    actualize_plot_mode=~isequal(old_plot_mode,new_plot_mode) || ~isequal(old_plot_mode_vals,new_plot_mode_vals);
else
    actualize_plot_mode=1;
end


% keyboard
%% Get info
info=               input_data(1,1,1).info;
blpxl=              info.blpxl;
framerate=          info.framerate;
duration=           info.duration;
no_frames=          info.no_frames;
trials=             info.trials;
no_focals=          info.no_focals;
no_neighbors=          info.no_neighbors;
no_neighborsRAND=          info.no_neighborsRAND;
no_focalsRAND=          info.no_neighborsRAND; % The same as neighbors!
% disp([mfilename ': Line 74, CHANGE!!!'])
no_groups=        info.no_groups;
no_subsets=        info.no_subsets;
no_trials=         max(info.no_trials(:));


%% Some standard parameters
act_method=def_options.act_method;
randomized_calcs=options.random_data;
max_no_focals=max(max(max(no_focals)));
max_no_neighbors=max(max(max(no_neighbors)));
max_no_neighborsRAND=max(max(max(no_neighborsRAND)));

max_no_frames=max(max(max(no_frames)));

timeintervals_in_min=options.timeintervals_in_min;
if isempty(timeintervals_in_min) % || timeintervals_in_min> min(min(floor(min(no_frames./framerate)/60)));
    timeintervals_in_min=duration;
end;
no_frames_part_array=floor(timeintervals_in_min.*framerate*60)-1;
no_parts=2*floor(no_frames./no_frames_part_array)+1;
max_no_parts=max(max(max(no_parts)));
if max_no_parts ==3 && skip_end_parts
    skip_end_parts = true;
end
min_no_parts=min(min(min(no_parts)));
discard_bottom_percent=options.discard_bottom_percent;

if ~isfield(functionInfo,'input_params_scaling') || isempty(functionInfo.input_params_scaling)
    functionInfo.input_params_scaling=cell(size(functionInfo.input_params,1),1);
    for k=1:size(functionInfo.input_params,1)
        functionInfo.input_params_scaling{k}=ones(no_groups,max(no_subsets),no_trials);
    end
end

if isfield(functionInfo,'input_params_scaling') && ~isempty(functionInfo.input_params_scaling)
    for k=1:size(functionInfo.input_params,1)
        if all(size(functionInfo.input_params_scaling{k})==1)
            functionInfo.input_params_scaling{k}=repmat(functionInfo.input_params_scaling{k},[no_groups,max(no_subsets),no_trials]);
        end
    end
end


%% Execute function
char_input_params=cell(size(functionInfo.input_params,1),1);
trajargin_idx=false(size(functionInfo.input_params,1),1);
trajargin=false(size(functionInfo.input_params,1),1);

%                         fldmdargin=false(size(functionInfo.input_params,1),1);
trajargin_idx(cellfun(@(x) ischar(x), functionInfo.input_params))=true;
char_input_params(trajargin_idx)=functionInfo.input_params(trajargin_idx);
char_input_params(~trajargin_idx)={''};
trajargin(...
    ismember(char_input_params, {'trajectory'}))=true;

if (all(cellfun(@(x) ischar(x),{input_data(:).movementdata})) && ~all(cellfun(@(x) exist(x,'file')==2,{input_data(:).movementdata})) && any(trajargin)) % The last condition checks if trajectories are temporarily saved and needed
    error('Could not find files in temporary folder.')
end


filterfields=strfind(fieldnames(options),'filter');
filterfields=cellfun(@(x)~isempty(x),filterfields);
optionfieldnames = fieldnames(options);
applyfilter = any(cellfun(@(x) ~isempty(options.(x)),optionfieldnames(filterfields)));

name_outfuncs=functionInfo.output2function(~strcmpi(functionInfo.output2function,'temp') & cellfun(@(x) ~isempty(x),functionInfo.output2function));

function_handle = functionInfo.handle;
if  actualize || manual_actualize
    disp([act_method ': Updating results...'])
    %     if randomized_calcs
    no_outfuncs=sum(~strcmpi(functionInfo.output2function,'temp') & cellfun(@(x) ~isempty(x),functionInfo.output2function));
    %     else
    %         no_outfuncs=sum(cellfun(@(x) isempty(x),strfind(functionInfo.output2function,'RAND')) & ...
    %             ~strcmpi(functionInfo.output2function,'temp') & cellfun(@(x) ~isempty(x),functionInfo.output2function));
    %     end
    
    
    
    
    
    fldoptargin=false(size(functionInfo.input_params,1),1);
    fldinfoargin=false(size(functionInfo.input_params,1),1);
    fldoptargin(...
        ismember(char_input_params, strcat('options.',fieldnames(options))))=true;
    fldinfoargin(...
        ismember(char_input_params, strcat('info.',fieldnames(input_data(1,1,1).info))))=true;
    %                         fldmdargin(...
    %                             ismember(char_input_params, 'movementdata') & ~ismember(char_input_params, 'movementdata.'))=true;
    
    argin=functionInfo.input_params;
    %                         mdargin=cellfun(@(x) act_md.(x),strrep(char_input_params(trajargin),'movementdata.',''),'UniformOutput',false);
    
    if any(fldoptargin)
        argin(fldoptargin)=cellfun(@(x) options.(x),strrep(char_input_params(fldoptargin),'options.',''),'UniformOutput', false);
    end
    
    %     if ~randomized_calcs
    %         good_idx=cellfun(@(x) ~isempty(x)&& ~strcmpi(x,'temp')&& isempty(strfind(x,'RAND')) ,functionInfo.output2function);
    %     else
    good_idx=cellfun(@(x) ~isempty(x)&& ~strcmpi(x,'temp'),functionInfo.output2function);
    %     end
    input_data(1,1,1).(core_act_method).options=cell2struct(options_vals,function_option_fields);
    
    comb_count=1;
    idx_combs=[];
    
    for group=1:no_groups
        for subset=1:no_subsets(group)
            if size(trials{group,subset},2)==1; trlist = trials{group,subset}'; else trlist = trials{group,subset}; end
            for trial= trlist
                no_frames=info.no_frames(group,subset,trial);
                for part=1:no_parts(group,subset,trial)
                    if ~skip_end_parts || part==2
                        idx1=max((part-2)*floor(no_frames_part_array(group,subset,trial)/2)+1,1);
                        idx2=min(floor(no_frames_part_array(group,subset,trial)/2)*part+1,no_frames);
                        if idx2>idx1
                            idx_combs=vertcat(idx_combs,[group subset trial part]);
                            comb_count=comb_count+1;
                        end
                    end
                end
            end
        end
    end
    
    
    no_combs=size(idx_combs,1);
    info_no_frames=NaN(no_combs,1);
    info_framerate=NaN(no_combs,1);
    info_bodylength=NaN(no_combs,1);
    info_circroi=NaN(no_combs,4);
    info_no_frames_part_array=NaN(no_combs,1);
    info_no_focals=NaN(no_combs,1);
    order_nb_trial=NaN(no_combs,1);
    order_nb_function=NaN(no_combs,1);
    %     function_handle=cell(no_combs,1);
    argin_scaled=cell(no_combs,1);
    rand_argin_scaled=cell(no_combs,1);
    comb_count=1;
    for group=1:no_groups
        for subset=1:no_subsets(group)
            if size(trials{group,subset},2)==1; trlist = trials{group,subset}'; else trlist = trials{group,subset}; end
            for trial= trlist
                
                no_frames=info.no_frames(group,subset,trial);
                for part=1:no_parts(group,subset,trial)
                    if ~skip_end_parts || part==2
                        order_nb_trial(comb_count)=input_data(1,1,1).options{group}{subset}{trial}.order_neighbors;
                        order_nb_function(comb_count) = order_nb_function_single && ~input_data(1,1,1).options{group}{subset}{trial}.order_neighbors;
                        
                        idx1=max((part-2)*floor(no_frames_part_array(group,subset,trial)/2)+1,1);
                        idx2=min(floor(no_frames_part_array(group,subset,trial)/2)*part+1,no_frames);
                        
                        if idx2>idx1
                            info_no_focals(comb_count)=info.no_focals(group,subset,trial);
                            info_no_frames(comb_count)=info.no_frames(group,subset,trial);
                            info_framerate(comb_count)=info.framerate(group,subset,trial);
                            info_bodylength(comb_count)=info.blpxl(group,subset,trial);
                            if isfield(info,'circ_roi')
                                info_circroi(comb_count,:)=squeeze(info.circ_roi(group,subset,trial,:));
                            end
                            info_no_frames_part_array(comb_count)=no_frames_part_array(group,subset,trial);
                            
                            argin_parfor=argin;
                            if any(trajargin)
                                argin_parfor{trajargin}=input_data(group,subset,trial).movementdata;
                                
                            end
                            if any(fldinfoargin)
                                fldstr=strrep(char_input_params,'info.','');
                                for k=find(fldinfoargin)'
                                    
                                    argin_parfor{k}=squeeze(input_data(1,1,1).info.(fldstr{k})(group,subset,trial,:));
                                end
                                %                             argin_parfor(fldinfoargin)=cellfun(@(x) input_data(1,1,1).info.(x)(group,subset,trial),strrep(char_input_params(fldinfoargin),'info.',''),'UniformOutput', false);
                            end
                            
                            argin_scaled{comb_count}=argin_parfor;
                            
                            argin_scaled{comb_count}(cellfun(@(x) isnumeric(x),argin_parfor))=...
                                cellfun(@(x,y) x.*squeeze(y(group,subset,trial,:,:,:))',argin_parfor(cellfun(@(x) isnumeric(x),argin_parfor)),...
                                functionInfo.input_params_scaling(cellfun(@(x) isnumeric(x),argin_parfor)),'UniformOutput',false);
                            
                            if randomized_calcs
                                rand_argin_parfor=argin_parfor;
                                rand_argin_parfor{trajargin}=input_data(group,subset,trial).rand_movementdata;
                                rand_argin_scaled{comb_count}=rand_argin_parfor;
                                rand_argin_scaled{comb_count}(cellfun(@(x) isnumeric(x),rand_argin_parfor))=...
                                    cellfun(@(x,y) x.*squeeze(y(group,subset,trial,:,:,:))',rand_argin_parfor(cellfun(@(x) isnumeric(x),rand_argin_parfor)),...
                                    functionInfo.input_params_scaling(cellfun(@(x) isnumeric(x),rand_argin_parfor)),'UniformOutput',false);
                                
                            end
                            
                                                     
                            comb_count=comb_count+1;
                        end
                    end
                end
            end
        end
    end
    
    %% Execute function
    
    
    
    no_funcs=size(functionInfo.output2function,1);
    %     function_handle=functionInfo.handle;
    limit_frames = NaN(no_combs,2);
    if use_parallel && no_combs>min_number_loop && numel(input_data) > min_data_lim4parallel
        try
            if ~verLessThan('matlab', '8.3.0')
                p = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(p)
                    no_cores_open = 0;
                else
                    no_cores_open = p.NumWorkers;
                end
                
                if no_cores_open~=psize && (psize~=Inf || no_cores_open~=feature('numCores'))
                    if no_cores_open~=0
                        delete(gcp('nocreate'));
                    end
                    if psize==Inf
                        parpool; % Use default configuration
                    else
                        parpool('local',psize)
                    end
                end
            else
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
            end
            
            
        catch poolerror
            disp([mfilename ': Could not open parallel pool. Continue without parallel processing.'])
            disp('Error message is:')
            disp(poolerror)
        end
        
    end
    
    % Execute function once to estimate output size:
    idx_count=min(2,no_combs);
    
    argin_act=argin_scaled{idx_count};
    if ischar(argin_scaled{idx_count}{trajargin}) % && exist(argin_scaled{idx_count}{trajargin},'file')==2
        
        tr=load(argin_scaled{idx_count}{trajargin});
        opt=tr.tr.options;
        tr=tr.tr.trajectory(opt.start_frame:min(size(tr.tr.trajectory,1),opt.end_frame),:,:,:);
        if applyfilter
            tr = idSocial_filters3D(tr,...
                options,...
                info_bodylength(idx_count),...
                info_framerate(idx_count),...
                info_circroi(idx_count,:));
        end
        
        if order_nb_function(idx_count)
            % I cannot load the already
            % existing ordered trajectory if I want
            % to apply filters before!
            % Skip if ordering has been applied before.
            tr=idSocial_auxiliaries_nearestneighbours(tr);
        end
        argin_act{trajargin}=tr;
    end
    alloutput=cell(no_funcs,1);
    try
        [alloutput{:}]=feval(function_handle,argin_act{:});
    catch
        keyboard
    end
    output_estimate=whos('alloutput');
    output_estimate = output_estimate.bytes*9.53674e-7;
    no_temp_saves = ceil(no_combs*output_estimate/mb_per_temp_save);
    no_combs_per_temp_save = ceil(no_combs/no_temp_saves);
    
    idx_conversion=NaN(no_combs,3);
    
    
    % The following loop is only there to create
    % 'limit_frames'. There sure is a shorter, nicer way!
    for temp_save = 1 : no_temp_saves
        temp_idx = (temp_save-1)*no_combs_per_temp_save+1 : min(temp_save*no_combs_per_temp_save,no_combs);
        idx_combs_temp = idx_combs(temp_idx,:);
        no_frames_temp = info_no_frames(temp_idx,:);
        info_no_frames_part_array_temp = info_no_frames_part_array(temp_idx,:);
        limit_frames_temp = NaN(min(no_combs_per_temp_save,numel(temp_idx)),2);
        for idx_count=  1:numel(temp_idx)
            idces_act=idx_combs_temp(idx_count,:);
            
            group=idces_act(1);
            subset=idces_act(2);
            trial=idces_act(3);
            part=idces_act(4);
            no_frames=no_frames_temp(idx_count,:);
            
            
            idx1=max((part-2)*floor(info_no_frames_part_array_temp(idx_count,:)/2)+1,1);
            idx2=min(floor(info_no_frames_part_array_temp(idx_count,:)/2)*part+1,no_frames);
            if idx2>idx1
                limit_frames_temp(idx_count,:)=[input_data(1,1,1).options{group}{subset}{trial}.start_frame, ...
                    input_data(1,1,1).options{group}{subset}{trial}.end_frame];
                
            end
            limit_frames(temp_idx,:) = limit_frames_temp;
        end
        
    end
    for temp_save = 1 : no_temp_saves
        output_parfor=cell(no_combs_per_temp_save,sum(good_idx));
        temp_idx = (temp_save-1)*no_combs_per_temp_save+1 : min(temp_save*no_combs_per_temp_save,no_combs);
        idx_combs_temp = idx_combs(temp_idx,:);
        no_frames_temp = info_no_frames(temp_idx,:);
        info_no_frames_part_array_temp = info_no_frames_part_array(temp_idx,:);
        info_framerate_temp = info_framerate(temp_idx);
        argin_scaled_temp = argin_scaled(temp_idx);
        %         limit_frames_temp = NaN(min(no_combs_per_temp_save,numel(temp_idx)),2);
        info_bodylength_temp = info_bodylength(temp_idx);
        order_nb_function_temp = order_nb_function(temp_idx);
        idx_conversion(temp_idx,1)=temp_idx;
        
        idx_conversion(temp_idx,2)=1:min(no_combs_per_temp_save,numel(temp_idx));
        
        idx_conversion(temp_idx,3)=temp_save;
        
        if ~use_remote || no_temp_saves == 1
            %par...
            for idx_count=  1:numel(temp_idx)
                
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
                        
                        tr=load(argin_scaled_temp{idx_count}{trajargin});
                        opt=tr.tr.options;
                        %                         limit_frames_temp(idx_count,:)=[opt.start_frame opt.end_frame];
                        tr=tr.tr.trajectory(opt.start_frame:min(size(tr.tr.trajectory,1),opt.end_frame),:,:,:);
                        tr=tr(idx1:idx2,:,:,:);
                        if applyfilter
                            tr = idSocial_filters3D(tr,...
                                options,...
                                info_bodylength_temp(idx_count),...
                                info_framerate_temp(idx_count),...
                                info_circroi(idx_count,:));
                            worst_indiv = min(min(sum(~isnan(tr(:,:,:,1)),1)./size(tr,1)));
                            fprintf('\t%.2f%% of frames left for worst focal-neighbor pair.\n',worst_indiv*100)
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
%                                                     if subset==10 && trial==1 && part==2; keyboard; end
                    alloutput=cell(no_funcs,1);
                    %                     if worst_indiv>0
                    try
                        
                        [alloutput{:}]=feval(function_handle,argin_act{:});
                        
                        
                        
                    catch exception
                        %                             keyboard
                        
                        msgString = getReport(exception);
                        warning([mfilename ': Execution failed. Error message reads:'])
                        disp(msgString)
                    end
                    if order_nb_function_temp(idx_count)
                        alloutput_new = alloutput;
                        for outc = 1:size(alloutput,1)
                            for ff = 1:size(alloutput{outc},1)
                                for nf = 1:size(alloutput{outc},2)
                                    if nf<ff
                                        alloutput_new{outc}(ff,nf,:,:,:,:,:) = alloutput{outc}(ff,nf,:,:,:,:,:);
                                    elseif nf>=ff && nf<size(alloutput{outc},2)
                                        alloutput_new{outc}(ff,nf,:,:,:,:,:) = alloutput{outc}(ff,nf+1,:,:,:,:,:);
                                    elseif  nf==size(alloutput{outc},2)
                                        alloutput_new{outc}(ff,nf,:,:,:,:,:) = alloutput{outc}(ff,ff,:,:,:,:,:);
                                    end
                                    
                                end
                            end
                        end
                        alloutput = alloutput_new;
                        %                         clear alloutput_new
                    end
                    output_parfor(idx_count,:)=alloutput(good_idx);
                    %                     end
                end
                
            end
            if no_temp_saves>1
                
                ofbufferTic=tic;
                for of = 1:no_outfuncs
                    disp(['Buffering ' 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat' ' ...'])
                    tempOutput =  output_parfor(:,of);
                    savefast([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'tempOutput');
                end
                %                                savefast([temp_savepath_root 'wrapper_buffer_' num2str(temp_save) '.mat'],'output_parfor');
                
                disp(['Done (' num2str(toc(ofbufferTic),'%.1f') 's)'])
                clear output_parfor tempOutput
            end
            %         limit_frames(temp_idx,:) = limit_frames_temp;
        else % Remote
            paramstruct.idx_combs_temp = idx_combs_temp;
            paramstruct.no_frames_temp = no_frames_temp;
            paramstruct.info_bodylength_temp = info_bodylength_temp;
            paramstruct.info_no_frames_part_array_temp = info_no_frames_part_array_temp;
            paramstruct.argin_scaled_temp = argin_scaled_temp;
            paramstruct.temp_save = temp_save;
            paramstruct.function_handle = function_handle;
            paramstruct.slave_output_folder = temp_savepath_slave_root;
            paramstruct.master_output_folder = temp_savepath_root;
            paramstruct.no_temp_saves = no_temp_saves;
            paramstruct.temp_idx = temp_idx;
            paramstruct.info_framerate_temp = info_framerate(temp_idx);
            paramstruct.trajargin=trajargin;
            paramstruct.applyfilter=applyfilter;
            paramstruct.order_nb_function_temp = order_nb_function_temp;
            paramstruct.no_funcs = no_funcs;
            paramstruct.no_outfuncs = no_outfuncs;
            paramstruct.good_idx = good_idx;
            paramstruct.no_combs_per_temp_save = no_combs_per_temp_save;
            
            act_remote = mod(temp_save-1,no_remotes)+1;
            
            temp_save_delivered = false;
            while ~temp_save_delivered
                if exist([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'],'file')==2
                    disp(['Starting round ' num2str(temp_save) ' on Computer ' num2str(act_remote)])
                    delete([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'])
                    paramstruct.slave_id = act_remote;
                    parampath_save = [temp_savepath_slave_root{1} 'inputparams_no_' num2str(temp_save) '.mat'];
                    savefast(parampath_save,'paramstruct');
                    parampath_load = [temp_savepath_slave_root{act_remote} 'inputparams_no_' num2str(temp_save) '.mat'];
                    if act_remote == 1 % Local
                        %                         system([' "C:\Program Files\Matlab\R2012a\bin\matlab.exe" "-nodisplay" "-nosplash" "-nodesktop" "-r" "idSocial_function_wrapper_ExecParfor ' parampath '"']);
                        system([' "C:\Program Files\Matlab\R2012a\bin\matlab.exe"  -nodisplay -nosplash -nodesktop -minimize -r "idSocial_function_wrapper_ExecParfor ' parampath_load ';exit;"']);
                        
                    else
                        %                         keyboard
                        system(['"C:\Program Files (x86)\PSTools\psexec" ' remote_computers{act_remote} ' "C:\Program Files\Matlab\R2012a\bin\matlab.exe"  "-nodisplay" "-nosplash" "-nodesktop" "-minimize" "-r" "idSocial_function_wrapper_ExecParfor ' parampath_load ' exit;" ']);
                    end
                    
                    temp_save_delivered = true;
                else
                    act_remote = mod(act_remote + 1,no_remotes)+1;
                end
            end
            
        end
    end
    
    all_done=false;
    done_array = false(no_temp_saves,no_outfuncs);
    if use_remote && no_temp_saves>1
        while ~all_done
            for temp_save = 1 : no_temp_saves
                for of = 1:no_outfuncs
                    if exist([temp_savepath_slave_root{1} 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'file')==2
                        done_array(temp_save,of)=true;
                    end
                end
            end
            pause(.01)
            all_done = all(done_array(:));
        end
    end
    
    for act_remote = 1: no_remotes
        if exist([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'],'file')==2
            delete([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'])
        end
    end
    
    
    for outfun=1:no_outfuncs
        if strcmpi(name_outfuncs{outfun}(end-3:end),'RAND') && ~isempty(options.internal_random_controls)
            if isnumeric(options.internal_random_controls)
                no_neighbors_act = NaN(size(no_neighbors));
                no_neighbors_act(~isnan(no_neighbors)) = options.internal_random_controls;
            elseif islogical(options.internal_random_controls) && options.internal_random_controls
                no_neighbors_act = no_neighbors;
            end
            
        else
            no_neighbors_act = no_neighbors;
        end
        output_new=cell(no_groups,max(no_subsets),no_trials,max_no_parts,max_no_focals,max_no_neighbors);
        temp_save = 0;
        for idx_count=1:no_combs
            
            temp_save_act = idx_conversion(idx_count,3);
            idx_count_temp = idx_conversion(idx_count,2);
            
            if temp_save_act ~= temp_save && no_temp_saves>1
                disp(['Loading temporary buffer ' 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat' ' ...'])
                %                 keyboard
                load([temp_savepath_root 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat'])
                output_parfor_act = tempOutput;
                disp('Done.')
                temp_save = temp_save_act;
            end
            if no_temp_saves==1
                output_parfor_act = output_parfor(:,outfun);
            end
            
            
            idces_act=idx_combs(idx_count,:);
            group=idces_act(1);
            subset=idces_act(2);
            trial=idces_act(3);
            part=idces_act(4);
            
            try
                % Note: In the following, it is always
                % output_parfor{idx_count,1} instead of
                % output_parfor{idx_count,outfun}, because
                % the first dimension will be deleted after
                % the loop in order to free memory.
                if isa(output_parfor_act{idx_count_temp},'double') && ...
                        size(output_parfor_act{idx_count_temp},1)==no_focals(group,subset,trial) && ...
                        size(output_parfor_act{idx_count_temp},2)==no_neighbors_act(group,subset,trial)
                    for ff=1:no_focals(group,subset,trial)
                        for nf=1:no_neighbors_act(group,subset,trial)
                            siz_op=size(output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:));
                            siz_op=siz_op(3:end);
                            if isempty(siz_op) && length(siz_op)>1
                                output_new{group,subset,trial,part,ff,nf}=reshape(output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:),....
                                    siz_op);
                            else
                                output_new{group,subset,trial,part,ff,nf}=squeeze(output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:));
                            end
                        end
                    end
                    
                    % The following is true for maps (cell-of-cells):
                elseif isa(output_parfor_act{idx_count_temp},'cell') && ...
                        size(output_parfor_act{idx_count_temp},1)==no_focals(group,subset,trial) && ...
                        size(output_parfor_act{idx_count_temp},2)==no_neighbors_act(group,subset,trial)
                    for ff=1:no_focals(group,subset,trial)
                        for nf=1:no_neighbors_act(group,subset,trial)
                            siz_op=size(output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:));
                            siz_op=siz_op(3:end);
                            if ~isempty(siz_op)
                                output_new{group,subset,trial,part,ff,nf}=reshape(output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:),....
                                    [siz_op 1]); % ..1] in case length(siz_op)==1
                            else
                                output_new{group,subset,trial,part,ff,nf}=output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:);
                            end
                        end
                    end
                    % The following is true for 'index lists' (cell-of-vectors):
                elseif isa(output_parfor_act{idx_count_temp},'cell') && ...
                        all(size(output_parfor_act{idx_count_temp})==1) && ...
                        size(output_parfor_act{idx_count_temp}{1},1)==no_focals(group,subset,trial) && ...
                        size(output_parfor_act{idx_count_temp}{1},2)==no_neighbors_act(group,subset,trial)
                    for ff=1:no_focals(group,subset,trial)
                        for nf=1:no_neighbors_act(group,subset,trial)
                            siz_op=size(output_parfor_act{idx_count_temp}{1}(ff,nf,:,:,:,:,:,:));
                            siz_op=siz_op(3:end);
                            output_new{group,subset,trial,part,ff,nf}=reshape(output_parfor_act{idx_count_temp}{1}(ff,nf,:,:,:,:,:,:),....
                                [siz_op 1]); % ..1] in case length(siz_op)==1
                        end
                    end
                    
                    
                else
                    error([mfilename ': BadOutput'],[mfilename ': Function ' name_outfuncs{outfun} ' produces wrong output format.'])
                end
            catch
                keyboard
            end
            
        end
        %         output_parfor(:,1) = [];
        
        % If output is very big, make a hard disk copy and load it when needed in order
        % to save some memory.
        mem_size=whos('output_new');
        user = memory;
        
        if mem_size.bytes/user.MemAvailableAllArrays>0
            
            if isfield(input_data(1,1,1),name_outfuncs{outfun}) && ...
                    isfield(input_data(1,1,1).(name_outfuncs{outfun}),'output') && ...
                    ~isempty(input_data(1,1,1).(name_outfuncs{outfun}).output) && ...
                    ischar(input_data(1,1,1).(name_outfuncs{outfun}).output)
                %                 if exist(input_data(1,1,1).(name_outfuncs{outfun}).output,'file')==2
                %                     delete(input_data(1,1,1).(name_outfuncs{outfun}).output);
                %                 end
                input_data(1,1,1).(name_outfuncs{outfun}).output=[];
            end
            strpath=[temp_savepath_root name_outfuncs{outfun} datestr(now, 30) '.mat'];
            
            try
                input_data(1,1,1).(name_outfuncs{outfun}).output=strpath;
                disp(['Buffering output to ' strpath ' ...'])
                
                tic;idSocial_auxiliaries_save(strpath,'output_new'); toc
                % Check if file was saved
                if exist(strpath,'file')==2
                    disp('Done.')
                else
                    warning([mfilename ': Buffering failed.'])
                end
            catch
                delete(strpath)
                warning([mfilename ': Buffering to hard disk failed. Using memory to store variable. ' ...
                    strpath 'deleted.'])
                input_data(1,1,1).(name_outfuncs{outfun}).output=output_new;
            end
            
        else
            input_data(1,1,1).(name_outfuncs{outfun}).output=output_new;
            
        end
        clear output_new
        
    end
    clear argin_scaled output_parfor;
    for temp_save=1 : no_temp_saves
        for of = 1:no_outfuncs
            if exist([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'file')==2
                delete([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'])
            end
        end
    end
    %%  The same for random data
    
    if randomized_calcs
        %%%%%%%%%%%%%%%%%%%
        % Execute function once to estimate output size:
        idx_count=min(2,no_combs);
        limit_frames_act = limit_frames(idx_count,:);
        rand_argin_act=rand_argin_scaled{idx_count};
        if ischar(rand_argin_scaled{idx_count}{trajargin}) %&& ischar(argin_scaled{idx_count}{trajargin}) % && exist(rand_argin_scaled{idx_count}{trajargin},'file')==2
            rand_tr=load(rand_argin_scaled{idx_count}{trajargin});
            rand_tr=rand_tr.rand_tr(limit_frames_act(1):min(size(rand_tr.rand_tr,1),limit_frames_act(2)),:,:,:);
            %             rand_tr=rand_tr(:,:,:,:);
            if applyfilter
                rand_tr = idSocial_filters3D(rand_tr,...
                    options,...
                    info_bodylength(idx_count),...
                    info_framerate(idx_count),...
                    info_circroi(idx_count,:));
            end
            
            if order_nb_function(idx_count)
                % I cannot load the already
                % existing ordered trajectory if I want
                % to apply filters before!
                rand_tr=idSocial_auxiliaries_nearestneighbours(rand_tr);
            end
            
            rand_argin_act{trajargin}=rand_tr;
            
            rand_alloutput=cell(no_funcs,1);
            [rand_alloutput{:}]=feval(function_handle,rand_argin_act{:});
            output_estimate=whos('rand_alloutput');
            output_estimate = output_estimate.bytes*9.53674e-7;
            no_temp_saves = ceil(no_combs*output_estimate/mb_per_temp_save);
            no_combs_per_temp_save = ceil(no_combs/no_temp_saves);
        end
        
        %%%%%%%%%%%%%%%%%
        if use_remote
            for rm = 1:numel(temp_savepath_slave_full)
                dummy = 'Dummy';
                savefast([temp_savepath_root 'REMOTE_' num2str(rm) '_READY.mat'],'dummy');
            end
        end
        
        
        for temp_save = 1 : no_temp_saves
            rand_output_parfor=cell(no_combs_per_temp_save,sum(good_idx));
            temp_idx = (temp_save-1)*no_combs_per_temp_save+1 : min(temp_save*no_combs_per_temp_save,no_combs);
            idx_combs_temp = idx_combs(temp_idx,:);
            no_frames_temp = info_no_frames(temp_idx,:);
            no_focals_temp = info_no_focals(temp_idx,:);
            info_no_frames_part_array_temp = info_no_frames_part_array(temp_idx,:);
            info_framerate_temp = info_framerate(temp_idx);
            rand_argin_scaled_temp = rand_argin_scaled(temp_idx);
            limit_frames_temp = limit_frames(temp_idx,:);
            info_bodylength_temp = info_bodylength(temp_idx);
            order_nb_function_temp = order_nb_function(temp_idx);
            idx_conversion(temp_idx,1)=temp_idx;
            idx_conversion(temp_idx,2)=1:min(no_combs_per_temp_save,numel(temp_idx));
            idx_conversion(temp_idx,3)=temp_save;
            
            %             use_remote = false;
            if ~use_remote || no_temp_saves==1
                %par...
                for idx_count=  1:numel(temp_idx)
                    idces_act=idx_combs_temp(idx_count,:);
                    limit_frames_act = limit_frames_temp(idx_count,:);
                    group=idces_act(1);
                    subset=idces_act(2);
                    trial=idces_act(3);
                    part=idces_act(4);
                    no_frames=no_frames_temp(idx_count,:);
                    
                    idx1=max((part-2)*floor(info_no_frames_part_array_temp(idx_count,:)/2)+1,1);
                    idx2=min(floor(info_no_frames_part_array_temp(idx_count,:)/2)*part+1,no_frames);
                    if idx2>idx1
                        fprintf('Random, Group %i, subset %i, trial %i, time %.1f-%.1f min\n',group,subset,trial,idx1/info_framerate_temp(idx_count)/60,idx2/info_framerate_temp(idx_count)/60);
                        
                        
                        rand_argin_act=rand_argin_scaled_temp{idx_count};
                        if ischar(rand_argin_scaled_temp{idx_count}{trajargin}) %&& ischar(argin_scaled{idx_count}{trajargin}) % && exist(rand_argin_scaled{idx_count}{trajargin},'file')==2
                            %                     tr=load(argin_scaled{idx_count}{trajargin});
                            
                            rand_tr=load(rand_argin_scaled_temp{idx_count}{trajargin});
                            try
                                rand_tr=rand_tr.rand_tr(limit_frames_act(1):min(size(rand_tr.rand_tr,1),limit_frames_act(2)),:,:,:);
                            catch
                                keyboard
                            end
                            rand_tr=rand_tr(idx1:idx2,:,:,:);
                            % Skip/mark pairs 'random focal - random
                            % neighbor'. In any core function, we can skip
                            % calculation for these pairs.
                            try
                                rand_tr(:,no_focals_temp+1:end,no_focals_temp+1:end,1) = NaN;
                                rand_tr(:,no_focals_temp+1:end,no_focals_temp+1:end,2) = Inf;
                            catch
                                keyboard
                            end
                            if applyfilter
                                rand_tr = idSocial_filters3D(rand_tr,...
                                    options,...
                                    info_bodylength_temp(idx_count),...
                                    info_framerate_temp(idx_count),...
                                    info_circroi(idx_count,:));
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
                            
                            if order_nb_function_temp(idx_count)
                                alloutput_new = rand_alloutput;
                                for outc = 1:size(rand_alloutput,1)
                                    for ff = 1:size(rand_alloutput{outc},1)
                                        for nf = 1:size(rand_alloutput{outc},2)
                                            if nf<ff
                                                alloutput_new{outc}(ff,nf,:,:,:,:,:) = rand_alloutput{outc}(ff,nf,:,:,:,:,:);
                                            elseif nf>=ff && nf<size(rand_alloutput{outc},2)
                                                alloutput_new{outc}(ff,nf,:,:,:,:,:) = rand_alloutput{outc}(ff,nf+1,:,:,:,:,:);
                                            elseif  nf==size(rand_alloutput{outc},2)
                                                alloutput_new{outc}(ff,nf,:,:,:,:,:) = rand_alloutput{outc}(ff,ff,:,:,:,:,:);
                                            end
                                            
                                        end
                                    end
                                end
                                rand_alloutput = alloutput_new;
                                %                                 clear alloutput_new
                            end
                            
                            rand_output_parfor(idx_count,:)=rand_alloutput(good_idx);
                            
                        end
                    end
                    
                end
                
                if no_temp_saves>1
                    ofbufferTic=tic;
                    for of = 1:no_outfuncs
                        disp(['Buffering ' 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat' ' ...'])
                        
                        tempOutput =  rand_output_parfor(:,of);
                        savefast([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'tempOutput');
                    end
                    
                    disp(['Done (' num2str(toc(ofbufferTic),'%.1f') 's)'])
                    clear rand_output_parfor tempOutput
                end
                %% Remote goes here:
            else % Remote
                paramstruct.idx_combs_temp = idx_combs_temp;
                paramstruct.no_frames_temp = no_frames_temp;
                paramstruct.info_bodylength_temp = info_bodylength_temp;
                paramstruct.info_no_frames_part_array_temp = info_no_frames_part_array_temp;
                paramstruct.rand_argin_scaled_temp = rand_argin_scaled_temp;
                paramstruct.temp_save = temp_save;
                paramstruct.function_handle = function_handle;
                paramstruct.slave_output_folder = temp_savepath_slave_root;
                paramstruct.master_output_folder = temp_savepath_root;
                paramstruct.no_temp_saves = no_temp_saves;
                paramstruct.temp_idx = temp_idx;
                paramstruct.info_framerate_temp = info_framerate(temp_idx);
                paramstruct.trajargin=trajargin;
                paramstruct.applyfilter=applyfilter;
                paramstruct.order_nb_function_temp = order_nb_function_temp;
                paramstruct.no_funcs = no_funcs;
                paramstruct.no_outfuncs = no_outfuncs;
                paramstruct.good_idx = good_idx;
                paramstruct.no_combs_per_temp_save = no_combs_per_temp_save;
                
                act_remote = mod(temp_save-1,no_remotes)+1;
                
                temp_save_delivered = false;
                
                while ~temp_save_delivered
                    if exist([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'],'file')==2
                        disp(['Starting round ' num2str(temp_save) ' on Computer ' num2str(act_remote)])
                        delete([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'])
                        paramstruct.slave_id = act_remote;
                        parampath_save = [temp_savepath_slave_root{1} 'inputparams_no_' num2str(temp_save) '.mat'];
                        savefast(parampath_save,'paramstruct');
                        parampath_load = [temp_savepath_slave_root{act_remote} 'inputparams_no_' num2str(temp_save) '.mat'];
                        if act_remote == 1 % Local
                            %                         system([' "C:\Program Files\Matlab\R2012a\bin\matlab.exe" "-nodisplay" "-nosplash" "-nodesktop" "-r" "idSocial_function_wrapper_ExecParfor ' parampath '"']);
                            system([' "C:\Program Files\Matlab\R2012a\bin\matlab.exe"  -nodisplay -nosplash -nodesktop -minimize -r "idSocial_function_wrapper_ExecParforRAND ' parampath_load '";exit;']);
                            
                        else
                            %                         keyboard
                            system(['"C:\Program Files (x86)\PSTools\psexec" ' remote_computers{act_remote} ' "C:\Program Files\Matlab\R2012a\bin\matlab.exe"  "-nodisplay" "-nosplash" "-nodesktop" "-minimize" "-r" "idSocial_function_wrapper_ExecParforRAND ' parampath_load ' exit;" ']);
                        end
                        
                        temp_save_delivered = true;
                    else
                        act_remote = mod(act_remote + 1,no_remotes)+1;
                    end
                end
                
            end
            %%
        end
        
        %         if use_parallel; matlabpool close; end
    end
    
    
    all_done=false;
    done_array = false(no_temp_saves,no_outfuncs);
    if use_remote && no_temp_saves>1
        while ~all_done
            for temp_save = 1 : no_temp_saves
                for of = 1:no_outfuncs
                    if exist([temp_savepath_slave_root{1} 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'file')==2
                        done_array(temp_save,of)=true;
                    end
                end
            end
            pause(.01)
            all_done = all(done_array(:));
        end
    end
    
    for act_remote = 1: no_remotes
        if exist([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'],'file')==2
            delete([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'])
        end
    end
    
    if  ~isempty(rand_argin_scaled{1}) && ~ischar(rand_argin_scaled{1}{trajargin}) %&& ~ischar(argin_scaled{1}{trajargin})
        warning([mfilename ': No randomized trajectories found. Continue without.'])
        randomized_calcs=false;
    end
    clear rand_argin_scaled
    if randomized_calcs
        for outfun=1:no_outfuncs
            output_new=cell(no_groups,max(no_subsets),no_trials,max_no_parts,max_no_focals,max_no_neighborsRAND);
            
            
            
            temp_save = 0;
            for idx_count=1:no_combs
                
                temp_save_act = idx_conversion(idx_count,3);
                idx_count_temp = idx_conversion(idx_count,2);
                
                
                if temp_save_act ~= temp_save && no_temp_saves>1
                    disp(['Loading temporary buffer ' 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat' ' ...'])
                    %                 keyboard
                    load([temp_savepath_root 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat'])
                    rand_output_parfor_act = tempOutput;
                    disp('Done.')
                    temp_save = temp_save_act;
                end
                if no_temp_saves==1
                    rand_output_parfor_act = rand_output_parfor(:,outfun);
                end
                
                
                idces_act=idx_combs(idx_count,:);
                group=idces_act(1);
                subset=idces_act(2);
                trial=idces_act(3);
                part=idces_act(4);
                
                try
                    if isa(rand_output_parfor_act{idx_count_temp},'double') && ...
                            size(rand_output_parfor_act{idx_count_temp},1)==no_focalsRAND(group,subset,trial) && ...
                            size(rand_output_parfor_act{idx_count_temp},2)==no_neighborsRAND(group,subset,trial)
                        for ff=1:no_focals(group,subset,trial)
                            for nf=1:no_neighborsRAND(group,subset,trial)
                                siz_op=size(rand_output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:));
                                siz_op=siz_op(3:end);
                                
                                if isempty(siz_op) && length(siz_op)>1
                                    output_new{group,subset,trial,part,ff,nf}=reshape(rand_output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:),....
                                        siz_op);
                                else
                                    output_new{group,subset,trial,part,ff,nf}=squeeze(rand_output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:));
                                end
                            end
                        end
                        
                        
                        % The following is true for maps (cell-of-cells):
                    elseif isa(rand_output_parfor_act{idx_count_temp},'cell') && ...
                            size(rand_output_parfor_act{idx_count_temp},1)==no_focalsRAND(group,subset,trial) && ...
                            size(rand_output_parfor_act{idx_count_temp},2)==no_neighborsRAND(group,subset,trial)
                        for ff=1:no_focals(group,subset,trial)
                            for nf=1:no_neighborsRAND(group,subset,trial)
                                siz_op=size(rand_output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:));
                                siz_op=siz_op(3:end);
                                
                                output_new{group,subset,trial,part,ff,nf}=reshape(rand_output_parfor_act{idx_count_temp}(ff,nf,:,:,:,:,:,:),....
                                    [siz_op 1]); % ..1] in case length(siz_op)==1
                                
                            end
                        end
                        % The following is true for 'index lists' (cell-of-vectors):
                    elseif isa(rand_output_parfor_act{idx_count_temp},'cell') && ...
                            all(size(rand_output_parfor_act{idx_count_temp})==1) && ...
                            size(rand_output_parfor_act{idx_count_temp}{1},1)==no_focalsRAND(group,subset,trial) && ...
                            size(rand_output_parfor_act{idx_count_temp}{1},2)==no_neighborsRAND(group,subset,trial)
                        for ff=1:no_focals(group,subset,trial)
                            for nf=1:no_neighborsRAND(group,subset,trial)
                                siz_op=size(rand_output_parfor_act{idx_count_temp}{1}(ff,nf,:,:,:,:,:,:));
                                siz_op=siz_op(3:end);
                                
                                output_new{group,subset,trial,part,ff,nf}=reshape(rand_output_parfor_act{idx_count_temp}{1}(ff,nf,:,:,:,:,:,:),....
                                    [siz_op 1]); % ..1] in case length(siz_op)==1
                                
                            end
                        end
                        
                        
                    else
                        error([mfilename ': BadOutput'],[mfilename ': Function ' name_outfuncs{outfun} ' produces wrong output format.'])
                    end
                catch
                    keyboard
                end
            end
            %             rand_output_parfor(:,1)=[];
            % If output is very big, make a hard disk copy and load it when needed in order
            % to save some memory.
            mem_size=whos('output_new');
            user = memory;
            
            if mem_size.bytes/user.MemAvailableAllArrays>0
                if exist(temp_savepath_root,'dir')~=7
                    mkdir(temp_savepath_root)
                end
                
                
                
                % For random:
                
                if isfield(input_data(1,1,1),[name_outfuncs{outfun} 'RAND']) && ...
                        isfield(input_data(1,1,1).([name_outfuncs{outfun} 'RAND']),'output') && ...
                        ~isempty(input_data(1,1,1).([name_outfuncs{outfun} 'RAND'])) && ...
                        ischar(input_data(1,1,1).([name_outfuncs{outfun} 'RAND']).output)
                    delete(input_data(1,1,1).([name_outfuncs{outfun} 'RAND']).output);
                    input_data(1,1,1).([name_outfuncs{outfun} 'RAND']).output=[];
                end
                strpath=[temp_savepath_root name_outfuncs{outfun} 'RAND' datestr(now, 30) '.mat'];
                
                try
                    input_data(1,1,1).([name_outfuncs{outfun} 'RAND']).output=strpath;
                    disp(['Buffering output to ' strpath ' ...'])
                    idSocial_auxiliaries_save(strpath,'output_new')
                    % Check if file was saved
                    if exist(strpath,'file')==2
                        disp('Done.')
                    else
                        warning([mfilename ': Buffering failed.'])
                    end
                    clear output_new;
                    
                catch
                    delete(strpath)
                    warning([mfilename ': Buffering to hard disk failed. Using memory to store variable. ' ...
                        strpath 'deleted.'])
                    input_data(1,1,1).([name_outfuncs{outfun} 'RAND']).output=output_new;
                end
                
            else
                input_data(1,1,1).([name_outfuncs{outfun} 'RAND']).output=output_new;
                
            end
            
            
        end
        clear rand_output_parfor
        for temp_save=1 : no_temp_saves
            for of = 1:no_outfuncs
                if exist([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'file')==2
                    delete([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'])
                end
            end
        end
    end
    try
        if ~verLessThan('matlab', '8.3.0')
            p = gcp('nocreate');
            if isempty(p)
                no_cores_open = 0;
            else
                no_cores_open = p.NumWorkers;
            end
            if use_parallel  && no_combs>min_number_loop && no_cores_open~=0;
                delete(gcp('nocreate'));
                
            end
        else
            no_cores_open=matlabpool('size');
            if use_parallel  && no_combs>min_number_loop && no_cores_open~=0;
                matlabpool close;
            end
        end
    catch
    end
end
input_data(1,1,1).(act_method).plot_mode=plot_mode;
if actualize_plot_mode || manual_actualize_plot_mode
    
    if isfield(input_data(1,1,1),[act_method '4Statistics'])
        input_data(1,1,1).(act_method).output_plot=idSocial_auxiliaries_reduceCellArray(input_data(1,1,1).(act_method).output,...
            plot_mode,input_data(1,1,1).([act_method '4Statistics']).output);
    else
        input_data(1,1,1).(act_method).output_plot=idSocial_auxiliaries_reduceCellArray(input_data(1,1,1).(act_method).output,...
            plot_mode);
    end
end


if randomized_calcs || isfield(input_data(1,1,1),[act_method 'RAND'])
    if actualize_plot_mode || manual_actualize_plot_mode
        if isfield(input_data(1,1,1),[act_method '4StatisticsRAND'])
            input_data(1,1,1).(act_method).output_plotRAND=idSocial_auxiliaries_reduceCellArray(input_data(1,1,1).([act_method 'RAND']).output,...
                plot_mode,input_data(1,1,1).([act_method '4StatisticsRAND']).output);
        else
            input_data(1,1,1).(act_method).output_plotRAND=idSocial_auxiliaries_reduceCellArray(input_data(1,1,1).([act_method 'RAND']).output,...
                plot_mode);
        end
        
        if isfield(input_data(1,1,1).(act_method),'output_plot') && ~isempty(input_data(1,1,1).(act_method).output_plot)
            pm=input_data(1,1,1).(act_method).output_plot;
            % Merge plot_mode
            no_entries=size(input_data(1,1,1).(act_method).output_plotRAND.data_string,2);
            input_data(1,1,1).(act_method).output_plotRAND.data_string=...
                cellfun(@(x) strrep(x,x,[char(x(1)+no_entries) x(2:end)]),input_data(1,1,1).(act_method).output_plotRAND.data_string,'UniformOutput',false);
            input_data(1,1,1).(act_method).output_plotRAND.legendstring=...
                cellfun(@(x) strrep(x,x,[x ', rand.']),input_data(1,1,1).(act_method).output_plotRAND.legendstring,'UniformOutput',false);
            
            
            pmrand=input_data(1,1,1).(act_method).output_plotRAND;
            pm.data=[pm.data pmrand.data];
            pm.data_sign=[pm.data_sign pmrand.data_sign];
            pm.data_signMedian=[pm.data_signMedian pmrand.data_signMedian];
            pm.data_Median=[pm.data_Median pmrand.data_Median];
            pm.data_dev=[pm.data_dev pmrand.data_dev];
            pm.statistics_type=[pm.statistics_type pmrand.statistics_type];
            pm.statistics_on_idx=[pm.statistics_on_idx pmrand.statistics_on_idx];
            pm.data_no_datapoints=[pm.data_no_datapoints pmrand.data_no_datapoints];
            pm.data_string=[pm.data_string pmrand.data_string];
            pm.filterstring=[pm.filterstring pmrand.filterstring];
            pm.legendstring=[pm.legendstring pmrand.legendstring];
            pm.dim_names=[pm.dim_names pmrand.dim_names];
            pm.subplotstring=[pm.subplotstring pmrand.subplotstring];
            input_data(1,1,1).(act_method).output_plot=pm;
            input_data(1,1,1).(act_method)=rmfield(input_data(1,1,1).(act_method),'output_plotRAND');
        end
    end
end
if ~isempty(external_func_name)
    input_data(1,1,1).(external_func_name)=input_data(1,1,1).(act_method);
    act_method = external_func_name;
end
% Save input_data
if exist(temp_savepath_root,'dir')~=7
    mkdir(temp_savepath_root)
end

if actualize || manual_actualize || ...
        actualize_plot_mode || manual_actualize_plot_mode
    disp(['Buffering input data for ' act_method '...'])
    % save(temp_savepath,'input_data','-v7.3');
    idSocial_auxiliaries_save(temp_savepath,'input_data');
    disp(['Saved ' temp_savepath])
end
diary off
end