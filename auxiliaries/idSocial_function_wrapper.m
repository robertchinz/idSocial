function input_data=idSocial_function_wrapper(input_data, options, def_options_function, plot_mode, functionInfo, external_func_name)
if nargin < 6 || isempty(external_func_name)
    external_func_name = '';
end

try
    [~, progressBox] = findall(0,'Tag','idSocialMessageBox');
catch
    progressBox = [];
end
% plot_mode.display_mode='extern'; % For now, there will no automatic plots any more.
mb_per_temp_save=1000;
bytes_per_temp_save=mb_per_temp_save*1048576;
pmfields = fieldnames(plot_mode);
no_flds = size(pmfields,1);
include_end_parts = false;
for k=1:no_flds
    if iscell(plot_mode.(pmfields{k}))
        chr = cellfun(@(x) ischar(x),plot_mode.(pmfields{k}));
        if any(chr)
            for charidx = find(chr)
                include_end_parts = ~isempty(strfind([plot_mode.(pmfields{k}){charidx}],'Time')) || include_end_parts;
            end
        end
%         if include_end_parts; keyboard; end
    elseif ischar(plot_mode.(pmfields{k}))
        include_end_parts = ~isempty(strfind(plot_mode.(pmfields{k}),'Time')) || include_end_parts;
%          if include_end_parts; keyboard; end
    end
end
skip_end_parts = ~include_end_parts;
skip_end_partsRAND = ~include_end_parts;
% Number of workers for parallel computations:
psize=4;
manual_actualize=false;
manual_actualize_plot_mode = false;%true;

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
def_options.project_name='';
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
def_options.filter_neighbor_circularROI = []; 
def_options.filter_AllMembersPresent = -inf; def_options.filter_AllMembersPresentWorstFocal = -inf; 
def_options.filter_AllMembersPresentWorstFocal = -inf; 
def_options.filter_WorstIndividual = -inf;
def_options.filter_CustomFrames = [];
def_options.filter_focal_spatial_sectors = []; %'lateral' or 'frontal'
% Options for maps:
def_options.rowlabelfontsize=12;
% Options for territories:
def_options.discard_bottom_percent=10;
def_options.plot_bootstrap_repetitions=10000;
def_options.plot_bootstrap_statfunc='Mean';
def_options.plot_deviation=true;
def_options.temp_savepath='';
def_options.temp_savepath_slave='';
def_options.order_neighbors = false;
def_options.order_com = false;
def_options.internal_random_controls = [];
def_options.parallel_processing = true;

% DELETE THIS
def_options.filter_filter_WorstIndividual = def_options.filter_WorstIndividual;

if ~isfield(options,'parallel_processing')
    options.parallel_processing=true;
end
if ~isfield(options,'order_com')
    options.order_com=false;
end
%%% END DELETE THIS

% Get save path
if ~isfield(options,'temp_savepath_slave')
    options.temp_savepath_slave='';
end




% [~, project_name]= idSocial_recursiveGetOptionsFromOptionsCell(input_data.options,'project_name');
% if isempty(project_name) %%&& strcmpi(project_name(1),'<') && strcmpi(project_name(end),'>')
%     project_name=['NewProject_' datestr(now, 30)];
% end
% if ~isempty(project_name) && strcmpi(project_name,'New Project') 
%     project_name=['NewProject_' datestr(now, 30)];
% end

[~, temp_savepath]= idSocial_recursiveGetOptionsFromOptionsCell(input_data.options,'temp_savepath');
options.temp_savepath = temp_savepath;
% if isfield(options,'temp_savepath') && ~isempty(options.temp_savepath) && ~(strcmpi(options.temp_savepath(1),'<') && strcmpi(options.temp_savepath(end),'>'))
%     temp_savepath = options.temp_savepath;
% end
% if ~strcmpi(temp_savepath(end-3:end),'.mat') % temp_savepath is directory
if ~strcmp(temp_savepath(end),filesep)
    temp_savepath = [temp_savepath filesep];
end
temp_savepath_root = temp_savepath;
temp_savepath = [temp_savepath 'input_data.mat'];
% else % temp_savepath is file
%     delidx=strfind(temp_savepath,filesep);
%     temp_savepath_root = temp_savepath(1:max(delidx));
% end
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
            if ~strcmp(temp_savepath_slave{rm}(end),filesep)
                temp_savepath_slave{rm} = [temp_savepath_slave{rm} filesep];
            end
            temp_savepath_slave_root{rm} = temp_savepath_slave{rm};
            temp_savepath_slave{rm} = [temp_savepath_slave{rm} 'input_data.mat'];
        else % temp_savepath is file
            delidx=strfind(temp_savepath_slave{rm},filesep);
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
logfile = [temp_savepath_root 'log_' datestr(now, 30) '.txt'];
diary(logfile)

def_options_fct_fldnames=fieldnames(def_options_function);
for k=1:size(def_options_fct_fldnames,1)
    def_options.(def_options_fct_fldnames{k})=def_options_function.(def_options_fct_fldnames{k});
end

act_method=def_options.act_method;
core_act_method=strrep(func2str(functionInfo.handle),'idSocial_','');

if strcmp(act_method(end-10:end), '4Statistics')
    principal_method = act_method(1:end-11);
else
    principal_method = act_method;
end
already_exists =  isfield(input_data,principal_method);
set_already_calculated = -1;
act_datestr = ['id' datestr(now, 30)];

[input_data, options]=idSocial_readparams(input_data,options,def_options,act_method);

[~, filter_focal_list]= idSocial_recursiveGetOptionsFromOptionsCell(input_data.options,'filter_focal_list');
[~, filter_neighbor_list]= idSocial_recursiveGetOptionsFromOptionsCell(input_data.options,'filter_neighbor_list');

if ~isempty(filter_focal_list) && isfield(options,'filter_focal_list')
    options.filter_focal_list = unique([filter_focal_list,options.filter_focal_list]);
else
    options.filter_focal_list =filter_focal_list;
end
if ~isempty(filter_neighbor_list)&& isfield(options,'filter_neighbor_list')
    options.filter_neighbor_list = unique([filter_focal_list,options.filter_neighbor_list]);
    else
    options.filter_neighbor_list =filter_neighbor_list;
end
% DELETE THIS
% if isfield(options,'filter_filter_WorstIndividual') && ~isfield(options,'filter_WorstIndividual')
%     options.filter_WorstIndividual = options.filter_filter_WorstIndividual;
% end
% if isfield(options,'filter_WorstIndividual') && ~isfield(options,'filter_filter_WorstIndividual')
%     options.filter_filter_WorstIndividual = options.filter_WorstIndividual;
% end

%%% END DELETE THIS

order_nb_function_single = options.order_neighbors;
order_com_function_single = options.order_com;

function_option_fields=intersect(fieldnames(input_data.(act_method).options),fieldnames(def_options_function));
use_parallel = options.parallel_processing;
use_parallel = false;
if strcmpi(getenv('COMPUTERNAME'),'INRC-GPOLA-08') % Raul
    use_parallel = false;
end

plotfields=strfind(function_option_fields,'plot');
plotfields=cellfun(@(x)~isempty(x),plotfields);
savepathfield = strcmpi(function_option_fields,'temp_savepath');
function_option_fields(ismember(function_option_fields,'act_method') | plotfields | savepathfield  )=[];
options_vals=cellfun(@(x) input_data.(act_method).options.(x),function_option_fields,'UniformOutput',false);
% keyboard

% Get act_set:

if ~already_exists || ~isfield(input_data.(principal_method),'Set1') || isempty(input_data.(principal_method).Set1) || ~isfield(input_data.(principal_method).Set1,'output_plot')
    act_set = 1;
    maxSet = 1;
else
    fn = fieldnames(input_data.(principal_method));
    maxSet = -1;
    for k=1:size(fn,1)
        if ~isempty(strfind(fn{k},'Set'))
            maxSet = max(maxSet,str2double(strrep(fn{k},'Set','')));
        end
    end
    act_set = maxSet+1;
end



if isfield(input_data,core_act_method) && isfield(input_data.(core_act_method),'options') && ...
        already_exists
    %
    %
    %     % Check if base data has changed:
    %     new_data = false;
    %     if isfield(input_data.(core_act_method),'trajectory_ver') && ...
    %     end
    %     %%%%%%%
    %
    %
    % DELETE THIS!!!
%     if ~isfield(input_data.(core_act_method).options,'filter_WorstIndividual') && isfield(input_data.(core_act_method).options,'filter_filter_WorstIndividual')
%         input_data.(core_act_method).options.filter_WorstIndividual = input_data.(core_act_method).options.filter_filter_WorstIndividual;
%         disp('DELETE THIS!')
%         if ~isfield(input_data.(core_act_method).options,'filter_AllMembersPresentWorstFocal')
%             input_data.(core_act_method).options.filter_AllMembersPresentWorstFocal = -inf;
%         end
%     end
%     if ~isfield(input_data.(core_act_method).options,'significance_between_groups')
%         disp([mfilename ' Line 243: DELETE THIS!'])
%         
%         input_data.(core_act_method).options.significance_between_groups = true;
%         
%     end
%     if ~isfield(input_data.(core_act_method).options,'parallel_processing')
%         disp([mfilename ' Line 243: DELETE THIS!'])
%         
%         input_data.(core_act_method).options.parallel_processing = true;
%         
%     end
%     if ~isfield(input_data.(core_act_method).options,'filter_CustomFrames')
%         disp([mfilename ' Line 243: DELETE THIS!'])
%         
%         input_data.(core_act_method).options.filter_CustomFrames = [];
%         
%     end
%     if ~isfield(input_data.(core_act_method).options,'filter_focal_accelerationlimits_bl_per_s2')
%         disp([mfilename ' Line 243: DELETE THIS!'])
%         input_data.(core_act_method).options.filter_focal_accelerationlimits_bl_per_s2 = [];
%     end
%     if ~isfield(input_data.(core_act_method).options,'filter_neighbor_accelerationlimits_bl_per_s2')
%         disp([mfilename ' Line 243: DELETE THIS!'])
%         input_data.(core_act_method).options.filter_neighbor_accelerationlimits_bl_per_s2 = [];
%     end
%     if ~isfield(input_data.(core_act_method).options,'filter_focal_spatial_sectors')
%         disp([mfilename ' Line 243: DELETE THIS!'])
%         
%         input_data.(core_act_method).options.filter_focal_spatial_sectors = [];
%         
%     end
%      if ~isfield(input_data.(core_act_method).options,'filter_framesWithAllNeighborsOnly')
%         disp([mfilename ' Line 280: DELETE THIS!'])
%         
%         input_data.(core_act_method).options.filter_framesWithAllNeighborsOnly = false;
%         
%      end
%      if ~isfield(input_data.(core_act_method).options,'temp_savepath_slave')
%         disp([mfilename ' Line 280: DELETE THIS!'])
%         
%         input_data.(core_act_method).options.temp_savepath_slave = [];
%         
%     end
%     if ~isfield(input_data.(core_act_method).options,'order_com')
%         disp([mfilename ' Line 280: DELETE THIS!'])
%         
%         input_data.(core_act_method).options.order_com = false;
%         
%     end
    %%% END
%     actualize = false;
%     for nset = 1:maxSet
%         if isfield(input_data.(act_method),['Set' num2str(nset)]) && isfield(input_data.(act_method).(['Set' num2str(nset)]),'options')
%             func_opts_valsSet = cellfun(@(x) input_data.(act_method).(['Set' num2str(nset)]).options.(x),function_option_fields,'UniformOutput',false);
%                     
%             changed_opts = cellfun(@(x,y) ~isequal(x,y),options_vals,func_opts_valsSet);
%             if ~any(changed_opts)
%                 set_already_calculated = nset;
%             end
%         end
%     end
%     if set_already_calculated<0
%         actualize = true;
%     end
    
    actualize = false;
    if isfield(input_data.(act_method),'output') && ~isempty(input_data.(act_method).output)
        outfn = fieldnames(input_data.(act_method).output);
        for nset = 1:numel(outfn)
            func_opts_valsSet = cellfun(@(x) input_data.(act_method).output.(outfn{nset}).options.(x),function_option_fields,'UniformOutput',false);
            
            changed_opts = cellfun(@(x,y) ~isequal(x,y),options_vals,func_opts_valsSet);
            if ~any(changed_opts)
                set_already_calculated = nset;
                act_datestr = outfn{set_already_calculated};
            end
        end
    end
    if set_already_calculated<0
        actualize = true;
%         act_datestr = ['id' datestr(now, 30)];
    end
    
    func_opts_vals=cellfun(@(x) input_data.(core_act_method).options.(x),function_option_fields,'UniformOutput',false);
%     actualize=~isequal(options_vals,func_opts_vals);
    if actualize
        idSocial_auxiliaries_message('Previous data found. The following options have changed:',progressBox);
        
        changed_opts = cellfun(@(x,y) ~isequal(x,y),options_vals,func_opts_vals);
        %         no_changed = sum(changed_opts);
        for k=find(changed_opts)'
            idSocial_auxiliaries_message([upper(function_option_fields{k}) ': '],progressBox);
            idSocial_auxiliaries_message('Old value: ',progressBox); display(func_opts_vals{k});
            idSocial_auxiliaries_message('New value: ',progressBox); display(options_vals{k});
            

        end
        
    end
else
    if ~isfield(input_data,core_act_method)
        idSocial_auxiliaries_message('No previous data found. Start calculations...',progressBox);
        idSocial_auxiliaries_message('No previous data found. Start calculations...',progressBox);
    end
    actualize=1;
    
end
input_data.(act_method).output.(act_datestr).options = input_data.(act_method).options;
%     if actualize==1; keyboard; end



actualize_plot_mode=1;


%% Get info
info=               input_data.info;
% blpxl=              info.blpxl;
framerate=          info.framerate;
duration=           info.duration;
no_frames=          info.no_frames;
trials=             info.trials;
no_focals=          info.no_focals;
no_neighbors=          info.no_neighbors;
no_neighborsRAND=          info.no_neighborsRAND;
no_focalsRAND=          info.no_neighborsRAND; % The same as neighbors!
% disp([mfilename ': Line 74, CHANGE!!!'])
no_groups=         info.no_groups;
no_subsets=        info.no_subsets;
no_trials=         max(info.no_trials(:));

if ~isfield(info,'filter_AllMembersPresent')
        idSocial_auxiliaries_message('info.filter_AllMembersPresent not found. Set to inf.',progressBox);
        info.filter_AllMembersPresent=inf(no_groups,no_subsets,no_trials);
end
if ~isfield(info,'filter_WorstIndividual')
    idSocial_auxiliaries_message('info.filter_WorstIndividual not found. Set to inf.',progressBox);
    
    info.filter_WorstIndividual=inf(no_groups,no_subsets,no_trials);
end
if ~isfield(options,'filter_AllMembersPresent')
    options.filter_AllMembersPresent=-inf(no_groups,no_subsets,no_trials);
end
if ~isfield(options,'filter_WorstIndividual')
    options.filter_WorstIndividual=-inf(no_groups,no_subsets,no_trials);
end
if ~isfield(options,'filter_AllMembersPresentWorstFocal')
        options.filter_AllMembersPresentWorstFocal=-inf(no_groups,no_subsets,no_trials);
end
if ~isfield(info,'filter_AllMembersPresentWorstFocal')
    idSocial_auxiliaries_message('info.filter_AllMembersPresentWorstFocal not found. Set to inf.',progressBox);
    info.filter_AllMembersPresentWorstFocal=inf(no_groups,no_subsets,no_trials);
end
%% Some standard parameters
act_method=def_options.act_method;
randomized_calcs=options.random_data;
max_no_focals=max(max(max(no_focals)));
max_no_neighbors=max(max(max(no_neighbors)));
max_no_neighborsRAND=max(max(max(no_neighborsRAND)));


timeintervals_in_min=options.timeintervals_in_min;
if isempty(timeintervals_in_min) || isinf(timeintervals_in_min) % || timeintervals_in_min> min(min(floor(min(no_frames./framerate)/60)));
    timeintervals_in_min=duration;
end;
no_frames_part_array=floor(timeintervals_in_min.*framerate*60)-1;
no_frames_part_arrayRAND=no_frames_part_array;

no_parts=2*floor(no_frames./no_frames_part_array)+1;
no_partsRAND=no_parts;

max_no_parts=max(max(max(no_parts)));
if max_no_parts ==3 && skip_end_parts
    skip_end_parts = true;
    skip_end_partsRAND = true;
end

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

if (all(cellfun(@(x) ischar(x),input_data.movementdata(:))) && ~all(cellfun(@(x) exist(x,'file')==2,input_data.movementdata(:))) && any(trajargin)) % The last condition checks if trajectories are temporarily saved and needed
    error('Could not find files in temporary folder.')
end


filterfields=strfind(fieldnames(options),'filter');
filterfields=cellfun(@(x)~isempty(x),filterfields);
optionfieldnames = fieldnames(options);
applyfilter = any(cellfun(@(x) ~isempty(options.(x)),optionfieldnames(filterfields)));


% Check if there is a custom filter from other function
% output
CustomFilterFrame = [];
if isfield(options,'filter_CustomFrames') && ~isempty(options.filter_CustomFrames)
    if ischar(options.filter_CustomFrames) && isfield(input_data,options.filter_CustomFrames) && ...
            exist(input_data.(options.filter_CustomFrames).output,'file') == 2
        load(input_data.(options.filter_CustomFrames).output);
        CustomFilterFrame = output_new;
        clear output_new;
    end
else
   CustomFilterFrame = -1; 

end

name_outfuncs=functionInfo.output2function(~strcmpi(functionInfo.output2function,'temp') & cellfun(@(x) ~isempty(x),functionInfo.output2function));

if isfield(input_data,core_act_method) && isfield(input_data.(core_act_method),'info') && ~actualize && ~manual_actualize
    funcinfo = input_data.(core_act_method).info;
else
    funcinfo = [];
end

parts_and_slices = -1 * ones(no_groups, max(no_subsets(:)), max(no_trials(:))); % ==0 if parts, but no slices, == 1 if slices, but no parts, == 2 if both parts and slices  
parts_and_slicesRAND = -1 * ones(no_groups, max(no_subsets(:)), max(no_trials(:))); % ==0 if parts, but no slices, == 1 if slices, but no parts, == 2 if both parts and slices  

function_handle = functionInfo.handle;
if  actualize || manual_actualize
    
    idSocial_auxiliaries_message([act_method ': Updating results...'],progressBox);
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
        ismember(char_input_params, strcat('info.',fieldnames(input_data.info))))=true;
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
    input_data.(core_act_method).options=cell2struct(options_vals,function_option_fields);
    
    comb_count=1;
    comb_countRAND=1;

    idx_combs=[];
    idx_combsRAND=[];
    
    for group=1:no_groups
        for subset=1:no_subsets(group)
            try
            if size(trials{group,subset},2)==1; trlist = trials{group,subset}'; else trlist = trials{group,subset}; end
            catch
                keyboard
            end
            for trial= trlist
                if info.filter_AllMembersPresent(group,subset,trial)>=options.filter_AllMembersPresent*100 && ...
                        info.filter_WorstIndividual(group,subset,trial)>=options.filter_WorstIndividual*100 && ...
                        info.filter_AllMembersPresentWorstFocal(group,subset,trial)>=options.filter_AllMembersPresentWorstFocal*100
                    no_frames=info.no_frames(group,subset,trial);
                    no_framesSlice = info.slice_framerange{group,subset,trial}(1,2) - info.slice_framerange{group,subset,trial}(1,2) + 1;
                    if ~isempty(info.slice_framerangeRAND{group,subset,trial})
                        no_framesSliceRAND = info.slice_framerangeRAND{group,subset,trial}(1,2) - info.slice_framerangeRAND{group,subset,trial}(1,2) + 1;
                    else
                        no_framesSliceRAND = 0;
                    end
                    if info.no_slices(group,subset,trial) >1 && no_parts(group,subset,trial) == 3 % No parts, but sliced -> Generate "parts" from slices
                        parts_and_slices(group,subset,trial) = 1;
                        no_parts_act = info.no_slices(group,subset,trial);
                        no_parts(group,subset,trial) = no_parts_act;
                        info.no_parts(group,subset,trial) = no_parts_act;
                        
                        skip_end_parts = false; % Better include
                        input_data.info = info;
                        
                    elseif info.no_slices(group,subset,trial) >1 && no_parts(group,subset,trial) > 1 % Parts and slices -> Convert "slices" into "parts" of the same length as "parts" (if possible)
                        parts_and_slices(group,subset,trial) = 2;
                        if no_frames_part_array(group,subset,trial) > no_framesSlice
                            warning([mfilename ': No. of frames per part exceeds no. of frames per slice.'])
                            disp(['No. of frames per part will be set to ' num2str(no_framesSlice)])
                            no_frames_part_array(group,subset,trial) = no_framesSlice;
                            no_parts(group,subset,trial) = info.no_slices(group,subset,trial);
                            info.no_parts(group,subset,trial) = no_parts(group,subset,trial);
                        end
                        no_parts_act = no_parts(group,subset,trial);
                        skip_end_parts = false;
                        input_data.info = info;
                        
                    elseif info.no_slices(group,subset,trial) == 1
                        parts_and_slices(group,subset,trial) = 0;
                        % Do nothing ?!
                        no_parts_act = no_parts(group,subset,trial);
                    end
                    
                    
                    for part=1:no_parts_act
                        if ~skip_end_parts || (part~=1 && part ~=no_parts(group,subset,trial))
                            if parts_and_slices(group,subset,trial) == 1 % slices, no parts. Do not use overlapping bins, use slices as parts.
                                idx1=info.slice_framerange{group,subset,trial}(part,1);
                                idx2=info.slice_framerange{group,subset,trial}(part,2);
                                if idx2>idx1
                                    idx_combs=vertcat(idx_combs,[group subset trial part]);
                                    comb_count=comb_count+1;
                                end
                            elseif parts_and_slices(group,subset,trial) == 0 % Only parts, use overlapping bins, no elaborate merging of slices necessary.
                                idx1=max((part-2)*floor(no_frames_part_array(group,subset,trial)/2)+1,1);
                                idx2=min(floor(no_frames_part_array(group,subset,trial)/2)*part+1,no_frames);
                                if idx2>idx1
                                    idx_combs=vertcat(idx_combs,[group subset trial part]);
                                    comb_count=comb_count+1;
                                end
                            end
                        end
                    end
                    
                    % Same for RAND
                    if info.no_slicesRAND(group,subset,trial) >1 && no_partsRAND(group,subset,trial) == 3 % No parts, but sliced -> Generate "parts" from slices
                        parts_and_slicesRAND(group,subset,trial) = 1;
                        no_partsRAND_act = info.no_slicesRAND(group,subset,trial);
                        no_partsRAND(group,subset,trial) = no_partsRAND_act;

                        info.no_partsRAND(group,subset,trial) = no_partsRAND_act;
                        skip_end_partsRAND = false; % Better include
                        input_data.info = info;
                        
                    elseif info.no_slicesRAND(group,subset,trial) >1 && no_partsRAND(group,subset,trial) > 1 % Parts and slices -> Convert "slices" into "parts" of the same length as "parts" (if possible)
                        parts_and_slicesRAND(group,subset,trial) = 2;
                        if no_frames_part_arrayRAND(group,subset,trial) > no_framesSliceRAND
                            warning([mfilename ': No. of frames per part for randomization exceeds no. of frames per slice.'])
                            disp(['No. of frames per part will be set to ' num2str(no_framesSliceRAND)])
                            no_frames_part_array(group,subset,trial) = no_framesSliceRAND;
                            no_partsRAND(group,subset,trial) = info.no_slicesRAND(group,subset,trial);
                            info.no_partsRAND(group,subset,trial) = no_partsRAND(group,subset,trial);
                        end
                        no_partsRAND_act = no_partsRAND(group,subset,trial);
                        skip_end_partsRAND = false;
                        input_data.info = info;
                        
                    elseif info.no_slicesRAND(group,subset,trial) == 1
                        parts_and_slicesRAND(group,subset,trial) = 0;
                        % Do nothing ?!
                        no_partsRAND_act = no_partsRAND(group,subset,trial);
                    end

                    for part=1:no_partsRAND_act
                        if ~skip_end_partsRAND || (part~=1 && part ~=no_partsRAND(group,subset,trial))
                            if parts_and_slicesRAND(group,subset,trial) == 1 % slices, no parts. Do not use overlapping bins, use slices as parts.
                                idx1=info.slice_framerangeRAND{group,subset,trial}(part,1);
                                idx2=info.slice_framerangeRAND{group,subset,trial}(part,2);
                                if idx2>idx1
                                    idx_combsRAND=vertcat(idx_combsRAND,[group subset trial part]);
                                    comb_countRAND=comb_countRAND+1;
                                end
                                
                            elseif parts_and_slicesRAND(group,subset,trial) == 0 % Only parts, use overlapping bins, no elaborate merging of slices necessary.
                                idx1=max((part-2)*floor(no_frames_part_array(group,subset,trial)/2)+1,1);
                                idx2=min(floor(no_frames_part_array(group,subset,trial)/2)*part+1,no_frames);
                                if idx2>idx1
                                    idx_combsRAND=vertcat(idx_combsRAND,[group subset trial part]);
                                    comb_count=comb_count+1;
                                end
                            end
                        end
                    end
                    
                    
                    
                end
                if info.filter_AllMembersPresent(group,subset,trial)<options.filter_AllMembersPresent*100
                    idSocial_auxiliaries_message(['Omitting (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') because of ''filter_AllMembersPresent'''],progressBox);
                end
                if info.filter_WorstIndividual(group,subset,trial)<options.filter_WorstIndividual*100
                    idSocial_auxiliaries_message(['Omitting (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') because of ''filter_WorstIndividual'''],progressBox);
                end
                if info.filter_AllMembersPresentWorstFocal(group,subset,trial)<options.filter_AllMembersPresentWorstFocal*100
                    idSocial_auxiliaries_message(['Omitting (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') because of ''filter_AllMembersPresentWorstFocal'''],progressBox);
                end
            end
        end
    end
    max_no_parts=max(max(max(no_parts))); % Update max_no_parts in case we use slices.
    
    no_combs=size(idx_combs,1);
    info_no_frames=NaN(no_combs,1);
    info_framerate=NaN(no_combs,1);
    info_bodylength=NaN(no_combs,1);
    info_circroi=NaN(no_combs,4);
    info_no_frames_part_array=NaN(no_combs,1);
    no_combsRAND=size(idx_combsRAND,1);
    info_no_frames_part_arrayRAND=NaN(no_combsRAND,1);

    info_no_focals=NaN(no_combs,1);
    order_nb_trial=NaN(no_combs,1);
    order_nb_function=NaN(no_combs,1);
    order_com_function=NaN(no_combs,1);

    %     function_handle=cell(no_combs,1);
    argin_scaled=cell(no_combs,1);
    rand_argin_scaled=cell(no_combs,1);
    CustomFilterCell = cell(no_combs,1);
%     comb_count=1;
    for group=1:no_groups
        for subset=1:no_subsets(group)
            if size(trials{group,subset},2)==1; trlist = trials{group,subset}'; else trlist = trials{group,subset}; end
            for trial= trlist
                
                no_frames=info.no_frames(group,subset,trial);
                
                % Parts for real data
                for part=1:no_parts(group,subset,trial)
                    [lg,ridx] = ismember([group,subset,trial,part],idx_combs,'rows');
                    if lg
                        comb_count = ridx;
                        if ~skip_end_parts || (part~=1 && part ~=no_parts(group,subset,trial))
                            order_nb_trial(comb_count)=input_data.options{group}{subset}{trial}.order_neighbors;
                            order_nb_function(comb_count) = order_nb_function_single && ~input_data.options{group}{subset}{trial}.order_neighbors;
                            order_com_function(comb_count) = order_com_function_single && ~(isfield(input_data.options{group}{subset}{trial},'order_com') && input_data.options{group}{subset}{trial}.order_com);

                            
%                             idx1=max((part-2)*floor(no_frames_part_array(group,subset,trial)/2)+1,1);
%                             idx2=min(floor(no_frames_part_array(group,subset,trial)/2)*part+1,no_frames);
                            if parts_and_slices(group,subset,trial) == 1 % slices, no parts. Do not use overlapping bins, use slices as parts.
                                idx1=info.slice_framerange{group,subset,trial}(part,1);
                                idx2=info.slice_framerange{group,subset,trial}(part,2);
                                tr_act = input_data.movementdata{group,subset,trial}{part};
                                if randomized_calcs
                                    tr_actRAND = input_data.rand_movementdata{group,subset,trial}{part};
                                end
                                no_frames_part_array(group,subset,trial) = idx2 - idx1 +1;

                            elseif parts_and_slices(group,subset,trial) == 0 % Only parts, use overlapping bins, no elaborate merging of slices necessary.
                                idx1=max((part-2)*floor(no_frames_part_array(group,subset,trial)/2)+1,1);
                                idx2=min(floor(no_frames_part_array(group,subset,trial)/2)*part+1,no_frames);
                                tr_act = input_data.movementdata{group,subset,trial};
                                if randomized_calcs
                                    tr_actRAND = input_data.rand_movementdata{group,subset,trial};
                                end
                            end
                            
                            
                            
                            if idx2>idx1
                                info_no_focals(comb_count)=info.no_focals(group,subset,trial);
                                info_no_frames(comb_count)=info.no_frames(group,subset,trial);
                                info_framerate(comb_count)=info.framerate(group,subset,trial);
                                info_bodylength(comb_count)=info.bodylength_in_pixels(group,subset,trial);
%                                 if  ~(isnumeric(CustomFilterFrame)&& CustomFilterFrame==-1)
%                                     CustomFilterCell{comb_count} = CustomFilterFrame{group,subset,trial};
%                                 else
%                                     CustomFilterCell{comb_count} = -1;
%                                 end
                                if isfield(info,'circ_roi')
                                    info_circroi(comb_count,:)=squeeze(info.circ_roi(group,subset,trial,:));
                                end
                                info_no_frames_part_array(comb_count)=no_frames_part_array(group,subset,trial);
                                
                                argin_parfor=argin;
                                if any(trajargin)
                                    argin_parfor{trajargin}=tr_act;
                                    
                                end
                                if any(fldinfoargin)
                                    fldstr=strrep(char_input_params,'info.','');
                                    for k=find(fldinfoargin)'
                                        
                                        argin_parfor{k}=squeeze(input_data.info.(fldstr{k})(group,subset,trial,:));
                                    end
                                    %                             argin_parfor(fldinfoargin)=cellfun(@(x) input_data.info.(x)(group,subset,trial),strrep(char_input_params(fldinfoargin),'info.',''),'UniformOutput', false);
                                end
                                
                                argin_scaled{comb_count}=argin_parfor;
                                
                                argin_scaled{comb_count}(cellfun(@(x) isnumeric(x),argin_parfor))=...
                                    cellfun(@(x,y) x.*squeeze(y(group,subset,trial,:,:,:))',argin_parfor(cellfun(@(x) isnumeric(x),argin_parfor)),...
                                    functionInfo.input_params_scaling(cellfun(@(x) isnumeric(x),argin_parfor)),'UniformOutput',false);
                                
%                                 if randomized_calcs
%                                     rand_argin_parfor=argin_parfor;
%                                     rand_argin_parfor{trajargin}=tr_actRAND;
%                                     rand_argin_scaled{comb_count}=rand_argin_parfor;
%                                     rand_argin_scaled{comb_count}(cellfun(@(x) isnumeric(x),rand_argin_parfor))=...
%                                         cellfun(@(x,y) x.*squeeze(y(group,subset,trial,:,:,:))',rand_argin_parfor(cellfun(@(x) isnumeric(x),rand_argin_parfor)),...
%                                         functionInfo.input_params_scaling(cellfun(@(x) isnumeric(x),rand_argin_parfor)),'UniformOutput',false);
%                                     
%                                 end
%                                 
                                
%                                 comb_count=comb_count+1;
                            end
                        end
                    end
                end % Parts for real data
                
                % Parts for random. data
                if randomized_calcs
                for part=1:no_partsRAND(group,subset,trial)
                    [lg,ridx] = ismember([group,subset,trial,part],idx_combsRAND,'rows');
                    if lg
                        comb_count = ridx;
                        if ~skip_end_partsRAND || (part~=1 && part ~=no_partsRAND(group,subset,trial))
                          
                            if parts_and_slicesRAND(group,subset,trial) == 1 % slices, no parts. Do not use overlapping bins, use slices as parts.
                                idx1=info.slice_framerangeRAND{group,subset,trial}(part,1);
                                idx2=info.slice_framerangeRAND{group,subset,trial}(part,2);
                                
                                tr_actRAND = input_data.rand_movementdata{group,subset,trial}{part};
                                no_frames_part_arrayRAND(group,subset,trial) = idx2 - idx1 +1;

                            elseif parts_and_slicesRAND(group,subset,trial) == 0 % Only parts, use overlapping bins, no elaborate merging of slices necessary.
                                idx1=max((part-2)*floor(no_frames_part_arrayRAND(group,subset,trial)/2)+1,1);
                                idx2=min(floor(no_frames_part_arrayRAND(group,subset,trial)/2)*part+1,no_frames);
                                tr_actRAND = input_data.rand_movementdata{group,subset,trial};
                            end
                            
                            
                            
                            if idx2>idx1
                                
                                
                                info_no_frames_part_arrayRAND(comb_count)=no_frames_part_array(group,subset,trial);
                                
                                
                                argin_parfor=argin;
                                if any(fldinfoargin)
                                    fldstr=strrep(char_input_params,'info.','');
                                    for k=find(fldinfoargin)'
                                        
                                        argin_parfor{k}=squeeze(input_data.info.(fldstr{k})(group,subset,trial,:));
                                    end
                                    %                             argin_parfor(fldinfoargin)=cellfun(@(x) input_data.info.(x)(group,subset,trial),strrep(char_input_params(fldinfoargin),'info.',''),'UniformOutput', false);
                                end
                                rand_argin_parfor=argin_parfor;
                                rand_argin_parfor{trajargin}=tr_actRAND;
                                rand_argin_scaled{comb_count}=rand_argin_parfor;
                                rand_argin_scaled{comb_count}(cellfun(@(x) isnumeric(x),rand_argin_parfor))=...
                                    cellfun(@(x,y) x.*squeeze(y(group,subset,trial,:,:,:))',rand_argin_parfor(cellfun(@(x) isnumeric(x),rand_argin_parfor)),...
                                    functionInfo.input_params_scaling(cellfun(@(x) isnumeric(x),rand_argin_parfor)),'UniformOutput',false);
                                
                                
                                
                                
                            end
                        end
                    end
                end % End: Parts for rand. data
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
            idSocial_auxiliaries_message([mfilename ': Could not open parallel pool. Continue without parallel processing.'],progressBox);
            idSocial_auxiliaries_message('Error message is:')
            idSocial_auxiliaries_message(poolerror,progressBox);
        end
        
    end
    
    % Execute function once to estimate output size:
    idx_count=min(2,no_combs);
    
    argin_act=argin_scaled{idx_count};
    if ischar(argin_scaled{idx_count}{trajargin}) % && exist(argin_scaled{idx_count}{trajargin},'file')==2
        
        tr=load(argin_scaled{idx_count}{trajargin});
        opt=tr.tr.options;
        tr = tr.tr;
        if isstruct(tr) && isfield(tr,'Tr') && isfield(tr,'Vel')&&  isfield(tr,'Acc')
            tr.Tr = tr.Tr(opt.start_frame:min(size(tr.Tr,1),opt.end_frame),:,:,:);
            tr.Vel = tr.Vel(opt.start_frame:min(size(tr.Tr,1),opt.end_frame),:,:,:);
            tr.Acc = tr.Acc(opt.start_frame:min(size(tr.Tr,1),opt.end_frame),:,:,:);
%             tr = rmfield(tr,'tr');
        elseif isstruct(tr) && isfield(tr,'Tr')
            tr=tr.Tr(opt.start_frame:min(size(tr.Tr,1),opt.end_frame),:,:,:);
        else
            tr=tr.tr(opt.start_frame:min(size(tr.tr.trajectory,1),opt.end_frame),:,:,:);
        end
       
        
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
        if order_com_function(idx_count)
            % I cannot load the already
            % existing ordered trajectory if I want
            % to apply filters before!
            % Skip if ordering has been applied before.
            tr=idSocial_auxiliaries_nearestCOM(tr);
        end
        argin_act{trajargin}=tr;
    end
    alloutput=cell(no_funcs,1);
    try
        [alloutput{:}]=feval(function_handle,argin_act{:});
        if isempty(funcinfo)
            funcinfoIdx = cellfun(@(x) isstruct(x)&&isfield(x,'Function')&&strcmpi(func2str(function_handle),x.Function),alloutput);
            if any(funcinfoIdx)
                funcinfo = alloutput{funcinfoIdx};
                input_data.(core_act_method).info = funcinfo;
                good_idx(funcinfoIdx) = [];
                no_outfuncs = no_outfuncs - 1;
            end
        end
    catch ME
        warning([mfilename ': Failed to execute ' func2str(function_handle) '. Error message is '])
        idSocial_auxiliaries_message(ME.identifier,progressBox);
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
                limit_frames_temp(idx_count,:)=[input_data.options{group}{subset}{trial}.start_frame, ...
                    input_data.options{group}{subset}{trial}.end_frame];
                
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
        order_com_function_temp = order_com_function(temp_idx);

        idx_conversion(temp_idx,1)=temp_idx;
        
        idx_conversion(temp_idx,2)=1:min(no_combs_per_temp_save,numel(temp_idx));
        
        idx_conversion(temp_idx,3)=temp_save;
        
        worst_indiv_temp = NaN(no_combs_per_temp_save,1);
         all_members_frames_temp = NaN(no_combs_per_temp_save,1);
         all_members_present_for_worst_focal_temp = NaN(no_combs_per_temp_save,1);
        
        if ~use_remote || no_temp_saves == 1
            %par...
            for idx_count=  1:numel(temp_idx)
                try
                idces_act=idx_combs_temp(idx_count,:);
                group=idces_act(1);
                subset=idces_act(2);
                trial=idces_act(3);
                part=idces_act(4);
                no_frames=no_frames_temp(idx_count,:);
                catch
                    keyboard
                end
                
                if parts_and_slices(group,subset,trial) == 1 % slices, no parts. Do not use overlapping bins, use slices as parts.
                    idx1=1;%info.slice_framerange{group,subset,trial}(part,1);
                    idx2=diff(info.slice_framerange{group,subset,trial}(part,:)) + 1;%info.slice_framerange{group,subset,trial}(part,2);
                    no_frames_part_array(group,subset,trial) = idx2 - idx1 +1;
                    idSocial_auxiliaries_message(sprintf('Group %i, subset %i, trial %i, time %.1f-%.1f min',group,subset,trial, ...
                            info.slice_framerange{group,subset,trial}(part,1)/info_framerate_temp(idx_count)/60, ...
                            info.slice_framerange{group,subset,trial}(part,2)/info_framerate_temp(idx_count)/60),progressBox);

                    
                elseif parts_and_slices(group,subset,trial) == 0 % Only parts, use overlapping bins, no elaborate merging of slices necessary.
                    idx1=max((part-2)*floor(info_no_frames_part_array_temp(idx_count,:)/2)+1,1);
                    idx2=min(floor(info_no_frames_part_array_temp(idx_count,:)/2)*part+1,no_frames);
                    idSocial_auxiliaries_message(sprintf('Group %i, subset %i, trial %i, time %.1f-%.1f min',group,subset,trial, ...
                            idx1/info_framerate_temp(idx_count)/60, ...
                            idx2/info_framerate_temp(idx_count)/60),progressBox);
                end
                
                
                if idx2>idx1
%                     fprintf('Group %i, subset %i, trial %i, time %.1f-%.1f min\n',group,subset,trial,idx1/info_framerate_temp(idx_count)/60,idx2/info_framerate_temp(idx_count)/60);
%                     idSocial_auxiliaries_message(sprintf('Group %i, subset %i, trial %i, time %.1f-%.1f min',group,subset,trial,idx1/info_framerate_temp(idx_count)/60,idx2/info_framerate_temp(idx_count)/60),progressBox);
                    %             fprintf('.')
                    argin_act=argin_scaled_temp{idx_count};
                    if ischar(argin_scaled_temp{idx_count}{trajargin}) % && exist(argin_scaled{idx_count}{trajargin},'file')==2
                        
                        tr=load(argin_scaled_temp{idx_count}{trajargin});
                        opt=tr.tr.options;
                        tr = tr.tr;
                        if isstruct(tr) && isfield(tr,'Vel') && isfield(tr,'Acc')
                            tr.Tr = tr.Tr(opt.start_frame:min(size(tr.Tr,1),opt.end_frame),:,:,:);
                            tr.Tr =  tr.Tr(idx1:idx2,:,:,:);
                            try
                            tr.Vel = tr.Vel(opt.start_frame:min(size(tr.Vel,1),opt.end_frame),:,:,:);
                            tr.Vel =  tr.Vel(idx1:idx2,:,:,:);
                            catch
                                keyboard
                            end
                            tr.Acc = tr.Acc(opt.start_frame:min(size(tr.Acc,1),opt.end_frame),:,:,:);
                            tr.Acc =  tr.Acc(idx1:idx2,:,:,:);
%                             tr = rmfield(tr,'tr');
                        elseif isstruct(tr) 
                            tr=tr.Tr(opt.start_frame:min(size(tr.Tr,1),opt.end_frame),:,:,:);
                             try
                             tr=tr(idx1:min(size(tr,1),idx2),:,:,:);
                             catch
                                 keyboard
                             end

                        else
                             tr=tr.Tr(opt.start_frame:min(size(tr.Tr,1),opt.end_frame),:,:,:);
                             try
                             tr=tr(idx1:min(size(tr,1),idx2),:,:,:);
                             catch
                                 keyboard
                             end

                        end
                        %                         limit_frames_temp(idx_count,:)=[opt.start_frame opt.end_frame];
                       
                        
                        
                        
                        
                        CustomFilterFrameArray = [];
                        if ~isempty(CustomFilterFrame)
                            if ~(isnumeric(CustomFilterFrame)&& CustomFilterFrame==-1) % Filter exists
                                if ~isempty(CustomFilterFrame(group,subset,trial,part,:,:)) && ...
                                        ~ all(cellfun(@(x) isempty(x),CustomFilterFrame(group,subset,trial,part,:)))
                                    no_filterFrames = cellfun(@(x) numel(x),squeeze(CustomFilterFrame(group,subset,trial,part,:,:)));
                                    no_filterFrames = max(no_filterFrames(:));
                                    if no_filterFrames> 0
                                        CustomFilterFrameArray = NaN(idx2-idx1+1,no_focals(group,subset,trial),no_neighbors(group,subset,trial));
                                        for ff = 1:no_focals(group,subset,trial)
                                            for nf = 1:no_neighbors(group,subset,trial)
                                                CustomFilterFrameArray(:,ff,nf) = CustomFilterFrame{group,subset,trial,part,ff,nf}(idx1:idx2);
                                            end
                                        end
                                    end
                                else % isempty
                                    CustomFilterFrameArray = false(idx2-idx1+1,no_focals(group,subset,trial),no_neighbors(group,subset,trial));
                                 end
                            else
                                CustomFilterFrameArray = -1;
                            end
                        end
                        
%                         if ~isempty(CustomFilterFrameArray) && isnumeric(CustomFilterFrameArray)
%                             CustomFilter = CustomFilterFrameArray;
%                         else
%                             CustomFilter=-1;
%                         end
                        if applyfilter 
                            
                            
                            
                            tr = idSocial_filters3D(tr,...
                                options,...
                                info_bodylength_temp(idx_count),...
                                info_framerate_temp(idx_count),...
                                info_circroi(idx_count,:), ...
                                CustomFilterFrameArray);
                            %                             worst_indiv_temp(idx_count) = min(min(sum(~isnan(tr(:,:,:,1)),1)./size(tr,1)))*100;
                            %
                            %
                            %
                            %                             trPerm=squeeze(permute(tr(:,:,:,1),[1,2,3,4]));
                            %                             trPerm = trPerm(:,logical(eye(size(tr,2),size(tr,2))));
                            %                             all_members_frames_temp(idx_count) = sum(sum(~isnan(trPerm),2)==size(tr,2))./size(tr,1)*100;
                            %
                            %
                            %                             all_members_present_for_worst_focal_temp(idx_count) = min(sum(squeeze(sum(~isnan(tr(:,:,:,1)),2)==size(tr,2)),1))/no_frames*100;
                            %
                            [worst_indiv_temp(idx_count), all_members_frames_temp(idx_count), all_members_present_for_worst_focal_temp(idx_count)] = ...
                                idSocial_auxiliaries_checkDataQuality(tr);
%                             fprintf('\t%.2f%% of frames left for worst focal-neighbor pair.\n',worst_indiv_temp(idx_count))
%                             fprintf('\t%.2f%% of frames with all individuals present after final filtering.\n',all_members_frames_temp(idx_count))
%                             fprintf('\t%.2f%% of frames with all individuals present for worst focal after final filtering.\n',all_members_present_for_worst_focal_temp(idx_count))
                            
                            idSocial_auxiliaries_message(sprintf('\t%.2f%% of frames left for worst focal-neighbor pair.',worst_indiv_temp(idx_count)),progressBox);
                            idSocial_auxiliaries_message(sprintf('\t%.2f%% of frames with all individuals present after final filtering.',all_members_frames_temp(idx_count)),progressBox);
                            idSocial_auxiliaries_message(sprintf('\t%.2f%% of frames with all individuals present for worst focal after final filtering.',all_members_present_for_worst_focal_temp(idx_count)),progressBox);

                            % Note: all_members_present_for_worst_focal_temp
                            % can be GREATER than all_members_frames_temp,
                           % because all_members_frames_temp considers only
                           % frames where all focals (the diagonals) are
                           % present, which can be filtered e.g. by the
                           % focal velocity filter.
                           % all_members_present_for_worst_focal_temp
                           % however looks if there is one focal with all
                           % its neighbors present, and individuals can 
                           % appear as neighbors even though they do not
                           % appear as focals.
                            
                        end
                        
                        if order_nb_function_temp(idx_count)
                            % I cannot load the already
                            % existing ordered trajectory if I want
                            % to apply filters before!
                            % Skip if ordering has been applied before.
                            tr=idSocial_auxiliaries_nearestneighbours(tr);
                        end
%                         if subset==7 && trial==1 && part==2; keyboard; end
                        if order_com_function_temp(idx_count)
                            % I cannot load the already
                            % existing ordered trajectory if I want
                            % to apply filters before!
                            % Skip if ordering has been applied before.
                            tr=idSocial_auxiliaries_nearestCOM(tr);
                        end
                        argin_act{trajargin}=tr;
                    end
                                                                        
                    alloutput=cell(no_funcs,1);
                    if all_members_frames_temp(idx_count) >= options.filter_AllMembersPresent*100 && ...
                            worst_indiv_temp(idx_count) >= options.filter_WorstIndividual*100	 && ...
                            all_members_present_for_worst_focal_temp(idx_count) >= options.filter_AllMembersPresentWorstFocal*100
                            %                     if worst_indiv>0
                        try
                            
                            if parts_and_slices(group,subset,trial) == 1
                                
                            end
                            
                            [alloutput{:}]=feval(function_handle,argin_act{:});
%                             if isempty(funcinfo)
%                                 funcinfoIdx = cellfun(@(x) isstruct(x)&&isfield(x,'Function')&&strcmpi(func2str(function_handle),x.Function),alloutput);
%                                 funcinfo = alloutput{funcinfoIdx};
%                                 alloutput(funcinfoIdx) = [];
%                                 good_idx(funcinfoIdx) = [];
%                             end
%                             if subset==7 && trial==1 && part==2; keyboard; end
                            
                        catch exception
                            %                             keyboard
                            
                            msgString = getReport(exception);
                            warning([mfilename ': Execution failed. Error message reads:'])
                            idSocial_auxiliaries_message(msgString,progressBox);
                        end
                        if order_nb_function_temp(idx_count) || order_com_function_temp(idx_count)
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
                                        %% Filter out all-nan dims
                                        try
                                        if iscell(alloutput_new{outc}) && ~any(iscell(vertcat(alloutput_new{outc}{ff,nf,:,:,:,:,:}))) && ...
                                                ~isempty(vertcat(alloutput_new{outc}{ff,nf,:,:,:,:,:})) && all(isnan(vertcat(alloutput_new{outc}{ff,nf,:,:,:,:,:})))
                                            alloutput_new{outc}(ff,nf,:,:,:,:,:) = {[]};
                                          
                                        end
                                        catch
                                            warning([mfilename ': Dont be shocked, this is the nan-filter failing, line 1055'])
                                            keyboard
                                        end
                                        
                                    end
                                end
                            end
                            alloutput = alloutput_new;
                            %                         clear alloutput_new
                        end
                         
                        output_parfor(idx_count,:)=alloutput(good_idx);
                    

                    end
                    if all_members_frames_temp(idx_count) < options.filter_AllMembersPresent*100
                         idSocial_auxiliaries_message(['Omitting (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') because of ''filter_AllMembersPresent'''],progressBox);
                    end
                    if worst_indiv_temp(idx_count) < options.filter_WorstIndividual*100
                        idSocial_auxiliaries_message(['Omitting (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') because of ''filter_WorstIndividual'''],progressBox);
                    end
                    if all_members_present_for_worst_focal_temp(idx_count) < options.filter_AllMembersPresentWorstFocal*100
                        idSocial_auxiliaries_message(['Omitting (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') because of ''filter_AllMembersPresentWorstFocal'''],progressBox);
                    end
                    %                     end
                end
                
            end
            if no_temp_saves>1
                
                ofbufferTic=tic;
                for of = 1:no_outfuncs
                    idSocial_auxiliaries_message(['Buffering ' 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat' ' ...'],progressBox);
                    tempOutput =  output_parfor(:,of);
                    savefast([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'tempOutput');
                end
                %                                savefast([temp_savepath_root 'wrapper_buffer_' num2str(temp_save) '.mat'],'output_parfor');
                
                idSocial_auxiliaries_message(['Done (' num2str(toc(ofbufferTic),'%.1f') 's)'],progressBox);
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
            paramstruct.order_com_function_temp = order_com_function_temp;
            paramstruct.no_funcs = no_funcs;
            paramstruct.no_outfuncs = no_outfuncs;
            paramstruct.good_idx = good_idx;
            paramstruct.no_combs_per_temp_save = no_combs_per_temp_save;
            
            act_remote = mod(temp_save-1,no_remotes)+1;
            
            temp_save_delivered = false;
            while ~temp_save_delivered
                if exist([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'],'file')==2
                    idSocial_auxiliaries_message(['Starting round ' num2str(temp_save) ' on Computer ' num2str(act_remote)],progressBox);
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
        worst_indiv = NaN(no_groups,max(no_subsets),no_trials);
        all_members_frames = NaN(no_groups,max(no_subsets),no_trials);
        all_members_present_for_worst_focal = NaN(no_groups,max(no_subsets),no_trials);
        
        temp_save = 0;
        for idx_count=1:no_combs
            
            temp_save_act = idx_conversion(idx_count,3);
            idx_count_temp = idx_conversion(idx_count,2);
            
            if temp_save_act ~= temp_save && no_temp_saves>1
                if exist('output_parfor_act','var')==1
                    clear output_parfor_act
                end
                idSocial_auxiliaries_message(['Loading temporary buffer ' 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat' ' ...'],progressBox);
                %                 keyboard
                load([temp_savepath_root 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat'])
                output_parfor_act = tempOutput;
                clear tempOutput
                idSocial_auxiliaries_message('Done.',progressBox);
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
            worst_indiv(group,subset,trial) = worst_indiv_temp(idx_count_temp); 
            catch
                keyboard
            end
            all_members_frames(group,subset,trial) = all_members_frames_temp(idx_count_temp);
            all_members_present_for_worst_focal(group,subset,trial) = all_members_present_for_worst_focal_temp(idx_count_temp);
            try
                % Note: In the following, it is always
                % output_parfor{idx_count,1} instead of
                % output_parfor{idx_count,outfun}, because
                % the first dimension will be deleted after
                % the loop in order to free memory.
                if ~isempty(output_parfor_act{idx_count_temp})
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
                end
            catch
                keyboard
            end
            
        end
        %         output_parfor(:,1) = [];
        
        % If output is very big, make a hard disk copy and load it when needed in order
        % to save some memory.
        mem_size=whos('output_new');
        user = [];
        try
        user = memory;
        catch
        end
        no_files = ceil(mem_size.bytes/(bytes_per_temp_save));
        if ~isempty(user) && mem_size.bytes/(bytes_per_temp_save)>0;%user.MemAvailableAllArrays>0
            
%             if isfield(input_data,name_outfuncs{outfun}) && ...
%                     isfield(input_data.(name_outfuncs{outfun}),'output') && ...
%                     isfield(input_data.(name_outfuncs{outfun}).output,act_datestr) && ...
%                     ~isempty(input_data.(name_outfuncs{outfun}).output.(act_datestr)) && ...
%                     ischar(input_data.(name_outfuncs{outfun}).output.(act_datestr))
%               
%                 input_data.(name_outfuncs{outfun}).output.(act_datestr)=[];
%             end
            try
                
                
                
                if no_files>1
                    strpath=[temp_savepath_root name_outfuncs{outfun} act_datestr '_.mat'];
                    input_data.(name_outfuncs{outfun}).output.(act_datestr).path=strpath;
                    idSocial_auxiliaries_message(['Buffering output to ' strpath '_1 to _' num2str(no_files) '...'],progressBox);
                    
                    siz_out = size(output_new(:),1);
                    pce_size = ceil(siz_out/no_files);
                    output_old = output_new;
                    clear output_new
                    for file_act = 1:no_files
                        
                        idx_act = (file_act-1)*pce_size+1:min(file_act*pce_size,siz_out);
                        outp.output = output_old(idx_act);
                        outp.size = size(output_old);
                        outp.no_files = no_files;
                        output_new = outp;
                        clear outp 
                        tic;idSocial_auxiliaries_save([strpath(1:end-4) '_' num2str(file_act) '.mat'] ,'output_new'); toc
                        % Check if file was saved
                        if exist([strpath(1:end-4) '_' num2str(file_act) '.mat'],'file')==2
                            idSocial_auxiliaries_message(['Done (' num2str(file_act) ')' ],progressBox);
                        else
                            warning([mfilename ': Buffering failed.'])
                        end
                    end
                    clear output_old 
                else
                    strpath=[temp_savepath_root name_outfuncs{outfun} act_datestr '.mat'];
                    input_data.(name_outfuncs{outfun}).output.(act_datestr).path=strpath;
                    idSocial_auxiliaries_message(['Buffering output to ' strpath ' ...'],progressBox);
                    
                   
                    tic;idSocial_auxiliaries_save(strpath,'output_new'); toc
                    % Check if file was saved
                if exist(strpath,'file')==2
                    idSocial_auxiliaries_message('Done.',progressBox);
                else
                    warning([mfilename ': Buffering failed.'])
                end
                end
                
            catch poolerror
                
                
                delete(strpath)
                warning([mfilename ': Buffering to hard disk failed. Using memory to store variable. ' ...
                    strpath 'deleted.'])
                idSocial_auxiliaries_message('Error message is:',progressBox);
                idSocial_auxiliaries_message(poolerror,progressBox);
                input_data.(name_outfuncs{outfun}).output.(act_datestr).path=output_new;
                keyboard
            end
            
        else
            input_data.(name_outfuncs{outfun}).output.(act_datestr).path=output_new;
            
        end
        input_data.(name_outfuncs{outfun}).output.(act_datestr).options = input_data.(act_method).options;

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
    clear CustomFilterCell CustomFilterFrame CustomFilterFrameArray alloutput output_parfor_act tr
    %%  The same for random data
    
    if randomized_calcs
        
        try
        if use_parallel && ~verLessThan('matlab', '8.3.0')
            p = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(p)
                no_cores_open = 0;
            else
                no_cores_open = p.NumWorkers;
                delete(gcp('nocreate'));
                no_cores_open = 0;
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
        elseif use_parallel
            matlabpool close
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
            idSocial_auxiliaries_message([mfilename ': Could not open parallel pool. Continue without parallel processing.'],progressBox);
            idSocial_auxiliaries_message('Error message is:',progressBox);
            idSocial_auxiliaries_message(poolerror,progressBox);
        end
        %%%%%%%%%%%%%%%%%%%
        % Execute function once to estimate output size:
        idx_count=min(2,no_combs);
        limit_frames_act = limit_frames(idx_count,:);
        rand_argin_act=rand_argin_scaled{idx_count};
        if ischar(rand_argin_scaled{idx_count}{trajargin}) %&& ischar(argin_scaled{idx_count}{trajargin}) % && exist(rand_argin_scaled{idx_count}{trajargin},'file')==2
            rand_tr=load(rand_argin_scaled{idx_count}{trajargin});
%             rand_tr = rand_tr.rand_tr;
            if isfield(rand_tr,'rand_tr') && isstruct(rand_tr.rand_tr)
                rand_tr.Tr = rand_tr.rand_tr.Tr(opt.start_frame:min(size(rand_tr.rand_tr.Tr,1),opt.end_frame),:,:,:);
                rand_tr.Vel = rand_tr.rand_tr.Vel(opt.start_frame:min(size(rand_tr.rand_tr.Tr,1),opt.end_frame),:,:,:);
                rand_tr.Acc = rand_tr.rand_tr.Acc(opt.start_frame:min(size(rand_tr.rand_tr.Tr,1),opt.end_frame),:,:,:);
                rand_tr = rmfield(rand_tr,'rand_tr');
            else
                rand_tr=rand_tr.rand_tr(limit_frames_act(1):min(size(rand_tr.rand_tr,1),limit_frames_act(2)),:,:,:);
            end
            
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
            if order_com_function(idx_count)
                % I cannot load the already
                % existing ordered trajectory if I want
                % to apply filters before!
                rand_tr=idSocial_auxiliaries_nearestCOM(rand_tr);
            end
            
            rand_argin_act{trajargin}=rand_tr;
            clear rand_tr
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
            idx_combs_temp = idx_combsRAND(temp_idx,:);
            no_frames_temp = info_no_frames(temp_idx,:);
            no_focals_temp = info_no_focals(temp_idx,:);
            info_no_frames_part_array_temp = info_no_frames_part_array(temp_idx,:);
            info_framerate_temp = info_framerate(temp_idx);
            rand_argin_scaled_temp = rand_argin_scaled(temp_idx);
            limit_frames_temp = limit_frames(temp_idx,:);
            info_bodylength_temp = info_bodylength(temp_idx);
            order_nb_function_temp = order_nb_function(temp_idx);
            order_com_function_temp = order_com_function(temp_idx);
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
                    
                    
                    if parts_and_slicesRAND(group,subset,trial) == 1 % slices, no parts. Do not use overlapping bins, use slices as parts.
                        idx1=1;%info.slice_framerangeRAND{group,subset,trial}(part,1);
                        idx2=diff(info.slice_framerangeRAND{group,subset,trial}(part,:)) + 1;%info.slice_framerangeRAND{group,subset,trial}(part,2);
                        idSocial_auxiliaries_message(sprintf('Random, Group %i, subset %i, trial %i, time %.1f-%.1f min',group,subset,trial, ...
                            info.slice_framerangeRAND{group,subset,trial}(part,1)/info_framerate_temp(idx_count)/60, ...
                            info.slice_framerangeRAND{group,subset,trial}(part,2)/info_framerate_temp(idx_count)/60),progressBox);

                        
                    elseif parts_and_slicesRAND(group,subset,trial) == 0 % Only parts, use overlapping bins, no elaborate merging of slices necessary.
                        idx1=max((part-2)*floor(info_no_frames_part_array_temp(idx_count,:)/2)+1,1);
                        idx2=min(floor(info_no_frames_part_array_temp(idx_count,:)/2)*part+1,no_frames);
                        idSocial_auxiliaries_message(sprintf('Random, Group %i, subset %i, trial %i, time %.1f-%.1f min',group,subset,trial,idx1/info_framerate_temp(idx_count)/60,idx2/info_framerate_temp(idx_count)/60),progressBox);

                    end
                    if idx2>idx1
%                         fprintf('Random, Group %i, subset %i, trial %i, time %.1f-%.1f min\n',group,subset,trial,idx1/info_framerate_temp(idx_count)/60,idx2/info_framerate_temp(idx_count)/60);

                        
                        rand_argin_act=rand_argin_scaled_temp{idx_count};
                        if ischar(rand_argin_scaled_temp{idx_count}{trajargin}) %&& ischar(argin_scaled{idx_count}{trajargin}) % && exist(rand_argin_scaled{idx_count}{trajargin},'file')==2
                            %                     tr=load(argin_scaled{idx_count}{trajargin});
                            
                            rand_tr=load(rand_argin_scaled_temp{idx_count}{trajargin});
                            try
                                if isstruct(rand_tr.rand_tr) && isnumeric(rand_tr.rand_tr.Tr)
                                    rand_tr.Tr = rand_tr.rand_tr.Tr(opt.start_frame:min(size(rand_tr.rand_tr.Tr,1),opt.end_frame),:,:,:);
                                    rand_tr.Vel = rand_tr.rand_tr.Vel(opt.start_frame:min(size(rand_tr.rand_tr.Tr,1),opt.end_frame),:,:,:);
                                    rand_tr.Acc = rand_tr.rand_tr.Acc(opt.start_frame:min(size(rand_tr.rand_tr.Tr,1),opt.end_frame),:,:,:);
                                    rand_tr = rmfield(rand_tr,'rand_tr');
                                else
                                    rand_tr=rand_tr.rand_tr(limit_frames_act(1):min(size(rand_tr.rand_tr,1),limit_frames_act(2)),:,:,:);
                                end
                                %                                 rand_tr=rand_tr.rand_tr(limit_frames_act(1):min(size(rand_tr.rand_tr,1),limit_frames_act(2)),:,:,:);
                            catch
                                keyboard
                            end
                            rand_tr.Tr=rand_tr.Tr(idx1:min(idx2,size(rand_tr.Tr,1)),:,:,:);
                            rand_tr.Vel=rand_tr.Vel(idx1:min(idx2,size(rand_tr.Vel,1)),:,:,:);
                            rand_tr.Acc=rand_tr.Acc(idx1:min(idx2,size(rand_tr.Acc,1)),:,:,:);
                            % Skip/mark pairs 'random focal - random
                            % neighbor'. In any core function, we can skip
                            % calculation for these pairs.
                            try
                                rand_tr.Tr(:,no_focals_temp+1:end,no_focals_temp+1:end,1) = NaN;
                                rand_tr.Tr(:,no_focals_temp+1:end,no_focals_temp+1:end,2) = Inf;
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
                             
                            if order_com_function_temp(idx_count)
                                % I cannot load the already
                                % existing ordered trajectory if I want
                                % to apply filters before!
                                rand_tr=idSocial_auxiliaries_nearestCOM(rand_tr);
                            end
                            
                            rand_argin_act{trajargin}=rand_tr;
                            
                            rand_alloutput=cell(no_funcs,1);
                            try
                                [rand_alloutput{:}]=feval(function_handle,rand_argin_act{:});
                            catch exception
                                                            
                                
                                msgString = getReport(exception);
                                warning([mfilename ': Execution failed. Error message reads:'])
                                idSocial_auxiliaries_message(msgString,progressBox);
                                keyboard
                            end
                            if order_nb_function_temp(idx_count) || order_com_function_temp(idx_count)
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
                                            %% Filter out all-nan dims
                                            try
                                                if ~any(iscell(alloutput_new{outc}{ff,nf,:,:,:,:,:})) && ...
                                                        ~isempty(alloutput_new{outc}{ff,nf,:,:,:,:,:}) && all(isnan(alloutput_new{outc}{ff,nf,:,:,:,:,:}))
                                                    alloutput_new{outc}(ff,nf,:,:,:,:,:) = {[]};
                                                    
                                                end
                                            catch
                                                warning([mfilename ': Dont be shocked, this is the nan-filter failing, line 1603'])
                                                keyboard
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
                        idSocial_auxiliaries_message(['Buffering ' 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat' ' ...'],progressBox);
                        
                        tempOutput =  rand_output_parfor(:,of);
                        savefast([temp_savepath_root 'tempOutput_of_' num2str(of) '_no_' num2str(temp_save) '.mat'],'tempOutput');
                    end
                    
                    idSocial_auxiliaries_message(['Done (' num2str(toc(ofbufferTic),'%.1f') 's)'],progressBox);
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
                paramstruct.order_com_function_temp = order_com_function_temp;
                paramstruct.no_funcs = no_funcs;
                paramstruct.no_outfuncs = no_outfuncs;
                paramstruct.good_idx = good_idx;
                paramstruct.no_combs_per_temp_save = no_combs_per_temp_save;
                
                act_remote = mod(temp_save-1,no_remotes)+1;
                
                temp_save_delivered = false;
                
                while ~temp_save_delivered
                    if exist([temp_savepath_root 'REMOTE_' num2str(act_remote) '_READY.mat'],'file')==2
                        idSocial_auxiliaries_message(['Starting round ' num2str(temp_save) ' on Computer ' num2str(act_remote)],progressBox);
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
                    if exist('rand_output_parfor_act','var')==1
                        clear rand_output_parfor_act
                    end
                    idSocial_auxiliaries_message(['Loading temporary buffer ' 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat' ' ...'],progressBox);
                    %                 keyboard
                    try
                        load([temp_savepath_root 'tempOutput_of_' num2str(outfun) '_no_' num2str(temp_save_act) '.mat'])
                    catch ME
                        idSocial_auxiliaries_message(ME.identifier,progressBox);
                        keyboard
                    end
                    rand_output_parfor_act = tempOutput;
                    clear tempOutput
                    idSocial_auxiliaries_message('Done.',progressBox);
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
            user = [];
            try
                user = memory;
            catch
            end
            
            if ~isempty(user) && mem_size.bytes/user.MemAvailableAllArrays>0
                if exist(temp_savepath_root,'dir')~=7
                    mkdir(temp_savepath_root)
                end
                
                
                
                % For random:
                
                if isfield(input_data,[name_outfuncs{outfun} 'RAND']) && ...
                        isfield(input_data.([name_outfuncs{outfun} 'RAND']),'output') && ...
                        ~isempty(input_data.([name_outfuncs{outfun} 'RAND'])) && ...
                        ischar(input_data.([name_outfuncs{outfun} 'RAND']).output)
%                     delete(input_data.([name_outfuncs{outfun} 'RAND']).output);
                    input_data.([name_outfuncs{outfun} 'RAND']).output.(act_datestr).path=[];
                end
                
                %%%%
                no_files = ceil(mem_size.bytes/(bytes_per_temp_save));
%                 keyboard
                if no_files>1
                    strpath=[temp_savepath_root name_outfuncs{outfun} 'RAND' act_datestr '_.mat'];
                    input_data.([name_outfuncs{outfun} 'RAND']).output.(act_datestr).path=strpath;
                    idSocial_auxiliaries_message(['Buffering output to ' strpath '_1 to _' num2str(no_files) '...'],progressBox);
                    
                    siz_out = numel(output_new(:));
                    pce_size = ceil(siz_out/no_files);
                    output_old = output_new;
                    clear output_new

                    for file_act = 1:no_files
                        try
                            idx_act = (file_act-1)*pce_size+1:min(file_act*pce_size,siz_out);
                            
                            outp.output = output_old(idx_act);
                            outp.size = size(output_old);
                            outp.no_files = no_files;
                            output_new = outp;
                            clear outp
                            tic;idSocial_auxiliaries_save([strpath(1:end-4) '_' num2str(file_act) '.mat'] ,'output_new'); toc
                            % Check if file was saved
                            if exist([strpath(1:end-4) '_' num2str(file_act) '.mat'],'file')==2
                                idSocial_auxiliaries_message(['Done (' num2str(file_act) ')' ],progressBox);
                            else
                                warning([mfilename ': Buffering failed.'])
                            end
                        catch
                            keyboard
                        end
                    end
                    clear output_old
                else
                    %%%%
                    try
                        strpath=[temp_savepath_root name_outfuncs{outfun} 'RAND' act_datestr '.mat'];
                        input_data.([name_outfuncs{outfun} 'RAND']).output.(act_datestr).path=strpath;
                        idSocial_auxiliaries_message(['Buffering output to ' strpath ' ...'],progressBox);
%                         keyboard
                        idSocial_auxiliaries_save(strpath,'output_new')
                        % Check if file was saved
                        if exist(strpath,'file')==2
                            idSocial_auxiliaries_message('Done.',progressBox);
                        else
                            warning([mfilename ': Buffering failed.'])
                        end
                        clear output_new;
                        
                    catch
                        delete(strpath)
                        warning([mfilename ': Buffering to hard disk failed. Using memory to store variable. ' ...
                            strpath 'deleted.'])
                        input_data.([name_outfuncs{outfun} 'RAND']).output.(act_datestr).path=output_new;
                    end
                end
                
            else
                input_data.([name_outfuncs{outfun} 'RAND']).output.(act_datestr).path=output_new;
                
            end
            input_data.([name_outfuncs{outfun} 'RAND']).output.(act_datestr).options=input_data.(act_method).options;


            
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
        if use_parallel && ~verLessThan('matlab', '8.3.0')
            p = gcp('nocreate');
            if isempty(p)
                no_cores_open = 0;
            else
                no_cores_open = p.NumWorkers;
            end
            if use_parallel  && no_combs>min_number_loop && no_cores_open~=0;
                delete(gcp('nocreate'));
                
            end
        elseif use_parallel
            no_cores_open=matlabpool('size');
            if use_parallel  && no_combs>min_number_loop && no_cores_open~=0;
                matlabpool close;
            end
        end
    catch
    end
    %     input_data.(act_method).info.frames_WorstIndividual = worst_indiv;
    %     input_data.(act_method).info.frames_AllMembersPresent = all_members_frames;
    %     input_data.(act_method).info.frames_AllMembersPresentWorstFocal = all_members_present_for_worst_focal;
    %     input_data.(act_method).info.parts_and_slices = parts_and_slices;
    for outfun=1:no_outfuncs
        input_data.(name_outfuncs{outfun}).info.frames_WorstIndividual = worst_indiv;
        input_data.(name_outfuncs{outfun}).info.frames_AllMembersPresent = all_members_frames;
        input_data.(name_outfuncs{outfun}).info.frames_AllMembersPresentWorstFocal = all_members_present_for_worst_focal;
        input_data.(name_outfuncs{outfun}).info.parts_and_slices = parts_and_slices;
    end
end

% If results for same options has been found in previous set:
% if set_already_calculated>0
% %     input_data.(act_method).(['Set' num2str(act_set)]) = input_data.(act_method).(['Set' num2str(set_already_calculated)]);
%     input_data.(act_method).(['Set' num2str(act_set)]) = input_data.(act_method).output.(outfn{set_already_calculated}).path;
% 
% end

input_data.(principal_method).(['Set' num2str(act_set)]).plot_mode=plot_mode;
input_data.(principal_method).(['Set' num2str(act_set)]).options = input_data.(act_method).options;

if input_data.(act_method).info.parts_and_slices == 1
   plot_mode.parts_and_slices = 1;
end
try
    if actualize_plot_mode || manual_actualize_plot_mode
        
        
        
        if isfield(input_data,[act_method '4Statistics'])
            input_data.(principal_method).(['Set' num2str(act_set)]).output_plot=idSocial_auxiliaries_reduceCellArray(input_data.(act_method).output.(act_datestr).path,...
                plot_mode,input_data.([act_method '4Statistics']).output.(act_datestr).path,[],options);
        else
            input_data.(principal_method).(['Set' num2str(act_set)]).output_plot=idSocial_auxiliaries_reduceCellArray(input_data.(act_method).output.(act_datestr).path,...
                plot_mode,[],[],options);
        end
        
        %     if ~verLessThan('matlab', '14.0.1')
        if ~isempty(funcinfo)
            infofn = fieldnames(funcinfo);
            for k=1:size(infofn,1);
                if ~strcmpi(infofn{k},'Function') && ~strcmpi(infofn{k},'callerFunction')
                    if (strcmpi(infofn{k},'Xtick') || strcmpi(infofn{k},'XtickLabel')) && ~(strcmpi(plot_mode.xaxis{1},'EdgeX') || strcmpi(plot_mode.xaxis{1},'EdgesX'))
                    elseif (strcmpi(infofn{k},'Ytick') || strcmpi(infofn{k},'YtickLabel')) && ~(strcmpi(plot_mode.xaxis{1},'EdgeY') || strcmpi(plot_mode.xaxis{1},'EdgesY'))
                    else
                        input_data.(act_method).(['Set' num2str(act_set)]).output_plot.(infofn{k}) = funcinfo.(infofn{k});
                    end
                    
                    
                end
            end
        end
        
        [~,outTable] = idSocial_auxiliaries_getOutputTable(input_data.(principal_method).(['Set' num2str(act_set)]).output_plot);
        %     end
    end
    
    
    
    if isfield(input_data,[act_method 'RAND']) && (randomized_calcs || isfield(options,'internal_random_controls') && ~isempty(options.internal_random_controls)) 
        if actualize_plot_mode || manual_actualize_plot_mode
            if isfield(input_data,[act_method '4StatisticsRAND'])
                input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND=idSocial_auxiliaries_reduceCellArray(input_data.([act_method 'RAND']).output.(act_datestr).path,...
                    plot_mode,input_data.([act_method '4StatisticsRAND']).output.(act_datestr).path,[],options);
            else
                input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND=idSocial_auxiliaries_reduceCellArray(input_data.([act_method 'RAND']).output.(act_datestr).path,...
                    plot_mode,[],[],options);
            end
            
            if isfield(input_data.(principal_method).(['Set' num2str(act_set)]),'output_plot') && ~isempty(input_data.(principal_method).(['Set' num2str(act_set)]).output_plot)
                pm=input_data.(principal_method).(['Set' num2str(act_set)]).output_plot;
                % Merge plot_mode
                no_entries=size(input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND.data_string,2);
                input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND.data_string=...
                    cellfun(@(x) strrep(x,x,[char(x(1)+no_entries) x(2:end)]),input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND.data_string,'UniformOutput',false);
                input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND.legendstring=...
                    cellfun(@(x) strrep(x,x,[x ', rand.']),input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND.legendstring,'UniformOutput',false);
                
                
                pmrand=input_data.(principal_method).(['Set' num2str(act_set)]).output_plotRAND;
                %             if ~verLessThan('matlab', '14.0.1')
                [~,outTable] = idSocial_auxiliaries_getOutputTable(pm,pmrand);
                %             end
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
                input_data.(principal_method).(['Set' num2str(act_set)]).output_plot=pm;
                input_data.(principal_method).(['Set' num2str(act_set)])=rmfield(input_data.(principal_method).(['Set' num2str(act_set)]),'output_plotRAND');
            end
        end
    end
    % if ~verLessThan('matlab', '14.0.1')
    if actualize_plot_mode || manual_actualize_plot_mode
        if ~isfield(input_data.(principal_method),'results') || isempty(input_data.(principal_method).results)
            input_data.(principal_method).results = outTable;
        else % Update table: merge
            newSiz = size(outTable);
            oldSiz = size(input_data.(principal_method).results);
            newTable = cell(max(newSiz(1),oldSiz(1)),newSiz(2)+oldSiz(2)-1);
            newTable(1:oldSiz(1),1) = input_data.(principal_method).results(:,1);
            newTable(oldSiz(1)+1:newSiz(1),1) = outTable(oldSiz(1)+1:newSiz(1),1);
            
            newTable(1:oldSiz(1),1:oldSiz(2)) = input_data.(principal_method).results;
            newTable(1:newSiz(1),oldSiz(2)+1:oldSiz(2)+newSiz(2)-1) = outTable(:,2:newSiz(2));
            
            input_data.(principal_method).results = newTable;
        end
        % end
        input_data.(principal_method).log_file = logfile;
    end
catch me
    msgText = getReport(me);
    warning([mfilename ': Could not apply statistics (Error in idSocial_auxiliaries_reduceCellArray)'])
    idSocial_auxiliaries_message('Error message is:',progressBox);
    idSocial_auxiliaries_message(msgText,progressBox);
    keyboard
end

if ~isempty(external_func_name)
    input_data.(external_func_name)=input_data.(principal_method);
    if isfield(input_data,[principal_method 'RAND'])
        input_data.([external_func_name 'RAND'])=input_data.([principal_method 'RAND']);
    end
    act_method = external_func_name;
end
% Save input_data
if exist(temp_savepath_root,'dir')~=7
    mkdir(temp_savepath_root)
end

if actualize || manual_actualize || ...
        actualize_plot_mode || manual_actualize_plot_mode
    idSocial_auxiliaries_message(['Buffering input data for ' act_method '...'],progressBox);
    % save(temp_savepath,'input_data','-v7.3');
    idSocial_auxiliaries_save(temp_savepath,'input_data');
    idSocial_auxiliaries_message(['Saved ' temp_savepath],progressBox);
end
diary off

end