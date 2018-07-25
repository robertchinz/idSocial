function [input_data, options, info]=idSocial_loadData(trajectory, datosegmlist, options)
% Load trajectories and additional tracking information and
% apply filters and/or smoothing.
%

act_path=mfilename('fullpath');

% C = strsplit(act_path,filesep);
C = textscan(act_path,'%s','Delimiter',[filesep filesep]);
C=C{1};
idS = strfind(C(1:end-1),'idSocial');
idSocIdx = find(cellfun(@(x) ~isempty(x),idS));
act_path = strcat(C(1:idSocIdx),filesep);
act_path = [act_path{:}];

% idSocIdx=strfind(act_path,'\idSocial');
% act_path=act_path(1:idSocIdx+length('idSocial')+1);
if nargin>0
    try
        cd(act_path)
    catch le
        %     le=lasterror;
        if strcmp(le.identifier,'MATLAB:cd:NonExistentDirectory')
            error([mfilename ': Failed to find idSocial directory.' act_path ' nonexistent or not a directory.'])
        end
    end
    
    addpath(genpath(act_path))
end

%% Set some grafics options

set(0,'DefaultPatchLineSmoothing','On')
set(0,'DefaultLineLineSmoothing','On')
orig_warning_state = warning;
warning('off','MATLAB:audiovideo:VideoReader:fileNotFound')
warning('off','MATLAB:nonExistentField')
input_data=struct();
if nargin<3 || isempty(options)
    options=[];
end

%% Default options
minx = -inf;
maxx = inf;
miny = -inf;
maxy = inf;
maxr = inf;
center = [0,0,0];

act_method=mfilename;
def_options.act_method=act_method(10:end);
def_options.project_name='New Project';%'<path>'; % If this is an empty string, no figures will be saved. If it is path, figures will be saved to path\figures and path\figures_latex.
def_options.framerate=33;
def_options.bodylength_in_pixels=50;
def_options.end_min=inf;
def_options.start_min=0;
def_options.smooth_method= {'none','moving', 'lowess','loess','sgolay','rlowess','rloess',{'moving'}};%{'none','moving', 'lowess','loess','sgolay','rlowess','rloess','adaptive_moving','smoothing_spline','DouglasPeucker','iterative_end_point','split_and_merge','moving_separate'};
def_options.smooth_degree=30;
def_options.random_data=true;
def_options.no_RandomNeighbors = 2;
def_options.filter_focal_list = [];
def_options.filter_neighbors_list = [];
def_options.random_data_method={'shuffle_frames','shuffle_sequences'};
def_options.random_data_seq_length=60;
def_options.order_neighbors = false;
def_options.filter_focal_speedlimits_bl_per_s=[0 inf];
def_options.filter_neighbor_speedlimits_bl_per_s=[0 inf];
def_options.filter_neighbor_accelerationlimits_bl_per_s2 = [0 inf];
def_options.filter_focal_accelerationlimits_bl_per_s2  = [0 inf];
def_options.filter_distancelimits_bl=[0 inf];
def_options.filter_focal_rectangularROI = [minx maxx miny maxy];
def_options.filter_neighbor_rectangularROI = [minx maxx miny maxy];
% def_options.filter_focal_circularROI = [center maxr];
def_options.filter_focal_circularROI = {-Inf Inf {'BL';'Arena'}};

def_options.filter_neighbor_circularROI = {-Inf Inf {'BL';'Arena'}};
def_options.filter_minProbIdentityAssignment = .8;
def_options.filter_AllMembersPresent = 0; % In %!! Not decimals
if nargin < 1
    input_data = def_options;
    return;
end
% def_options.project_path=''; % !!! Repeated and this time it is empty!
def_options.memory_limit_trajectory_MB = 100;
def_options.start_frame=1;
def_options.end_frame=inf;
def_options.focalReconstruction_minProbIdentity4VelCalcs = [];
def_options.focalReconstruction_minProbIdentityFocal = [];
def_options.focalReconstruction_minProbIdentityNeighbor = [];
def_options.significance_between_groups=true;
def_options.smooth_adaptive_noise = 10;
def_options.median_filter=false;
def_options.median_filter_order=1;
def_options.interpolate_trajectories=false;
def_options.interpolation_mode=[];%'linear';
def_options.avg_transform2centerofmass=0;
def_options.smooth_max_deviation=1;
def_options.smooth_spline_degree = 2;

[~ ,def_options]=idSocial_readparams([],[],def_options,act_method);
% def_options.blpxl=def_options.bodylength_in_pixels; % For compability reasons.
def_options.temp_savepath='';




trajdepth = idSocial_recursiveCellSize(trajectory);
datosegmdepth = idSocial_recursiveCellSize(datosegmlist);

if ~isempty(datosegmlist) && trajdepth ~= datosegmdepth
    warning([mfilename ': Size of trajectory list does not match size of datosegm list.'])
end


options = idSocial_recursiveReadParams(options,input_data,def_options,def_options.act_method);

[~, project_name]= idSocial_recursiveGetOptionsFromOptionsCell(options,'project_name');
[~, no_neighborsRAND] = idSocial_recursiveGetOptionsFromOptionsCell(options,'no_RandomNeighbors');

if isempty(project_name)
    project_name=['NewProject_' datestr(now, 30)];
end
if ~isempty(project_name) && strcmpi(project_name,'New Project')
    project_name=['NewProject_' datestr(now, 30)];
end
%% Check input trajectory:
trcell=cell(1,1);
optionscell=cell(1,1);
treat_opt=false;

[~, memory_limit_trajectory_MB]= idSocial_recursiveGetOptionsFromOptionsCell(options,'memory_limit_trajectory_MB');
if isa(trajectory,'char')
    trcell{1}{1}{1}=trajectory;
elseif isa(trajectory,'cell') && size(trajectory,1)==1 && all(cellfun(@(x) isempty(x) || isa(x,'char') ,trajectory))
    % One group, one subgroup and a cell of chars containing trajectory
    % locations.
    trcell{1}{1}=trajectory;
    if numel(options)==numel(trajectory) && (isstruct(options) || all(cellfun(@(x) isstruct(x),options)))
        treat_opt=true;
        if isstruct(options)
            optionscell{1}{1}{1}=options;
        else
            optionscell{1}{1}=options;
        end
    end
elseif isa(trajectory,'cell')
    trajectory = idSocial_recursiveCorrectTrajectoryStructure(trajectory);
    trajectory = idSocial_recursiveCorrectTrajectoryDimensions(trajectory);
    allchar = cellfun(@(x) isa(x,'char'),trajectory);
    if all(allchar)
        trtemp{1} = trajectory;
        trajectory = trtemp;
    end
    allchar=true;
    for gr=1:numel(trajectory)
        if ~isempty(trajectory{gr})
            allchar=allchar & all(cellfun(@(x) isempty(x) || isa(x,'char') ,trajectory{gr}));
        end
    end
    allcell=true;
    for gr=1:numel(trajectory)
        if ~isempty(trajectory{gr})
            allcell=allcell & all(cellfun(@(x) isa(x,'cell') || isempty(x),trajectory{gr}));
        end
    end
    if allchar % No more cells: trajectory{1}='...mat';
        % Insert extra dimension 'group'
        siztr=numel(trajectory);
        for sb=1:siztr
            notr=numel(trajectory{sb});
            for tr=1:notr
                trcell{1}{sb}{tr}=trajectory{sb}{tr};
            end
        end
    elseif allcell
        allchar=inf;
        allchar=true;
        for gr=1:numel(trajectory)
            if ~isempty(trajectory{gr})
                for sb=1:numel(trajectory{gr})
                    if ~isempty(trajectory{gr}{sb})
                        allchar=allchar & all(cellfun(@(x) isempty(x) || isa(x,'char') ,trajectory{gr}{sb}));
                    end
                end
            end
        end
        if  allchar
            % Already in the desired format.
            trcell=trajectory;
        else
            error('loadData:WrongLocationFormat',[mfilename ': Trajectory location does not have the adequate format.'])
        end
        
    end
elseif all(cellfun(@(x) isa(x,'cell') || isempty(x),vertcat(trajectory{:}))) % More cells: trajectory{1}{1}='...mat';
    % Already in the desired format.
    
end

% The following deletes empty groups. This might be useful, but destroys
% the order in case of numbered groups.
% if any(cellfun(@(x) isempty(x),trcell))
%     emptCell=cellfun(@(x) isempty(x),trcell);

%     trcell(emptCell)=[];
% end
trajectory=trcell;
datosegmlist  = idSocial_recursiveDatosegmCellSub2DatosegmCell(trajectory,datosegmlist);

%%
% keyboard
if treat_opt; options = optionscell; end
if isa(options,'cell') && ~treat_opt
    options= idSocial_recursiveOptionsCellSub2OptionsCell(options,trcell);
end
if isa(options,'struct')
    optionscell = trajectory;
    options = idSocial_recursiveOptions2OptionsCell(optionscell,options);
end

clear input_data
options = idSocial_recursiveSetOptionsInOptionsCell(options,'act_method',[]);

[~, temp_savepath]= idSocial_recursiveGetOptionsFromOptionsCell(options,'temp_savepath');
[~, random_data]= idSocial_recursiveGetOptionsFromOptionsCell(options,'random_data');
%% If no path for temporary files is given...
conf_file_exists = false;
if isempty(temp_savepath)
    while conf_file_exists == false
        % Try to read folder information from idSocial.conf, and show menu if
        % it does not exist
        if exist([pwd filesep 'idSocial.conf'],'file')==2
            conf_file_exists = true;
            fileID = fopen([pwd filesep 'idSocial.conf'],'r');
            tempString = textscan(fileID,'%s','Delimiter','');
            fclose(fileID);
            if ~isempty(tempString{1})
                temp_savepath = tempString{1}{1};
            end
        end
        if exist([pwd filesep 'idSocial.conf'],'file')~=2 || isempty(temp_savepath)
            temp_savepath = uigetdir([pwd filesep], ...
                sprintf('Please select directory for temporary files\n(idSocial will save results temporarily to the hard disk\n in order to save memory)'));
            temp_savepath = [temp_savepath filesep];
        end
        % Open idSocial.conf and write folder information if the folder exists. If not
        if  ~isempty(temp_savepath) && exist(temp_savepath,'dir')==7
            fileID = fopen([pwd filesep 'idSocial.conf'],'w');
            fprintf(fileID,'%s',temp_savepath);
            fclose(fileID);
            conf_file_exists =true;
        end
        if exist([pwd filesep 'idSocial.conf'],'file')~=2 || isempty(temp_savepath) || exist(temp_savepath,'dir')~=7
            warning('Could not find folder for temporary files.')
            if exist([pwd filesep 'idSocial.conf'],'file')==2
                delete([pwd filesep 'idSocial.conf'])
            end
            conf_file_exists = false;
            
        end
    end
    %     options.temp_savepath = temp_savepath;
elseif ~isempty(temp_savepath) && ischar(temp_savepath)
    try
        if exist(temp_savepath,'dir')~=7
            mkdir(temp_savepath)
        end
    catch
        error([mfilename ': Could not create directory.'])
    end
end
if ~strcmp(temp_savepath(end),filesep)
    temp_savepath=[temp_savepath filesep];
end
% if ~strcmpi(temp_savepath(1),filesep);
%     temp_savepath = [temp_savepath filesep];
% end
temp_savepath = [temp_savepath project_name];
if ~strcmpi(temp_savepath(end),filesep);
    temp_savepath = [temp_savepath filesep];
end
options = idSocial_recursiveSetOptionsInOptionsCell(options,'temp_savepath',temp_savepath);

%%

no_groups=size(trajectory,2);
no_subsets=max(cellfun(@(x) numel(x),trajectory));
no_trials=max(cellfun(@(y) max(cellfun(@(x) numel(x),y)),trajectory(cellfun(@(x) ~isempty(x),trcell))));

% Set options:--------------------------------------------------------------------------


% Run:----------------------------------------------------------------------------------
% Load trajectories and datosegm:-------------------------------------------------------
% load('trayectorias_20130318_territoriality_rectgroup5dia1trial1.mat')
% load('datosegm_20130318_territoriality_rectgroup5dia1trial1.mat')
trajectorycell=cell(no_groups,no_subsets,no_trials);
probtrajectorycell=cell(no_groups,no_subsets,no_trials);
rand_trajectorycell=cell(no_groups,no_subsets,no_trials);
datosegms=cell(no_groups,no_subsets,no_trials);
movementdata=cell(no_groups,no_subsets,no_trials);
order_movementdata=cell(no_groups,no_subsets,no_trials);
rand_movementdata=cell(no_groups,no_subsets,no_trials);
randOrder_movementdata=cell(no_groups,no_subsets,no_trials);
segmentation_path=cell(no_groups,no_subsets,no_trials);


% info.no_groups=no_groups;
% info.no_trials=no_trials;
info.no_groups=0;
info.no_trials=zeros(no_groups,no_subsets);
info.no_frames=NaN(no_groups,no_subsets,no_trials);
info.duration=NaN(no_groups,no_subsets,no_trials);
info.no_focals=NaN(no_groups,no_subsets,no_trials);
info.no_neighbors=NaN(no_groups,no_subsets,no_trials);
info.no_neighborsRAND=NaN(no_groups,no_subsets,no_trials);
info.framerate=ones(no_groups,no_subsets,no_trials);
info.circ_roi=NaN(no_groups,no_subsets,no_trials,4);
info.bodylength_in_pixels=NaN(no_groups,no_subsets,no_trials);
% info.good_frames=NaN(no_groups,no_subsets,no_trials);
% info.filter_AllMembersPresent=NaN(no_groups,no_subsets,no_trials);
% info.filter_WorstIndividual=NaN(no_groups,no_subsets,no_trials);
% info.filter_AllMembersPresentWorstFocal=NaN(no_groups,no_subsets,no_trials);
% info.filter_AllMembersPresentBeforeFilter=NaN(no_groups,no_subsets,no_trials);
% info.filter_AllMembersPresentWorstFocalBeforeFilter=NaN(no_groups,no_subsets,no_trials);
info.filter_WorstIndividualBeforeFilter=NaN(no_groups,no_subsets,no_trials);
info.roi=NaN(no_groups,no_subsets,no_trials,2,2);
info.bodylength_mm=NaN(no_groups,no_subsets,no_trials);
info.arena_center=NaN(no_groups,no_subsets,no_trials,2);
info.arena_radius=NaN(no_groups,no_subsets,no_trials);
info.group_name=cell(no_groups,no_subsets);
info.no_slices=ones(no_groups,no_subsets,no_trials);
info.slice_framerange = cell(no_groups,no_subsets,no_trials);
info.no_slicesRAND=ones(no_groups,no_subsets,no_trials);
info.slice_framerangeRAND = cell(no_groups,no_subsets,no_trials);

disp('Looking for data...')
for group=1:no_groups
    no_subsets_act=numel(trajectory{group});
    for subset=1:no_subsets_act
        no_trials_act=numel(trajectory{group}{subset});
        for trial=1:no_trials_act
            
            if exist('datosegmlist','var')==1 && ~isempty(datosegmlist) && ~isempty(datosegmlist{group}{subset}) && ~isempty(datosegmlist{group}{subset}{trial}) && ...
                    ischar(datosegmlist{group}{subset}{trial}) && ...
                    exist(datosegmlist{group}{subset}{trial},'file')==2
                
                clear datosegm variable
                warning('off','MATLAB:nonExistentField')
                load(datosegmlist{group}{subset}{trial});
                disp(['Loading ' datosegmlist{group}{subset}{trial}])
                warning('on','MATLAB:nonExistentField')
                if exist('variable','var')==1 && exist('datosegm','var')~=1
                    if isnumeric(variable)
                        datosegms{group,subset,trial}=load_encrypt(datosegmlist{group}{subset}{trial},1);
                    elseif isstruct(variable)
                        datosegms{group,subset,trial}=variable;
                    end
                elseif exist('datosegm','var')==1
                    datosegms{group,subset,trial}=datosegm;
                end;
                
                if isfield(datosegms{group,subset,trial},'roi') && ~isempty(datosegms{group,subset,trial}.roi) && isequal(size(datosegms{group,subset,trial}.roi),[2 2]) % only rectangular roi
                    info.roi(group,subset,trial,:,:)=datosegms{group,subset,trial}.roi;
                elseif isfield(datosegms{group,subset,trial},'tam') && ~isempty(datosegms{group,subset,trial}.tam)
                    info.roi(group,subset,trial,:,:)=[0 0;...
                        datosegms{group,subset,trial}.tam];
                end
                if isfield(datosegms{group,subset,trial},'arena_center') && ~isempty(datosegms{group,subset,trial}.arena_center)
                    info.arena_center(group,subset,trial,:) = datosegms{group,subset,trial}.arena_center;
                end
                if isfield(datosegms{group,subset,trial},'arena_radius') && ~isempty(datosegms{group,subset,trial}.arena_radius)
                    info.arena_radius(group,subset,trial) = datosegms{group,subset,trial}.arena_radius;
                end
                if isfield(datosegms{group,subset,trial},'bodylength') && ~isempty(datosegms{group,subset,trial}.bodylength) && ...
                        (isempty(options{group}{subset}{trial}.bodylength_in_pixels) || ...
                        ischar(options{group}{subset}{trial}.bodylength_in_pixels) && strcmpi(options{group}{subset}{trial}.bodylength_in_pixels,'idTracker'))
                    info.bodylength_in_pixels(group,subset,trial)=nanmean(datosegms{group,subset,trial}.bodylength(:));
                else
                    %                     keyboard
                    disp(['Body length for (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') taken from options, set to ' num2str(options{group}{subset}{trial}.bodylength_in_pixels) ' pxl.'])
                    info.bodylength_in_pixels(group,subset,trial)=options{group}{subset}{trial}.bodylength_in_pixels;
                end
                if isfield(datosegms{group,subset,trial},'bodylength_mm') && ~isempty(datosegms{group,subset,trial}.bodylength_mm) && ...
                        (isempty(options{group}{subset}{trial}.bodylength_mm) || ...
                        ischar(options{group}{subset}{trial}.bodylength_mm) && strcmpi(options{group}{subset}{trial}.bodylength_in_pixels,'idTracker'))
                    info.bodylength_mm(group,subset,trial)=nanmean(datosegms{group,subset,trial}.bodylength_mm(:));
                else
                    %                     keyboard
                    disp(['Body length for (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') taken from options, set to ' num2str(options{group}{subset}{trial}.bodylength_in_pixels) ' pxl.'])
                    info.bodylength_mm(group,subset,trial)=options{group}{subset}{trial}.bodylength_mm;
                end
                if isfield(datosegms{group,subset,trial},'framerate') && ~isempty(datosegms{group,subset,trial}.framerate) && ...
                        isnumeric(datosegms{group,subset,trial}.framerate)
                    info.framerate(group,subset,trial)=nanmean(datosegms{group,subset,trial}.framerate);
                else
                    disp(['Frame rate for (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') taken from options, set to ' num2str(options{group}{subset}{trial}.framerate) ' pxl.'])
                    info.framerate(group,subset,trial)=nanmean(options{group}{subset}{trial}.framerate);
                end
                
                if exist([datosegmlist{group}{subset}{trial}(1:end-12) 'segm_1.mat'],'file')==2
                    segmentation_path{group,subset,trial}=datosegmlist{group}{subset}{trial}(1:end-12);
                end
            elseif (~(exist('datosegmlist','var')==1) || (isempty(datosegmlist) || isempty(datosegmlist{group}{subset}) || isempty(datosegmlist{group}{subset}{trial}))) &&...
                    ~isempty(trajectory{group}{subset}{trial}) && ischar(trajectory{group}{subset}{trial}) &&...
                    (~isempty(strfind(trajectory{group}{subset}{trial},'trayectorias')) || ...
                    ~isempty(strfind(trajectory{group}{subset}{trial},'trajectories')))
                
                datosegm_string=strrep(trajectory{group}{subset}{trial},'trayectorias','datosegm');
                datosegm_string=strrep(datosegm_string,'trajectories','datosegm');
                if exist(datosegm_string,'file')==2
                    clear datosegm variable
                    try
                        load(datosegm_string);
                    catch
                        whos
                        keyboard
                    end
                    if exist('variable','var')==1 && exist('datosegm','var')~=1
                        datosegms{group,subset,trial}=variable;
                    else
                        datosegms{group,subset,trial}=datosegm;
                    end;
                    
                    if isfield(datosegms{group,subset,trial},'psFitParams')
                        datosegms{group,subset,trial} = rmfield(datosegms{group,subset,trial},'psFitParams');
                    end
                    if isfield(datosegms{group,subset,trial},'psFitFunc')
                        datosegms{group,subset,trial} = rmfield(datosegms{group,subset,trial},'psFitFunc');
                    end
                    if isfield(datosegms{group,subset,trial},'ps_list')
                        datosegms{group,subset,trial} = rmfield(datosegms{group,subset,trial},'ps_list');
                    end
                    if isfield(datosegms{group,subset,trial},'day_list')
                        datosegms{group,subset,trial} = rmfield(datosegms{group,subset,trial},'day_list');
                    end
                    
                    if isfield(datosegms{group,subset,trial},'arena_center') && ~isempty(datosegms{group,subset,trial}.arena_center)
                        info.arena_center(group,subset,trial,:) = datosegms{group,subset,trial}.arena_center;
                    end
                    if isfield(datosegms{group,subset,trial},'arena_radius') && ~isempty(datosegms{group,subset,trial}.arena_radius)
                        info.arena_radius(group,subset,trial) = datosegms{group,subset,trial}.arena_radius;
                    end
                    if isfield(datosegms{group,subset,trial},'roi') && ~isempty(datosegms{group,subset,trial}.roi) && ...
                            isequal(size(datosegms{group,subset,trial}.roi),[2 2]) % Rectangular roi
                        info.roi(group,subset,trial,:,:)=datosegms{group,subset,trial}.roi;
                    elseif isfield(datosegms{group,subset,trial},'tam') && ~isempty(datosegms{group,subset,trial}.tam)
                        info.roi(group,subset,trial,:,:)=[0 0;...
                            datosegms{group,subset,trial}.tam];
                    end
                    if isfield(datosegms{group,subset,trial},'bodylength_mm') && ~isempty(datosegms{group,subset,trial}.bodylength_mm)
                        info.bodylength_mm(group,subset,trial)=nanmean(datosegms{group,subset,trial}.bodylength_mm);
                    else
                        info.bodylength_mm(group,subset,trial)=NaN;
                    end
                    if isfield(datosegms{group,subset,trial},'bodylength') && ~isempty(datosegms{group,subset,trial}.bodylength)
                        try
                            info.bodylength_in_pixels(group,subset,trial)=nanmean(datosegms{group,subset,trial}.bodylength(:));
                        catch
                            keyboard
                        end
                        disp([mfilename ': Body length taken from idTracker data file (datosegm.mat) for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(nanmean(datosegms{group,subset,trial}.bodylength)) ' pxl.'])
                        
                    elseif ~ischar(options{group}{subset}{trial}.bodylength_in_pixels)
                        disp([mfilename ': Body length taken from options for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(options{group}{subset}{trial}.bodylength_in_pixels) ' pxl.'])
                        info.bodylength_in_pixels(group,subset,trial)=options{group}{subset}{trial}.bodylength_in_pixels;
                    else
                        disp([mfilename ': Cannot find values for body length. Set to 1 pxl.'])
                        info.bodylength_in_pixels(group,subset,trial)=1;
                        
                    end
                    if isfield(datosegms{group,subset,trial},'framerate') && ~isempty(datosegms{group,subset,trial}.framerate)
                        info.framerate(group,subset,trial)=nanmean(datosegms{group,subset,trial}.framerate);
                        disp([mfilename ': Frame rate taken from idTracker data file (datosegm.mat) for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(nanmean(datosegms{group,subset,trial}.framerate)) ' fps.'])
                        
                    else
                        info.framerate(group,subset,trial)=options{group}{subset}{trial}.framerate;
                    end
                    if isfield(datosegms{group,subset,trial},'arena_center') && ~isempty(datosegms{group,subset,trial}.arena_center) && ...
                            isfield(datosegms{group,subset,trial},'arena_radius') && ~isempty(datosegms{group,subset,trial}.arena_radius)
                        if numel(datosegms{group,subset,trial}.arena_center)==2
                            info.circ_roi(group,subset,trial,:)=[datosegms{group,subset,trial}.arena_center 0 ...
                                datosegms{group,subset,trial}.arena_radius];
                        else
                            info.circ_roi(group,subset,trial,:)=[datosegms{group,subset,trial}.arena_center ...
                                datosegms{group,subset,trial}.arena_radius];
                        end
                    end
                    if exist([datosegm_string(1:end-12) 'segm_1.mat'],'file')==2
                        segmentation_path{group,subset,trial}=datosegm_string(1:end-12);
                    end
                else
                    
                    try
                        info.bodylength_in_pixels(group,subset,trial)=options{group}{subset}{trial}.bodylength_in_pixels;
                        info.framerate(group,subset,trial)=options{group}{subset}{trial}.framerate;
                    catch
                        warning([mfilename ': Could not find file ' datosegm_string])
                        %                         keyboard
                    end
                end
            else
                try
                    info.roi(group,subset,trial,:,:)=NaN(2,2);
                    if isfield(options{group}{subset}{trial},'framerate') && ~isempty(options{group}{subset}{trial}.framerate) && isnumeric(options{group}{subset}{trial}.framerate)
                        info.framerate(group,subset,trial)=options{group}{subset}{trial}.framerate;
                        disp([mfilename ': Frame rate taken from options for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Frame rate = ' num2str(options{group}{subset}{trial}.framerate) ' pxl.'])
                        
                    else
                        info.framerate(group,subset,trial)=NaN;
                    end
                    if ischar(options{group}{subset}{trial}.bodylength_in_pixels)
                        info.bodylength_in_pixels(group,subset,trial)=NaN;
                        
                    else
                        info.bodylength_in_pixels(group,subset,trial)=options{group}{subset}{trial}.bodylength_in_pixels;
                        disp([mfilename ': Body length taken from options for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(options{group}{subset}{trial}.bodylength_in_pixels) ' pxl.'])
                    end
                    
                    if ~isfield(options{group}{subset}{trial},'bodylength_mm') || ischar(options{group}{subset}{trial}.bodylength_mm)
                        info.bodylength_mm(group,subset,trial)=NaN;
                        
                    elseif isfield(options{group}{subset}{trial},'bodylength_mm')
                        info.bodylength_mm(group,subset,trial)=options{group}{subset}{trial}.bodylength_mm;
                        disp([mfilename ': Body length in mm taken from options for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(options{group}{subset}{trial}.bodylength_mm) ' mm.'])
                    end
                    
                    
                catch
                    keyboard
                end
            end
        end
    end
end
disp('Done.')

if exist(temp_savepath,'dir')~=7
    try
        mkdir(temp_savepath)
    catch
        error([mfilename ': Could not create ' temp_savepath '.'])
        
    end
end
if exist([temp_savepath 'trajectories' filesep],'dir')~=7
    try
        mkdir([temp_savepath 'trajectories' filesep])
    catch
        error([mfilename ': Could not create ' [temp_savepath 'trajectories' filesep] '.'])
        
    end
end
%% Try to find out how and where to find the trajectorycell
info.no_groups=no_groups;
no_slices = 1;
no_sliceRAND = 1;
for group=1:no_groups
    if ~isempty(trajectory{group})
        no_subsets_act=numel(trajectory{group});
        info.no_subsets(group)=no_subsets_act;
        for subset=1:no_subsets_act
            if ~isempty(trajectory{group}{subset})
                no_trials_act=numel(trajectory{group}{subset});
                info.no_trials(group,subset)=no_trials_act;
                for trial=1:no_trials_act
                    clear trayectorias trajectories probtrayectorias probtrajectories
                    if ~isempty(trajectory{group}{subset}{trial}) || (~isempty(segmentation_path{group,subset,trial}) && exist([segmentation_path{group,subset,trial} 'trayectorias.mat'],'file')==2)
                        if ischar(trajectory{group}{subset}{trial}) && (exist(trajectory{group}{subset}{trial},'dir')==2 )
                            if ~strcmp(trajectory{group}{subset}{trial}(end),filesep); trajectory{group}{subset}{trial}=[trajectory{group}{subset}{trial} filesep]; end;
                            if exist([trajectory{group}{subset}{trial} 'trajectories.mat'],'file')==2; load([trajectory{group}{subset}{trial} 'trajectories.mat']); end;
                            kaka
                        elseif ischar(trajectory{group}{subset}{trial}) && ...
                                (exist(trajectory{group}{subset}{trial},'file')==2 || exist([trajectory{group}{subset}{trial} '.mat'],'file')==2) ...
                                && exist(trajectory{group}{subset}{trial},'dir')~=2
                            
                            %% Checking contents of file for possible trajectory
                            if strcmpi(trajectory{group}{subset}{trial}(end-3:end),'.mat')
                                act_file=trajectory{group}{subset}{trial};
                            else
                                act_file=[trajectory{group}{subset}{trial} '.mat'];
                            end
                            details = whos('-file',act_file);
                            good_idx = cellfun(@(x) length(x)==3 && x(3)>=2,{details.size}); % Good idx when variable has ndims == 3 (no_frames, no_fish, no_dims) and no_dims >= 2.
                            if sum(good_idx)==1 % Exactly one good idx
                                
                                
                                good_name = details(good_idx).name;
                                load(trajectory{group}{subset}{trial});
                                eval(['trayectorias= ' good_name ';']);
                                disp(['Loading '  trajectory{group}{subset}{trial} '.'])
                            elseif sum(good_idx)>=1
                                warning([mfilename ': File ' act_file ' contains more than one possible trajectory array. The one with a greater number (or the first one, if all no. of frames is same for all) of frames will be chosen automatically.'])
                                no_frames_in_file =  cellfun(@(x) x(1),{details.size});
                                good_idx = good_idx && no_frames_in_file == max(no_frames_in_file);
                                good_name = details(good_idx).name;
                                load(trajectory{group}{subset}{trial});
                                eval(['trayectorias= ' good_name(1) ';']);
                            elseif sum(good_idx)==0
                                error([mfilename ': File ' act_file ' does not contain any possible trajectory array.'])
                            end
                            
                            %                             % Check variable size
                            %                             varsize_MB = prod(details.size)*details.size(2)*8* 1e-6; % "Future" size, when neighbor dimension will be added.
                            %                             if varsize_MB > memory_limit_trajectory_MB
                            %                                 no_slices = ceil(varsize_MB/memory_limit_trajectory_MB);
                            %                                 info.slices{group,subset} = no_slices;
                            %                                 disp([mfilename ': Trajectory is big (' num2str(varsize_MB) 'MB after adding neighbor dimension).'])
                            %                                 disp(['We will slice them into ' num2str(no_slices)  ' pieces of ' num2str(memory_limit_trajectory_MB) 'MB each.'])
                            %                             end
                            
                            good_probtr_idx = cellfun(@(x) length(x)==2 && size(trayectorias,1)==x(1) && ...
                                size(trayectorias,2)==x(2),{details.size});
                            if sum(good_probtr_idx)==1 % Exactly one good idx
                                good_name = details(good_probtr_idx).name;
                                eval(['probtrayectorias= ' good_name ';']);
                            elseif sum(good_probtr_idx)>=1
                                eval(['probtrayectorias= ' good_probtr_name(1) ';']);
                            elseif sum(good_idx)==0
                                disp([mfilename ': File ' act_file ' does not contain any possible array for assignation probabilities.'])
                            end
                        elseif isempty(trajectory{group}{subset}{trial}) && (~isempty(segmentation_path{group,trial}) && exist([segmentation_path{group,trial} 'trajectories.mat'],'file')==2)
                            %                             details = whos('-file',[segmentation_path{group,trial} 'trajectories.mat']);
                            %                             % Check variable size
                            %                             varsize_MB = prod(details.size)*details.size(2)*8* 1e-6; % "Future" size, when neighbor dimension will be added.
                            %                             if varsize_MB > memory_limit_trajectory_MB
                            %                                 no_slices = ceil(varsize_MB/memory_limit_trajectory_MB);
                            %                                 info.slices{group,subset} = no_slices;
                            %                                 disp([mfilename ': Trajectory is big (' num2str(varsize_MB) 'MB after adding neighbor dimension).'])
                            %                                 disp(['We will slice them into ' num2str(no_slices)  ' pieces of ' num2str(memory_limit_trajectory_MB) 'MB each.'])
                            %                             end
                            load([segmentation_path{group,trial} 'trajectories.mat']);
                            disp(['Loading '  segmentation_path{group,trial} 'trajectories.mat'.'])
                            
                        elseif ischar(trajectory{group}{subset}{trial}) && ...
                                exist(trajectory{group}{subset}{trial},'file')~=2
                            error([mfilename ': File ' trajectory{group}{subset}{trial} ' not found!'])
                        else
                            trayectorias=[];%trajectory{group}{subset}{trial};
                            trajectory{group}{subset}{trial}=[];
                        end
                    end
                    if ~isempty(trajectory{group}{subset}{trial})
                        if ~isempty(options{group}{subset}{trial}.start_frame)
                            start_frame=max(options{group}{subset}{trial}.start_frame,1);
                        else
                            start_frame=1;
                        end
                        if ~isempty(options{group}{subset}{trial}.end_frame)
                            end_frame=min(options{group}{subset}{trial}.end_frame,size(trayectorias,1));
                        else
                            end_frame=size(trayectorias,1);
                        end
                        if ~isempty(options{group}{subset}{trial}.start_min) && start_frame == 1
                            start_frame=(options{group}{subset}{trial}.start_min*60*info.framerate(group,subset,trial))+1;
                        else
                            start_frame=1;
                        end
                        if ~isempty(options{group}{subset}{trial}.end_min) && (end_frame == inf || end_frame == size(trayectorias,1))
                            end_frame=min(options{group}{subset}{trial}.end_min*60*info.framerate(group,subset,trial),size(trayectorias,1));
                        else
                            end_frame=size(trayectorias,1);
                        end
                        
                        trPreOrig = trayectorias(start_frame:end_frame,:,:);
                        
                        % Treat probtrajectories
                        if ~exist('probtrayectorias','var')==1 || isempty(probtrayectorias)
                            probtrayectorias = ones(end_frame-start_frame+1,size(trayectorias,2));
                        else
                            probtrayectorias = probtrayectorias(start_frame:end_frame,:);
                        end
                        
                        disp(['Preprocessing  (' num2str(group) ',' num2str(subset) ',' num2str(trial) ')'])
                        
                        try
                            probIdCheck = isfinite(options{group}{subset}{trial}.filter_minProbIdentityAssignment) && any(all(isnan(probtrayectorias) | ...
                                probtrayectorias < options{group}{subset}{trial}.filter_minProbIdentityAssignment,1));
                            
                            tr=idSocial_prepareTrajectories3D(trPreOrig,...
                                probtrayectorias,...
                                options{group}{subset}{trial},info.bodylength_in_pixels(group,subset,trial));
%                             trPreOrig=tr.trajectory;
                            if probIdCheck
                                warning([mfilename ': Assignment probability of at least one individual below threshold (filter_minProbIdentityAssignment) during whole trial.'])
                            end
                        catch
                            keyboard
                        end
                        
                        save([temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) 'Orig.mat'],'tr')
                        trPre = tr;
                        clear tr;
                        % END: Treat probtrajectories
                        
                        no_dim = size(trPre.Tr,3);
                        % Check variable size
                        varsize_MB = size(trPre.Tr,1)*size(trPre.Tr,2)*size(trPre.Tr,2)*no_dim*8 * 1e-6; % "Future" size, when neighbor dimension will be added.
                        if varsize_MB > memory_limit_trajectory_MB
                            no_slices = ceil(varsize_MB/memory_limit_trajectory_MB);
                            info.no_slices(group,subset,trial) = no_slices;
                            disp([mfilename ': Trajectory is big (' num2str(varsize_MB) 'MB after adding neighbor dimension).'])
                            disp(['We will slice them into ' num2str(no_slices)  ' pieces of ' num2str(memory_limit_trajectory_MB) 'MB each.'])
                        end
                        
                        strpathCell = cell(info.no_slices(group,subset,trial),1);
                        if info.no_slices(group,subset,trial)>1
                            strpathCell = cell(info.no_slices(group,subset,trial),1);
                            no_frames_slice = ceil((end_frame-start_frame+1)/info.no_slices(group,subset,trial));
                            info.slice_framerange{group,subset,trial} = NaN(info.no_slices(group,subset,trial),2);
                            for trslice = 1:info.no_slices(group,subset,trial)
                                %                                 if trslice<no_slices
                                %                                     fprintf('%d,',trslice)
                                %
                                %                                 elseif trslice==no_slices
                                %                                     fprintf('%d\n',trslice)
                                %
                                %                                 end
%                                 info.slice_framerange{group,subset,trial}(trslice,:) =  ...
%                                     [(trslice-1)*no_frames_slice+1 min(trslice*no_frames_slice+3,end_frame)];
%                                 frRange = (trslice-1)*no_frames_slice+1 : min(trslice*(no_frames_slice+2),end_frame) ; % Overlap of 2 so vel and acc can be calculated "seamingless"(?)
                                % No overlap necessary with structs,
                                % because Vel and Acc are already there:
                                frRange = (trslice-1)*no_frames_slice+1 : min(trslice*(no_frames_slice),end_frame) ;
                                info.slice_framerange{group,subset,trial}(trslice,:) =  ...
                                    [frRange(1) frRange(end)];
                                tr.Tr  = trPre.Tr(frRange,:,:); 
                                tr.options = trPre.options;
                                tr.Vel = trPre.Vel(frRange,:,:);
                                tr.Acc = trPre.Acc(frRange,:,:);

                                strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trslice) '.mat'];
                                savefast(strpath,'tr');
                                trajectorycell{group,subset,trial}=strpath;
                                strpathCell{trslice} = strpath;
                            end
                            trajectorycell{group,subset,trial} = strpathCell;
                            
                        else
                            tr  = trPre;
                            tr.options = trPre.options;
                            strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                            strpathCell{1} = strpath;
                            savefast(strpath,'tr');
                            trajectorycell{group,subset,trial}=strpath;
                        end
                        
                        if isstruct(trPre)
                            info.no_focals(group,subset,trial)=size(trPre.Tr,2);
                            info.no_neighbors(group,subset,trial)=size(trPre.Tr,2);
                        else
                            info.no_focals(group,subset,trial)=size(trPre,2);
                            info.no_neighbors(group,subset,trial)=size(trPre,2);
                        end
                        clear tr trPre
                        
                        
                        info.no_frames(group,subset,trial)=end_frame-start_frame+1;
                        
                        if info.no_slices(group,subset,trial)>1
                            strpathCellProb = cell(info.no_slices(group,subset,trial),1);
                            for trslice = 1:info.no_slices(group,subset,trial)
                                no_frames_slice = ceil((end_frame-start_frame+1)/info.no_slices(group,subset,trial));
%                                 frRange = (trslice-1)*no_frames_slice+1 : min(trslice*(no_frames_slice+2),end_frame) ; % Overlap of 2 so vel and acc can be calculated "seamingless"(?)
                                % No overlap necessary with structs,
                                % because Vel and Acc are already there:
                                frRange = (trslice-1)*no_frames_slice+1 : min(trslice*(no_frames_slice),end_frame) ;
                                trProb =probtrayectorias(frRange ,:);
                                strpathProb=[temp_savepath 'trajectories' filesep 'trProb_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trslice) '.mat'];
                                savefast(strpathProb,'trProb');
                                strpathCellProb{trslice} = strpathProb;
                                savefast(strpathProb,'trProb');
                            end
                            probtrajectorycell{group,subset,trial}=strpathCellProb;
                            
                        else
                            strpathProb=[temp_savepath 'trajectories' filesep 'trProb_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                            trProb = probtrayectorias;
                            savefast(strpathProb,'trProb');
                            probtrajectorycell{group,subset,trial}=strpathProb;
                            info.slice_framerange{group,subset,trial} = [1 info.no_frames(group,subset,trial)];
                            
                        end
                        clear trProb probtrayectorias
                        
                        %                         info.no_focals(group,subset,trial)=size(trayectorias(start_frame:end_frame,:,:),2);
                        
                        info.duration(group,subset,trial)=info.no_frames(group,subset,trial)./info.framerate(group,subset,trial)/60;
                        %             if trial==1; info.no_groups=info.no_groups+1; end
                        %             info.no_trials(group,subsets)=info.no_trials(group)+1;
                        if all(isnan(info.roi(group,subset,trial,:,:)))
                            minx=floor(min(min(trayectorias(start_frame:end_frame,:,1))));
                            miny=floor(min(min(trayectorias(start_frame:end_frame,:,2))));
                            maxx=round(floor(max(max(trayectorias(start_frame:end_frame,:,1)))));
                            maxy=round(floor(max(max(trayectorias(start_frame:end_frame,:,2)))));
                            info.roi(group,subset,trial,:,:) = [miny minx; maxy maxx];
                        end
                        if all(isnan(info.circ_roi(group,subset,trial,:))) && ...
                                ~all(isnan(info.roi(group,subset,trial,:)))
                            info.circ_roi(group,subset,trial,:) = [info.roi(group,subset,trial,1,2)+(info.roi(group,subset,trial,2,2)-info.roi(group,subset,trial,1,2))/2 ...
                                info.roi(group,subset,trial,1,1)+(info.roi(group,subset,trial,2,1)-info.roi(group,subset,trial,1,1))/2 ...
                                0 ...
                                min((info.roi(group,subset,trial,2,2)-info.roi(group,subset,trial,1,2))/2,(info.roi(group,subset,trial,2,1)-info.roi(group,subset,trial,1,1))/2)];
                            if all(isnan(info.arena_center(group,subset,trial,:)))
                                info.arena_center(group,subset,trial,:) = info.circ_roi(group,subset,trial,1:2);
                            end
                            if all(isnan(info.arena_radius(group,subset,trial,:)))
                                info.arena_radius(group,subset,trial,:) = info.circ_roi(group,subset,trial,4);
                            end
                            
                            disp([mfilename ': Circular roi calculated from extreme coordinates for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': [center X, Y, Z, radius] = [' num2str(info.circ_roi(group,subset,trial,1)) ',' num2str(info.circ_roi(group,subset,trial,2)) ',' num2str(info.circ_roi(group,subset,trial,3)) ',' num2str(info.circ_roi(group,subset,trial,4)) ']'])
                            
                        end
                        
                    end
                end
                if nargin>2 && ~isempty(options{group}{subset}{trial}) && isfield(options{group}{subset}{trial},'group_names') && ~isempty(options{group}{subset}{trial}.group_names) && ...
                        size(options{group}{subset}{trial}.group_names,2)>=no_groups
                    info.group_name{group}=options{group}{subset}{trial}.group_names{group};
                else
                    info.group_name{group}=['Gr. ' num2str(group)];
                end
                %                 keyboard
                %                 info.trials{group,subset}=find(cellfun(@(x) ~isempty(x),trajectory{group}(subset)));
                
                
                try
                    if ~isempty(trajectory{group}{subset})
                        tr_idces=find(cellfun(@(x) ~isempty(x),trajectory{group}{subset}));
                        if size(tr_idces,2)>=1
                            info.trials{group,subset}=tr_idces';
                        elseif size(tr_idces,1)>1
                            info.trials{group,subset}=tr_idces;
                        else
                            info.trials{group,subset} = [];
                        end
                        
                    else
                        info.trials{group,subset}=[];
                    end
                catch
                    keyboard
                end
                
            end % ~isempty(trajectory{group}{subset))
        end
    end % if ~isempty(trajectory{group})
end


%% Filtering
info.good_frames=NaN(no_groups,no_subsets,no_trials,no_slices);
info.filter_AllMembersPresent=NaN(no_groups,no_subsets,no_trials,no_slices);
info.filter_WorstIndividual=NaN(no_groups,no_subsets,no_trials,no_slices);
info.filter_AllMembersPresentWorstFocal=NaN(no_groups,no_subsets,no_trials,no_slices);
info.filter_AllMembersPresentBeforeFilter=NaN(no_groups,no_subsets,no_trials,no_slices);
info.filter_AllMembersPresentWorstFocalBeforeFilter=NaN(no_groups,no_subsets,no_trials,no_slices);
filter_on = true;
if filter_on
    for group=1:no_groups
        if ~isempty(trajectory{group})
            no_subsets_act=numel(trajectory{group});
            for subset=1:no_subsets_act
                if ~isempty(trajectory{group}{subset})
                    no_trials_act=numel(trajectory{group}{subset});
                    for trial=1:no_trials_act
                        if ~isempty(trajectory{group}{subset}{trial})
                            
                            if no_slices>1
                                strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                                movementdata{group,subset,trial}=strpath;
                            else
                                strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                                movementdata{group,subset,trial}=strpath;
                            end
                            order_strpathCell = cell(info.no_slices(group,subset,trial),1);
                            
                            for trslice = 1:no_slices
                                
                                if no_slices>1
                                    
                                    strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trslice) '.mat'];
                                    order_strpath=[temp_savepath 'trajectories' filesep 'trORDER_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trslice) '.mat'];
                                    strpathProb=[temp_savepath 'trajectories' filesep 'trProb_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trslice) '.mat'];
%                                     movementdata{group,subset,trial}=strpathCell;

                                else
                                    strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                                    order_strpath=[temp_savepath 'trajectories' filesep 'trORDER_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                                    strpathProb=[temp_savepath 'trajectories' filesep 'trProb_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
%                                     movementdata{group,subset,trial}=strpath;
                                end
                                order_strpathCell{trslice} = order_strpath;
                                strpathCell{trslice} = strpath;
                                if ~isempty(trajectorycell{group,subset,trial})
                                    
                                    
                                    
                                    load(strpathProb)
                                    tr = load(strpath);
                                    tr = tr.tr; % Get back to struct tr.Tr, tr.options
                                    
                                    
                                    % Apply filters:
                                    if ~isempty(info.bodylength_in_pixels(group,subset,trial)) && ~isnan(info.bodylength_in_pixels(group,subset,trial))
                                        bodylength = info.bodylength_in_pixels(group,subset,trial);
                                    elseif ~isempty(options{group}{subset}{trial}.bodylength_in_pixels) && ~isnan(options{group}{subset}{trial}.bodylength_in_pixels)
                                        bodylength = options{group}{subset}{trial}.bodylength_in_pixels;
                                    else
                                        bodylength = [];
                                    end
                                    
                                    if ~isempty(info.framerate(group,subset,trial)) && ~isnan(info.framerate(group,subset,trial))
                                        framerate = info.framerate(group,subset,trial);
                                    elseif ~isempty(options{group}{subset}{trial}.framerate) && ~isnan(options{group}{subset}{trial}.framerate)
                                        framerate = options{group}{subset}{trial}.framerate;
                                    else
                                        framerate = [];
                                    end
                                    
                                    try
                                        [good_frames_before_filter, all_members_frames_before_filter, all_members_present_for_worst_focal_before_filter] = idSocial_auxiliaries_checkDataQuality(tr.Tr);
                                    catch
                                        keyboard
                                    end
                                    opts = tr.options;
                                    tr  = idSocial_filters3D(tr,options{group}{subset}{trial},bodylength,framerate,squeeze(info.circ_roi(group,subset,trial,:)));
                                    tr.options = opts;
                                    [good_frames, all_members_frames,all_members_present_for_worst_focal ] = idSocial_auxiliaries_checkDataQuality(tr.Tr);
                                    savefast(strpath,'tr')
                                    if no_slices>1
                                        if trslice<no_slices
                                            fprintf('%d,',trslice)
                                            
                                        elseif trslice==no_slices
                                            fprintf('%d\n',trslice)
                                            
                                        end
                                    end
                                    
                                    info.filter_AllMembersPresentBeforeFilter(group,subset,trial,trslice) = all_members_frames_before_filter;
                                    info.filter_WorstIndividualBeforeFilter(group,subset,trial,trslice) = good_frames_before_filter;
                                    info.filter_AllMembersPresentWorstFocalBeforeFilter(group,subset,trial,trslice) = all_members_present_for_worst_focal_before_filter;
                                    info.filter_AllMembersPresent(group,subset,trial,trslice) = all_members_frames;
                                    info.filter_WorstIndividual(group,subset,trial,trslice) = good_frames;
                                    info.filter_AllMembersPresentWorstFocal(group,subset,trial,trslice) = all_members_present_for_worst_focal;
                                    info.good_frames(group,subset,trial)=good_frames;
                                    
                                    % End: Apply filters
                                    
                                    
                                    
                                else
                                    empty_options=options{group}{subset}{trial};
                                    empty_options.interpolate_trajectories=0;
                                    movementdata{group,subset,trial}=strpath;
                                    order_movementdata{group,subset,trial}=order_strpath;
                                    tr=...
                                        idSocial_prepareTrajectories3D(NaN(end_frame-start_frame+1,2,2),...
                                        NaN(end_frame-start_frame+1,2),empty_options);
                                    trOrder=tr;
                                    savefast(strpath,'tr')
                                    savefast(strpath,'trOrder')
                                    if random_data && info.no_focals(group,subset,trial)>2
                                        rand_movementdata{group,subset,trial}=rand_strpath;
                                        randOrder_movementdata{group,subset,trial}=randOrder_strpath;
                                        clear rand_tr;
                                        rand_tr=...
                                            idSocial_prepareTrajectories3D(NaN(end_frame-start_frame+1,2,2),...
                                            NaN(end_frame-start_frame+1,2),...
                                            empty_options);
                                        savefast(rand_strpath,'rand_tr')
                                        randtrOrder=rand_tr;
                                        savefast(randOrder_strpath,'randtrOrder')
                                        clear rand_tr;
                                    end
                                    clear tr;
                                    info.no_frames(group,subset,trial)=0;
                                    info.bodylength_in_pixels(group,subset,trial)=-1;
                                    info.no_focals(group,subset,trial)=-1;
                                    info.no_neighbors(group,subset,trial)=-1;
                                    info.no_neighborsRAND(group,subset,trial)=-1;
                                    info.duration(group,subset,trial)=0;
                                    info.no_slice(group,subset,trial)=0;
                                end
                                
                            end
                            if no_slices>1
                                movementdata{group,subset,trial}=strpathCell;
                            else
                                movementdata{group,subset,trial}=strpath;
                            end
                            
                        end
                    end
                    
                    
                    if nanmean(info.filter_AllMembersPresent(group,subset,trial,:))<options{group}{subset}{trial}.filter_AllMembersPresent*100
                        warning(['Omitting (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') because of ''filter_AllMembersPresent'''])
                        movementdata{group,subset,trial}=[];
                        order_movementdata{group,subset,trial}=[];
                        rand_movementdata{group,subset,trial}=[];
                        randOrder_movementdata{group,subset,trial}=[];
                        info.no_frames(group,subset,trial)=0;
                        info.bodylength_in_pixels(group,subset,trial)=-1;
                        info.no_focals(group,subset,trial)=-1;
                        info.no_neighbors(group,subset,trial)=-1;
                        info.no_neighborsRAND(group,subset,trial)=-1;
                        info.duration(group,subset,trial)=0;
                    else
                        if good_frames<70
                            warning(['% of remaining frames after filtering (Min. for any individual): ' num2str(nanmean(info.filter_WorstIndividual(group,subset,trial,:)))])
                        else
                            disp(['% of remaining frames after filtering (Min. for any individual): ' num2str(nanmean(info.filter_WorstIndividual(group,subset,trial,:)))])
                        end
                        if all_members_frames<70
                            warning(['% of frames in which all individuals are present: ' num2str(nanmean(info.filter_AllMembersPresent(group,subset,trial,:)),'%.1f')])
                        else
                            disp(['% of frames in which all individuals are present: ' num2str(nanmean(info.filter_AllMembersPresent(group,subset,trial,:)),'%.1f')])
                        end
                        for trslice = 1:no_slices
%                             if no_slices>1
                                strpath = strpathCell{trslice};
%                             end
                            %                         rand_strpath = rand_strpathCell{trslice};
                            order_strpath = order_strpathCell{trslice};
                            %                         randOrder_strpath = randOrder_strpathCell{trslice};
                            
                            if no_slices > 1
                                order_movementdata{group,subset,trial}=order_strpathCell;
                            else
                                order_movementdata{group,subset,trial}=order_strpath;
                            end
                            
                            load(strpath)
                            
                            trOrder = idSocial_auxiliaries_nearestneighbours(tr);
                            
                            
                            
%                             if options{group}{subset}{trial}.order_neighbors % Substitute original trajectories by ordered ones
%                                 tr = trOrder;
%                             end
                            savefast(order_strpath,'trOrder')

                            clear trOrder;
                            
%                             savefast(strpath,'tr')
                            clear tr;
                            
                        end
                    end
                    
                end %  ~isempty(trajectory{group}{subset})
            end
        end %  ~isempty(trajectory{group})
    end
end
%% Generate randomized trajectories


if random_data
    for group=1:no_groups
        if ~isempty(trajectory{group})
            no_subsets_act=numel(trajectory{group});
            for subset=1:no_subsets_act
                if ~isempty(trajectory{group}{subset})
                    no_trials_act=numel(trajectory{group}{subset});
                    
                    for trial=1:no_trials_act
                        no_neighborsRAND = max(no_neighborsRAND,info.no_focals(group, subset, trial));
                        options{group}{subset}{trial}.no_RandomNeighbors = no_neighborsRAND;
                        if info.no_focals(group,subset,trial)>1
                            
                            if ~isempty(options{group}{subset}{trial}.start_frame)
                                start_frame=max(options{group}{subset}{trial}.start_frame,1);
                            else
                                start_frame=1;
                            end
                            if ~isempty(options{group}{subset}{trial}.end_frame)
                                end_frame=min(options{group}{subset}{trial}.end_frame,info.no_frames(group,subset,trial));
                            else
                                end_frame=size(trayectorias,1);
                            end
                            if ~isempty(options{group}{subset}{trial}.start_min) && start_frame == 1
                                start_frame=(options{group}{subset}{trial}.start_min*60*info.framerate(group,subset,trial))+1;
                            else
                                start_frame=1;
                            end
                            if ~isempty(options{group}{subset}{trial}.end_min) && (end_frame == inf || end_frame == info.no_frames(group,subset,trial))
                                end_frame=min(options{group}{subset}{trial}.end_min*60*info.framerate(group,subset,trial),info.no_frames(group,subset,trial));
                            else
                                end_frame=info.no_frames(group,subset,trial);
                            end
                            
                            
                            % Check variable size
                            varsize_MB = max(info.no_frames(:))*no_neighborsRAND*no_neighborsRAND*8*no_dim* 1e-6; % "Future" size, when neighbor dimension will be added.
                            if varsize_MB > memory_limit_trajectory_MB
                                info.no_slicesRAND(group,subset,trial) = ceil(varsize_MB/memory_limit_trajectory_MB);
%                                 info.no_slicesRAND(group,subset,trial) = no_slicesRAND;

                                disp([mfilename ': Rand. trajectory will be big (' num2str(varsize_MB) 'MB after adding neighbor dimension).'])
                                disp(['We will slice them into ' num2str(info.no_slicesRAND(group,subset,trial))  ' pieces of ' num2str(memory_limit_trajectory_MB) 'MB each.'])
                            end
                            rand_strpathCell = cell(info.no_slicesRAND(group,subset,trial),1);
                            randOrder_strpathCell = cell(info.no_slicesRAND(group,subset,trial),1);
                            
                            % Load 'trPreOrig' (no_frames,no_fish,no_dims)
                            load([temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) 'Orig.mat'])
                            disp(['Preprocessing randomization (' num2str(group) ',' num2str(subset) ',' num2str(trial) ')'])
                            info.slice_framerangeRAND{group,subset,trial} = NaN(info.no_slicesRAND(group,subset,trial),2);
                            for trsliceRAND = 1:info.no_slicesRAND(group,subset,trial)
                                if info.no_slicesRAND(group,subset,trial)>1
                                    
                                    rand_strpath=[temp_savepath 'trajectories' filesep 'trRAND_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trsliceRAND) '.mat'];
                                    randOrder_strpath=[temp_savepath 'trajectories' filesep 'trRANDORDER_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trsliceRAND) '.mat'];
                                    
                                else
                                    rand_strpath=[temp_savepath 'trajectories' filesep 'trRAND_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                                    randOrder_strpath=[temp_savepath 'trajectories' filesep 'trRANDORDER_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                                end
                                rand_strpathCell{trsliceRAND} = rand_strpath;
                                randOrder_strpathCell{trsliceRAND} = randOrder_strpath;
                                
                                
                                
                                clear rand_tr;
                                if info.no_slicesRAND(group,subset,trial)>1
                                    if trsliceRAND<info.no_slicesRAND(group,subset,trial)
                                        fprintf('%d,',trsliceRAND)
                                        
                                    elseif trsliceRAND==info.no_slicesRAND(group,subset,trial)
                                        fprintf('%d\n',trsliceRAND)
                                        
                                    end
                                end
                                
                                
                                no_frames_slice = ceil((end_frame-start_frame+1)/info.no_slicesRAND(group,subset,trial));
                                frame_range = (trsliceRAND-1)*no_frames_slice+1 : min(trsliceRAND*no_frames_slice+3,end_frame);
                                info.slice_framerangeRAND{group,subset,trial}(trsliceRAND,:) =  ...
                                    [frame_range(1) frame_range(end)];
                                %                                 rand_tr = trPreOrig((trslice-1)*no_frames_slice+1 : min(trslice*(no_frames_slice+2),end_frame) ,:,:);
                                structOut = true;
                                rand_tr=...
                                    idSocial_randomPrepareTrajectories3D(tr, ...
                                    options{group}{subset}{trial}, ...
                                    frame_range,structOut);
                                
                                
                                randtrOrder=...
                                    idSocial_auxiliaries_nearestneighbours(rand_tr);
                                
                                savefast(randOrder_strpath,'randtrOrder')
                               
                                info.no_neighborsRAND(group,subset,trial)=size(rand_tr.Tr,3);
                                if options{group}{subset}{trial}.order_neighbors % Substitute original trajectories by ordered ones
                                    rand_tr = randtrOrder;
                                end
                                savefast(rand_strpath,'rand_tr')
                                
%                                 rand_movementdata{group,subset,trial}=rand_strpath;
%                                 randOrder_movementdata{group,subset,trial}=randOrder_strpath;
                                clear rand_tr;
                                clear randtrOrder;
                                
                            end
                            if info.no_slicesRAND(group,subset,trial) > 1
                                randOrder_movementdata{group,subset,trial}=randOrder_strpathCell;
                                rand_movementdata{group,subset,trial}=rand_strpathCell;
                            else
                                randOrder_movementdata{group,subset,trial}=randOrder_strpath;
                                rand_movementdata{group,subset,trial}=rand_strpath;

                            end
                        end
                        %                     tr = idSocial_filters3D(tr,options{group}{subset}{trial},bodylength,framerate,squeeze(info.circ_roi(group,subset,trial,:)));
                        
                        
                    end
                    
                end
            end
        end
    end
%     fprintf('\n')
end


%
% if info.no_slices(group,subset,trial)>1
%     strpathCell = cell(info.no_slices(group,subset,trial),1);
%     for trslice = 1:info.no_slices(group,subset,trial)
%         no_frames_slice = ceil((end_frame-start_frame+1)/info.no_slices(group,subset,trial));
%         trPre  = trPreOrig( (trslice-1)*no_frames_slice+1 : min(trslice*(no_frames_slice+2),end_frame) ,:,:); % Overlap of 2 so vel and acc can be calculated "seamingless"(?)
%         strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '_' num2str(trslice) '.mat'];
%         savefast(strpath,'trPre');
%         trajectorycell{group,subset,trial}=strpath;
%         strpathCell{trslice} = strpath;
%     end
%     trajectorycell{group,subset,trial} = strpathCell;
%
% else
%     trPre  = trPreOrig;
%     strpath=[temp_savepath 'trajectories' filesep 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
%     savefast(strpath,'trPre');
%     trajectorycell{group,subset,trial}=strpath;
% end

%%% End: Generate randomized trajectory
%%
for group=1:no_groups
    if ~isempty(trajectory{group})
        no_subsets_act=numel(trajectory{group});
        for subset=1:no_subsets_act
            if ~isempty(trajectory{group}{subset})
                no_trials_act=numel(trajectory{group}{subset});
                for trial=1:no_trials_act
                    %             input_data(group,subset,trial).trayectorias=trajectorycell{group,subset,trial};
                    if ~isempty(probtrajectorycell{group,subset,trial})
                        input_data.probtrayectorias{group,subset,trial}=probtrajectorycell{group,subset,trial};
                    else
                        input_data.probtrayectorias{group,subset,trial}=ones(size(trajectorycell{group,subset,trial}));
                    end
                    input_data.movementdata{group,subset,trial}=movementdata{group,subset,trial};
                    input_data.rand_movementdata{group,subset,trial}=rand_movementdata{group,subset,trial};
                    input_data.order_movementdata{group,subset,trial}=order_movementdata{group,subset,trial};
                    input_data.randOrder_movementdata{group,subset,trial}=randOrder_movementdata{group,subset,trial};
                    %                     input_data.options{group}{subset}{trial}.segmentation_path=segmentation_path{group,subset,trial};
                    input_data.options.segmentation_path{group,subset,trial}=segmentation_path{group,subset,trial};
                    
                end
            end
        end
    end
end
input_data(1,1,1).info=info;
input_data(1,1,1).info.trajectory_ver = datestr(now,30);
input_data(1,1,1).info.trajectory_origLocation = trajectory;
input_data(1,1,1).options=options;
savefast([temp_savepath 'input_data.mat'],'input_data')
warning(orig_warning_state)
disp('Preprocessing done.')