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
def_options.project_path=''; % If this is an empty string, no figures will be saved. If it is path, figures will be saved to path\figures and path\figures_latex.
def_options.blpxl=40;
def_options.interpolate_trajectories=false;
def_options.avg_transform2centerofmass=0;
def_options.filter_focal_speedlimits_bl_per_s=[-inf inf];
def_options.filter_neighbor_speedlimits_bl_per_s=[-inf inf];
def_options.filter_neighbor_accelerationlimits_bl_per_s2 = [-inf inf];
def_options.filter_focal_accelerationlimits_bl_per_s2  = [-inf inf];
def_options.filter_distancelimits_bl=[-inf inf];
def_options.filter_focal_rectangularROI = [minx maxx; miny maxy];
def_options.filter_neighbor_rectangularROI = [minx maxx; miny maxy];
def_options.filter_focal_circularROI = [center maxr];
def_options.filter_neighbor_circularROI = [center maxr];
def_options.filter_minProbIdentityAssignment = .8;
def_options.filter_AllMembersPresent = -inf;
def_options.focalReconstruction_minProbIdentity4VelCalcs = [];
def_options.focalReconstruction_minProbIdentityFocal = [];
def_options.focalReconstruction_minProbIdentityNeighbor = [];
def_options.order_neighbors = false;
def_options.framerate=1;
def_options.significance_between_groups=true;
def_options.start_frame=1;
def_options.end_frame=inf;
def_options.group_names=[];
def_options.end_min=inf;
def_options.start_min=0;
def_options.interpolation_mode=[];%'linear';
def_options.smooth_method='moving';
def_options.smooth_degree=1;
def_options.smooth_max_deviation=1;
def_options.smooth_spline_degree = 2;
def_options.smooth_adaptive_noise = 10;
def_options.median_filter=false;
def_options.median_filter_order=1;
def_options.random_data=false;
def_options.random_data_seq_length=60;
def_options.random_data_method='shuffle_frames';
def_options.temp_savepath='';
def_options.no_RandomNeighbors = [];

if nargin < 1
    input_data = def_options;
    return;
end

trajdepth = idSocial_recursiveCellSize(trajectory);
datosegmdepth = idSocial_recursiveCellSize(datosegmlist);

if ~isempty(datosegmlist) && trajdepth ~= datosegmdepth
    warning([mfilename ': Size of trajectory list does not match size of datosegm list.'])
end


options = idSocial_recursiveReadParams(options,input_data,def_options,def_options.act_method);

%% Check input trajectory:
trcell=cell(1,1);
optionscell=cell(1,1);
treat_opt=false;


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
if any(cellfun(@(x) isempty(x),trcell))
    trcell{cellfun(@(x) isempty(x),trcell)}={};
end
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
if isempty(temp_savepath)
    temp_savepath = uigetdir([pwd filesep '~temp'], ...
        sprintf('Please select directory for temporary files\n(idSocial will save results temporarily to the hard disk\n in order to save memory)'));
    temp_savepath = [temp_savepath filesep];
    %     options.temp_savepath = temp_savepath;
elseif ~isempty(temp_savepath) && ischar(temp_savepath)
    try
        if exist(temp_savepath,'dir')~=2
            mkdir(temp_savepath)
        end
    catch
        error([mfilename ': Could not create directory.'])
    end
end
if ~strcmp(temp_savepath(end),filesep)
    temp_savepath=[temp_savepath filesep];
    %         options.temp_savepath = temp_savepath;
    options = idSocial_recursiveSetOptionsInOptionsCell(options,'temp_savepath',temp_savepath);
    
end

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
info.good_frames=NaN(no_groups,no_subsets,no_trials);
info.no_neighbors=NaN(no_groups,no_subsets,no_trials);
info.no_neighborsRAND=NaN(no_groups,no_subsets,no_trials);
info.framerate=ones(no_groups,no_subsets,no_trials);
info.circ_roi=NaN(no_groups,no_subsets,no_trials,4);
info.blpxl=NaN(no_groups,no_subsets,no_trials);
info.filter_AllMembersPresent=NaN(no_groups,no_subsets,no_trials);
info.filter_WorstIndividual=NaN(no_groups,no_subsets,no_trials);
info.filter_AllMembersPresentWorstFocal=NaN(no_groups,no_subsets,no_trials);
info.filter_AllMembersPresentBeforeFilter=NaN(no_groups,no_subsets,no_trials);
info.filter_AllMembersPresentWorstFocalBeforeFilter=NaN(no_groups,no_subsets,no_trials);
info.filter_WorstIndividualBeforeFilter=NaN(no_groups,no_subsets,no_trials);
info.roi=NaN(no_groups,no_subsets,no_trials,2,2);
info.bodylength_mm=NaN(no_groups,no_subsets,no_trials);
info.arena_center=NaN(no_groups,no_subsets,no_trials,2);
info.arena_radius=NaN(no_groups,no_subsets,no_trials);
info.group_name=cell(no_groups,no_subsets);
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
                        (isempty(options{group}{subset}{trial}.blpxl) || ...
                        ischar(options{group}{subset}{trial}.blpxl) && strcmpi(options{group}{subset}{trial}.blpxl,'idTracker'))
                    info.blpxl(group,subset,trial)=nanmean(datosegms{group,subset,trial}.bodylength(:));
                else
                    %                     keyboard
                    disp(['Body length for (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') taken from options, set to ' num2str(options{group}{subset}{trial}.blpxl) ' pxl.'])
                    info.blpxl(group,subset,trial)=options{group}{subset}{trial}.blpxl;
                end
                if isfield(datosegms{group,subset,trial},'bodylength_mm') && ~isempty(datosegms{group,subset,trial}.bodylength-mm) && ...
                        (isempty(options{group}{subset}{trial}.bodylength_mm) || ...
                        ischar(options{group}{subset}{trial}.bodylength_mm) && strcmpi(options{group}{subset}{trial}.blpxl,'idTracker'))
                    info.bodylength_mm(group,subset,trial)=nanmean(datosegms{group,subset,trial}.bodylength_mm(:));
                else
                    %                     keyboard
                    disp(['Body length for (' num2str(group) ',' num2str(subset) ',' num2str(trial) ') taken from options, set to ' num2str(options{group}{subset}{trial}.blpxl) ' pxl.'])
                    info.blpxl(group,subset,trial)=options{group}{subset}{trial}.blpxl;
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
                            info.blpxl(group,subset,trial)=nanmean(datosegms{group,subset,trial}.bodylength(:));
                        catch
                            keyboard
                        end
                        disp([mfilename ': Body length taken from idTracker data file (datosegm.mat) for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(nanmean(datosegms{group,subset,trial}.bodylength)) ' pxl.'])
                        
                    elseif ~ischar(options{group}{subset}{trial}.blpxl)
                        disp([mfilename ': Body length taken from options for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(options{group}{subset}{trial}.blpxl) ' pxl.'])
                        info.blpxl(group,subset,trial)=options{group}{subset}{trial}.blpxl;
                    else
                        disp([mfilename ': Cannot find values for body length. Set to 1 pxl.'])
                        info.blpxl(group,subset,trial)=1;

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
                        info.blpxl(group,subset,trial)=options{group}{subset}{trial}.blpxl;
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
                    if ischar(options{group}{subset}{trial}.blpxl)
                        info.blpxl(group,subset,trial)=NaN;
                        
                    else
                        info.blpxl(group,subset,trial)=options{group}{subset}{trial}.blpxl;
                        disp([mfilename ': Body length taken from options for group ' num2str(group) ', subset ' num2str(subset) ', trial ' num2str(trial) ': Body length = ' num2str(options{group}{subset}{trial}.blpxl) ' pxl.'])
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
%% Try to find out how and where to find the trajectorycell
info.no_groups=no_groups;
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
                            
                            
                            %                             probtrayectorias=[];
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
                            load([segmentation_path{group,trial} 'trajectories.mat']);
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
                        
                        
                        trajectorycell{group,subset,trial}=trayectorias(start_frame:end_frame,:,:);
                        
                        if exist('probtrayectorias','var')==1 && ~isempty(probtrayectorias)
                            probtrajectorycell{group,subset,trial}=probtrayectorias(start_frame:end_frame,:,:);
                        end
                        if exist('probtrajectories','var')==1 && ~isempty(probtrajectories)
                            probtrajectorycell{group,subset,trial}=probtrajectories(start_frame:end_frame,:,:);
                        end
                        if isempty(probtrajectorycell{group,subset,trial})
                            probtrajectorycell{group,subset,trial}=...
                                ones(end_frame-start_frame+1,size(trayectorias,2));
                        end
                        info.no_frames(group,subset,trial)=size(trayectorias(start_frame:end_frame,:,:),1);
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
                    end
                end
                if nargin>2 && ~isempty(options{group}{subset}{trial}) && isfield(options{group}{subset}{trial},'group_names') && ~isempty(options{group}{subset}{trial}.group_names) && ...
                        size(options{group}{subset}{trial}.group_names,2)>=no_groups
                    info.group_name{group}=options{group}{subset}{trial}.group_names{group};
                else
                    info.group_name{group}=['Gr. ' num2str(group)];
                end
%                 keyboard
                info.trials{group,subset}=find(cellfun(@(x) ~isempty(x),trajectory{group}(subset)));
            end % ~isempty(trajectory{group}{subset))
        end
    end % if ~isempty(trajectory{group})
end

%% Check if user wants trajectory reconstruction from segmentation instead of center of mass-trajectorycell
if ~isempty(segmentation_path) && iscell(options{group}) && iscell(options{group}{subset})&& isstruct(options{group}{subset}{trial}) && isfield(options{group}{subset}{trial},'reconstruction_from_segmentation') && options{group}{subset}{trial}.reconstruction_from_segmentation
%     for group=1:no_groups
%         no_subsets_act=numel(trajectory{group});
%         for subset=1:no_subsets_act
%             no_trials_act=numel(trajectory{group}{subset});
%             for trial=1:no_trials_act
%                 input_data(group,subset,trial).options.segmentation_path=segmentation_path{group,subset,trial};
%                 
%                 %                 input_data(group,subset,trial).options{group}{subset}{trial}.segmentation_path=segmentation_path{group,subset,trial};
%             end
%         end
%     end
    input_data.options.segmentation_path=segmentation_path;
    
    input_data(1,1,1).info=info;
    disp('Trajectory reconstruction from segmentation. This may take a while.')
    input_data=idSocial_auxiliaries_analyseBlobShape(input_data);
    % idSocial_auxiliaries_analyseBlobShape will check automatically if calculation has
    % been done already or not. If so, it will simply load the corresponding files
    % instead of re-calculation.
    
    for group=1:no_groups
        no_subsets_act=numel(trajectory{group});
        for subset=1:no_subsets_act
            no_trials_act=numel(trajectory{group}{subset});
            for trial=1:no_trials_act
                trajectorycell{group,subset,trial}=...
                    squeeze(input_data(1,1,1).auxiliaries_analyseBlobShape.central_coordinate(group,subset,trial,start_frame:end_frame,:,1,:));
                bl= input_data(1,1,1).auxiliaries_analyseBlobShape.bodylength(group,subset,trial,:,:,1);
                wd= input_data(1,1,1).auxiliaries_analyseBlobShape.bodywidth(group,subset,trial,:,:,1);
                info.blpxl(group,subset,trial)= nanmean(bl(:));
                info.bwidthpxl(group,subset,trial)= nanmean(wd(:));
            end
        end
    end
end
%% Prepare movement data (calculate velocity, acceleration, ... )

if exist(temp_savepath,'dir')~=7
    mkdir(temp_savepath)
end
for group=1:no_groups
    if ~isempty(trajectory{group})
        no_subsets_act=numel(trajectory{group});
        for subset=1:no_subsets_act
            if ~isempty(trajectory{group}{subset})
                no_trials_act=numel(trajectory{group}{subset});
                for trial=1:no_trials_act
                    if ~isempty(trajectory{group}{subset}{trial})
                        strpath=[temp_savepath 'tr_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                        rand_strpath=[temp_savepath 'trRAND_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                        order_strpath=[temp_savepath 'trORDER_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                        randOrder_strpath=[temp_savepath 'trRANDORDER_' num2str(group) '_' num2str(subset) '_' num2str(trial) '.mat'];
                        
                        if ~isempty(trajectorycell{group,subset,trial})
                            disp(['Preprocessing  (' num2str(group) ',' num2str(subset) ',' num2str(trial) ')'])
                            movementdata{group,subset,trial}=strpath;
                            
                            
                            
                            try
                                probIdCheck = isfinite(options{group}{subset}{trial}.filter_minProbIdentityAssignment) && any(all(isnan(probtrajectorycell{group,subset,trial}) | ...
                                probtrajectorycell{group,subset,trial} < options{group}{subset}{trial}.filter_minProbIdentityAssignment,1));
                            
                                tr=idSocial_prepareTrajectories3D(trajectorycell{group,subset,trial},...
                                probtrajectorycell{group,subset,trial},...
                                options{group}{subset}{trial},info.blpxl(group,subset,trial));
                                if probIdCheck
                                    warning([mfilename ': Assignment probability of at least one individual below threshold (filter_minProbIdentityAssignment) during whole trial.'])
                                end
                            catch
                                keyboard
                            end
                            % Apply filters:
                            if ~isempty(info.blpxl(group,subset,trial)) && ~isnan(info.blpxl(group,subset,trial))
                                bodylength = info.blpxl(group,subset,trial);
                            elseif ~isempty(options{group}{subset}{trial}.blpxl) && ~isnan(options{group}{subset}{trial}.blpxl)
                                bodylength = options{group}{subset}{trial}.blpxl;
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
                            [good_frames_before_filter, all_members_frames_before_filter, all_members_present_for_worst_focal_before_filter] = idSocial_auxiliaries_checkDataQuality(tr.trajectory);
                            catch
                                keyboard
                            end
                            tr.trajectory = idSocial_filters3D(tr.trajectory,options{group}{subset}{trial},bodylength,framerate,squeeze(info.circ_roi(group,subset,trial,:)));
                            %                             good_frames=squeeze(sum(~isnan(tr.trajectory(:,:,:,1)),1)./size(tr.trajectory,1));
                            %                             good_frames=round(min(good_frames(:))*100);
                            %
                            %                             trPerm=squeeze(permute(tr.trajectory(:,:,:,1),[1,4,2,3]));
                            %                             trPerm = trPerm(:,logical(eye(size(tr.trajectory,2),size(tr.trajectory,2))));
                            %
                            %                             all_members_frames = sum(sum(~isnan(trPerm),2)==size(tr.trajectory,2))./size(tr.trajectory,1);
                            %                             all_members_frames = all_members_frames*100;
                            
                            [good_frames, all_members_frames,all_members_present_for_worst_focal ] = idSocial_auxiliaries_checkDataQuality(tr.trajectory);
                            if good_frames<70
                                warning(['% of remaining frames after filtering (Min. for any individual): ' num2str(good_frames)])
                            else
                                disp(['% of remaining frames after filtering (Min. for any individual): ' num2str(good_frames)])
                            end
                            if all_members_frames<70
                                warning(['% of frames in which all individuals are present: ' num2str(all_members_frames,'%.1f')])
                            else
                                disp(['% of frames in which all individuals are present: ' num2str(all_members_frames,'%.1f')])
                            end
                            
                            
                            info.filter_AllMembersPresentBeforeFilter(group,subset,trial) = all_members_frames_before_filter;
                            info.filter_WorstIndividualBeforeFilter(group,subset,trial) = good_frames_before_filter;
                            info.filter_AllMembersPresentWorstFocalBeforeFilter(group,subset,trial) = all_members_present_for_worst_focal_before_filter;
                            info.filter_AllMembersPresent(group,subset,trial) = all_members_frames;
                            info.filter_WorstIndividual(group,subset,trial) = good_frames;
                            info.filter_AllMembersPresentWorstFocal(group,subset,trial) = all_members_present_for_worst_focal;
                            % End: Apply filters
                            
                            order_movementdata{group,subset,trial}=order_strpath;
                            trOrder=tr;
                            trOrder.trajectory=idSocial_auxiliaries_nearestneighbours(tr.trajectory);
                            savefast(order_strpath,'trOrder')
                            
                            if isstruct(tr.trajectory)
                                info.no_focals(group,subset,trial)=size(tr.trajectory.Tr,2);
                                info.no_neighbors(group,subset,trial)=size(tr.trajectory.Tr,3);
                                info.good_frames(group,subset,trial)=good_frames;
                            else
                                info.no_focals(group,subset,trial)=size(tr.trajectory,2);
                                info.no_neighbors(group,subset,trial)=size(tr.trajectory,3);
                                info.good_frames(group,subset,trial)=good_frames;
                            end
                            
                            if random_data && info.no_focals(group,subset,trial)>1
                                rand_movementdata{group,subset,trial}=rand_strpath;
                                randOrder_movementdata{group,subset,trial}=randOrder_strpath;
                                clear rand_tr;
                                
                                rand_tr=...
                                    idSocial_randomPrepareTrajectories3D(tr.trajectory,options{group}{subset}{trial});
                                
                                randtrOrder=...
                                    idSocial_auxiliaries_nearestneighbours(rand_tr);
                                
                                
                                savefast(randOrder_strpath,'randtrOrder')
                                
                                info.no_neighborsRAND(group,subset,trial)=size(rand_tr,3);
                            end
                            
                            if options{group}{subset}{trial}.order_neighbors % Substitute original trajectories by ordered ones
                                tr = trOrder;
                                rand_tr = randtrOrder;
                            end
                            clear trOrder randtrOrder;
                            
                            savefast(strpath,'tr')
                            clear tr;
                            if random_data && info.no_focals(group,subset,trial)>1
                                savefast(rand_strpath,'rand_tr')
                                clear rand_tr;
                            end
                            
                        else
                            trajectorycell{group,subset,trial}=NaN(end_frame-start_frame+1,2,2);
                            rand_trajectorycell{group,subset,trial}=NaN(end_frame-start_frame+1,2,2);
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
                            info.blpxl(group,subset,trial)=-1;
                            info.no_focals(group,subset,trial)=-1;
                            info.no_neighbors(group,subset,trial)=-1;
                            info.no_neighborsRAND(group,subset,trial)=-1;
                            info.duration(group,subset,trial)=0;
                            
                        end
                    end
                end
                if nargin>2 && ~isempty(options{group}{subset}{trial}) && isfield(options{group}{subset}{trial},'group_names') && ~isempty(options{group}{subset}{trial}.group_names) && ...
                        size(options{group}{subset}{trial}.group_names,2)>=no_groups
                    info.group_name{group}=options{group}{subset}{trial}.group_names{group};
                else
                    info.group_name{group}=['Gr. ' num2str(group)];
                end
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
            end %  ~isempty(trajectory{group}{subset})
        end
    end %  ~isempty(trajectory{group})
end
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
input_data(1,1,1).info.trajectory_ver = datestr(now,30);
input_data(1,1,1).info=info;
input_data(1,1,1).options=options;
savefast([temp_savepath 'input_data'],'input_data')
warning(orig_warning_state)