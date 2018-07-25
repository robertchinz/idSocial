function propString = idSocial_auxiliaries_optstruct2string(input_data,funcstring,act_fig)
% act_set = input_data.(funcstring).(fn_figs{act_fig}).act_set;
try
if isfield(input_data.(funcstring).(['Fig' num2str(act_fig)]),'act_set')
    act_set = input_data.(funcstring).(['Fig' num2str(act_fig)]).act_set;
    act_set = unique(floor(act_set)); % Handle rand. controls. 

% propStringCell = cell(2,numel(act_set));

% for as = 1:numel(act_set)
%     if isfield(input_data.(funcstring),['Set' num2str(as)])
%         propStringCell{1,as} = evalc(['disp(input_data.' funcstring '.Set' num2str(as) '.options)']);
%         propStringCell{2,as} = evalc(['disp(input_data.' funcstring '.Set' num2str(as) '.plot_mode)']);
%     end
% end
% 
% propString = [];
% for as = 1:numel(act_set)
%     propString = [propString propStringCell{1,as}(:)' propStringCell{2,as}(:)'];
% end

%% Start new


allOptCell = cell(1,numel(act_set));
allPmCell = cell(1,numel(act_set));


for as = 1:numel(act_set)
    if isfield(input_data.(funcstring),['Set' num2str(act_set(as))])
        allOptCell{as} = input_data.(funcstring).(['Set' num2str(act_set(as))]).options;
        allPmCell{as} = input_data.(funcstring).(['Set' num2str(act_set(as))]).plot_mode;
    end
end

fnamesOpts = cellfun(@(x) fieldnames(x),allOptCell,'UniformOutput',false);
fnamesPm = cellfun(@(x) fieldnames(x),allPmCell,'UniformOutput',false);

allFnamesOpts = unique(vertcat(fnamesOpts{:}));
allFnamesPm = unique(vertcat(fnamesPm{:}));


def_optionsAll = eval(['idSocial_' funcstring]);
def_optionsFn = fieldnames(def_optionsAll);

allFnamesOpts = allFnamesOpts(~strcmp(allFnamesOpts,'act_method'));

allFnamesOpts = allFnamesOpts(ismember(allFnamesOpts,def_optionsFn));

no_opts = numel(allFnamesOpts);
optVals = cell(no_opts,numel(act_set));
for ao = 1:no_opts
    optVals{ao,1} = allFnamesOpts{ao};
    for as=1:numel(act_set)
        if isfield(input_data.(funcstring).(['Set' num2str(act_set(as))]).options,allFnamesOpts{ao})
            try
                optVals{ao,as+1} = num2str(input_data.(funcstring).(['Set' num2str(act_set(as))]).options.(allFnamesOpts{ao}));
            catch
                optVals{ao,as+1} = input_data.(funcstring).(['Set' num2str(act_set(as))]).options.(allFnamesOpts{ao});
            end
        end
    end
end


optDiffs = false(no_opts,1);

for fn = 1:no_opts
    if all(cellfun(@(x) isfield(x,allFnamesOpts{fn}),allOptCell)) % Field is in all opts: Check if values are different.
        act_vals = cellfun(@(x) x.(allFnamesOpts{fn}),allOptCell,'UniformOutput',false);
        optDiffs(fn) = ~all(cellfun(@(x) isequal(x,act_vals{1}),act_vals));
    else % Field is NOT in all opts
        optDiffs(fn) = true;
    end
end
    
optDiffsName = allFnamesOpts(optDiffs);
optDiffsVals = cell(numel(optDiffsName),numel(act_set));

for fnDiff = 1:numel(optDiffsName)
    optDiffsVals(fnDiff,:) = cellfun(@(x) x.(optDiffsName{fnDiff}),allOptCell,'UniformOutput',false);
    
end

no_pm = numel(allFnamesPm);
pmDiffs = false(no_pm,1);
for fn = 1:no_pm
    if all(cellfun(@(x) isfield(x,allFnamesPm{fn}),allPmCell)) % Field is in all Pm: Check if values are different.
        act_vals = cellfun(@(x) x.(allFnamesPm{fn}),allPmCell,'UniformOutput',false);
        pmDiffs(fn) = ~all(cellfun(@(x) isequal(x,act_vals{1}),act_vals));
    else % Field is NOT in all Pm
        pmDiffs(fn) = true;
    end
end

pmDiffsName = allFnamesPm(pmDiffs);
pmDiffsVals = cell(numel(pmDiffsName),numel(act_set));

for fnDiff = 1:numel(pmDiffsName)
    pmDiffsVals(fnDiff,:) = cellfun(@(x) x.(pmDiffsName{fnDiff}),allPmCell,'UniformOutput',false);
    
    % Special for data_filter:
    if strcmpi('data_filter',pmDiffsName{fnDiff})
        for as = 1:numel(act_set)
            pmDiffsVals{fnDiff,as} = input_data.(funcstring).(['Set' num2str(act_set(as))]).output_plot.data_string{1}(5:end);
        end
    end
    % End data_filter
end


% propString = vertcat(evalc('disp(vertcat([{''Set ''} num2cell(act_set)],[optDiffsName optDiffsVals ],[pmDiffsName pmDiffsVals ],optVals))'));

propString = vertcat(evalc('disp(vertcat([{''Set ''} num2cell(act_set)],[optDiffsName optDiffsVals ],[pmDiffsName pmDiffsVals ],optVals))'));

%% End new
else
    propString = 'Sorry, no info.';
end
catch exception
    msgText = getReport(exception);
    propString = ['There has been some error: ' msgText];
end