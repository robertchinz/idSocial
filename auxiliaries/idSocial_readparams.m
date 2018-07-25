function [input_data,options_out]=idSocial_readparams(input_data,options,def_options,act_method,post_complete)
% Sets the parameters for the current method:
% Parameters have the following sort of hierarchy:
% input_data(1,1).options contains options for all experiments,
% input_data(group,trial).options contains specific options for each trial (and can be empty),
% input_data(group,trial).(current_method).options contains options only necessary and
% valid for the current method, and
% input_data(group,trial).options overwrites input_data(1,1).options (where 'overwrite'
% means that the corresponding options are taken for the current calculations, the
% actual values stay untouched!). The values which are valid for the current method are
% then saved into input_data(group,trial).(current_method).options.
% Usually, each function contains some default values, which are also overwritten by any
% of the previous.
% NOTE that in theory for each group and trial, input_data(group,trial).(current_method).options
% is different, but most functions which have to act on data from various groups and trials
% take options from input_data(1,1).(current_method).options.
%
% def_options_flds=fieldnames(def_options);
% no_fields_def=size(def_options_flds,1);

if nargin <4
    act_method = '';
end

if nargin <5 || isempty(post_complete)
    post_complete = false;
end

if size(def_options,2)==2
    def_options_tips = def_options(2);
    def_options = def_options(1);
end

% For "combo box" options which in the def. options section of each
% function is defined as a cell, pick the first element as default.
if ~isempty(def_options)
    options_def_flds=fieldnames(def_options);
    no_fields_def=size(options_def_flds,1);
    for def_fld=1:no_fields_def
        if  iscell(def_options.(options_def_flds{def_fld})) && all(cellfun(@(x) isequal(class(x),class(def_options.(options_def_flds{def_fld}){1})),def_options.(options_def_flds{def_fld}))) && ...
                any(size(def_options.(options_def_flds{def_fld}))==1)
            def_options.(options_def_flds{def_fld}) = def_options.(options_def_flds{def_fld}){1};
        end
    end
end

% For "combo box" options which in the options section of each
% function is defined as a cell, pick the last element as (which stores the choice).
if ~isempty(options)
    options_flds=fieldnames(options);
    no_fields=size(options_flds,1);
    for def_fld=1:no_fields
        if  iscell(options.(options_flds{def_fld})) && ...
                ~isempty(options.(options_flds{def_fld})) && ...
                all(cellfun(@(x) isequal(class(x),class(options.(options_flds{def_fld}){1})),options.(options_flds{def_fld})(1:end-1))) && ...
                iscell(options.(options_flds{def_fld})(end)) && ...
                any(size(options.(options_flds{def_fld}))==1) && ...
                ~isempty(options.(options_flds{def_fld}){end}) && ...
                all(cellfun(@(x) isequal(class(x),class(options.(options_flds{def_fld}){end})),options.(options_flds{def_fld})(1:end-1)))
%             optNo = options.(options_flds{def_fld}){end}{1};
            options.(options_flds{def_fld}) =  options.(options_flds{def_fld}){end};%options.(options_flds{def_fld}){optNo};
        end
    end
end

if ~isempty(options) && ~isempty(def_options)
    %% Compare input_data(1,1).options with default options:
    % Take options from input_data(1,1).options, and if some option is missing take the
    % default options from def_options as given in the header of the function.
    options_def=def_options;
    options_def_flds=fieldnames(options_def);
    options_flds=fieldnames(options);
    no_fields_def=size(options_def_flds,1);
    no_fields=size(options_flds,1);
    options_all=options;
    
    options_out=options_def;
    
    
    
    for def_fld=1:no_fields_def
        if isfield(options_all,options_def_flds{def_fld}) %&& ~isempty(options_all.(options_def_flds{def_fld}))
            options_out.(options_def_flds{def_fld})=options_all.(options_def_flds{def_fld});
        else
            options_out.(options_def_flds{def_fld})=options_def.(options_def_flds{def_fld});
        end
    end
    
    if post_complete
        for fld=1:no_fields
            if ~isfield(options_out,options_flds{fld})
                options_out.(options_flds{fld})=options.(options_flds{fld});
            end
        end
    end
    
    % %% Compare input_data(group,trial).options and input_data(group,trial).(current_method).options with defaults/input_data(1,1).options
    % for group=1:no_groups
    %     for trial=1:no_trials
    %
    %
    %         options_grtr=input_data(group,trial).options;
    %         for def_fld=1:no_fields_def
    %             if isfield(options_grtr,options_def_flds{def_fld}) && ~isempty(options_grtr.(options_def_flds{def_fld}))%&& ...
    % %                         ~(isfield(options_out,def_options_flds{def_fld}) || isempty(options_out.def_options_flds{def_fld}))
    %                 options_out.(options_def_flds{def_fld})=options_grtr.(options_def_flds{def_fld});
    %             end
    %         end
    %
    % %         % Complete with input_data(group,trial).(act_method).options in case it provides more options which
    % %         % are not set by input_data(1,1).options, input_data(group,trial).options or
    % %         % def_options.
    % %         if isfield(input_data(group,trial),act_method) && isfield(input_data(group,trial).(act_method),'options')
    % %             options_fct=input_data(group,trial).(act_method).options;
    % %             for def_fld=1:no_fields_def
    % %                 if isfield(options_fct,def_options_flds{def_fld}) && ~isempty(options_fct.(def_options_flds{def_fld})) && ...
    % %                         ~(isfield(options_out,def_options_flds{def_fld}) || isempty(options_out.def_options_flds{def_fld}))
    % %                     options_out.(def_options_flds{def_fld})=options_fct.(def_options_flds{def_fld});
    % %                 end
    % %             end
    % %         end
    %         input_data(group,trial).(act_method).options=options_out;
    %     end
    % end
elseif isempty(options)
    options_out = def_options;
elseif isempty(def_options)
    options_out = options;
end
% options_out=input_data(1,1).options;
if ~isempty(input_data)
    input_data(1,1).(act_method).options=options_out;
else
    input_data=[];
end