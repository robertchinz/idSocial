function input_data = idSocial_auxiliaries_removeSet(input_data,act_method,act_set)

if strcmp(act_method(end-10:end), '4Statistics')
    principal_method = act_method(1:end-11);
else
    principal_method = act_method;
end

if isfield(input_data,principal_method)
    fn = fieldnames(input_data.(principal_method));
    maxSet = -1;
    for k=1:size(fn,1)
        if ~isempty(strfind(fn{k},'Set'))
            maxSet = max(maxSet,str2double(strrep(fn{k},'Set','')));
        end
    end
    
    if nargin < 3 || isempty(act_set)
        act_set = maxSet;
    end
    
    if ischar(act_set) && strcmpi(act_set,'all')
        act_set = 1:maxSet;
    end
    
    for as = act_set
        if as<=maxSet
            input_data.(principal_method) = rmfield(input_data.(principal_method),['Set' num2str(as)]);
        else
            warning([mfilename ': Set' num2str(as) ' does not exist.'])
        end
    end

else
    disp([mfilename ': Field ' principal_method ' does not exist.'])
end
