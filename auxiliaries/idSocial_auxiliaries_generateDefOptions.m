function opts = idSocial_auxiliaries_generateDefOptions(opts)

if ~isempty(opts)
    fn = fieldnames(opts);
    
    for k = 1:size(fn,1)

        if iscell(opts.(fn{k})) && all(all(cellfun(@(x) ischar(x),opts.(fn{k})))) && ...
                ~all(cellfun(@(x) isempty(x),opts.(fn{k})))
%             try
            opts.(fn{k}) = opts.(fn{k}){1};
%             catch
%                 keyboard
%             end
        end

    end
end