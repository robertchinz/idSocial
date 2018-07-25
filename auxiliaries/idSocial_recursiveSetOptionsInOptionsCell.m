function optionscell = idSocial_recursiveSetOptionsInOptionsCell(optionscell,optentry_field,optentry_value)

if isa(optionscell,'cell')
    for k = 1:numel(optionscell)
        if ~isempty(optionscell{k})
            optionscell{k} = idSocial_recursiveSetOptionsInOptionsCell(optionscell{k},optentry_field,optentry_value);
        end
    end
else
    
    optionscell.(optentry_field)=optentry_value;
    return
end
end
