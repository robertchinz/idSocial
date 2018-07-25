function [optionscell, value]= idSocial_recursiveGetOptionsFromOptionsCell(optionscell,optentry_field)
value=[];

if isempty(value) && isa(optionscell,'cell')
    for k = 1:numel(optionscell)
        if ~isempty(optionscell{k})
            [optionscell{k} value]= idSocial_recursiveGetOptionsFromOptionsCell(optionscell{k},optentry_field);
        end
    end
else
    if isfield(optionscell,optentry_field)
        value = optionscell.(optentry_field);
    else
        value = [];
    end
    return
end
end
