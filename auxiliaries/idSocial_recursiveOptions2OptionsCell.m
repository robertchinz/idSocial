function optionscell = idSocial_recursiveOptions2OptionsCell(optionscell,options)

if isa(optionscell,'cell')
    for k = 1:numel(optionscell)
 
        optionscell{k} = idSocial_recursiveOptions2OptionsCell(optionscell{k},options);
    end
else
    optionscell = options;
    return
end
end
