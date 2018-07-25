function [patterncell optionscell] = idSocial_recursiveOptionsCellSub2OptionsCell(optionscell,patterncell)
% This function takes the cell 'patterncell', which is a 'cell of cell' of the form patterncell{1}=cell(..,..);
% patterncell{1}{1}=cell(...,...); etc. and copies information
% from 'optionscell', which has the same structure but can
% have smaller depth downstream in a way that if
% optionscell{1}{1} = options, patterncell{1}{1}{:}=options.
% ATTENTION: The first output patterncell returns the new
% optionscell!


if all(cellfun(@(x) isa(x,'cell')|| isempty(x) ,patterncell))
    for k = 1:numel(patterncell)
        if ~isempty(patterncell{k})
            if isa(optionscell,'cell')
                
                [patterncell{k} optionscell{k}] = idSocial_recursiveOptionsCellSub2OptionsCell(optionscell{k},patterncell{k});
                
            else
                
                [patterncell{k} optionscell] = idSocial_recursiveOptionsCellSub2OptionsCell(optionscell,patterncell{k});
                
            end
        end
    end
else
    options = optionscell;
    for k=1:numel(patterncell)
        if iscell(options) && numel(options)==1
            patterncell{k}=options{1};
        elseif iscell(options) && numel(options)==numel(patterncell)
            patterncell{k}=options{k};
        elseif isstruct(options)
            patterncell{k}=options;
        end
    end
end


end


