function parameters = handleComboMenus(parameters)
% Handle combo menus
% parameters = gui.defPreOpts;
optNames = fieldnames(parameters(1));
optVals = struct2cell(parameters(1));
isCombo = cellfun(@(x) iscell(x)&size(x,2)>1,optVals);
optCombos = optVals(isCombo);
optComboName = optNames(isCombo);

for k=1:sum(isCombo)
%     all(cellfun(@(x) ischar(x),optCombos{k}(1:end-1))) && ...

    % All entries but the last one are of the same class, the last one is a
    % cell containing the chosen option. In this case set parameter to this
    % chosen option
    if all(cellfun(@(x) isequal(class(x),class(optCombos{k}{1})),optCombos{k}(1:end-1))) && ...
            iscell(optCombos{k}{end}) && all(cellfun(@(x) isequal(class(x),class(optCombos{k}{end}{1})),optCombos{k}(1:end-1)))
        
        optCombos{k} = optCombos{k}{end};%optCombos{k}{optCombos{k}{end}{1}};
        parameters(1).(optComboName{k}) = optCombos{k};
        
    % Not all entries (but last one) are the same, and there are some cells    
    elseif any(cellfun(@(x) iscell(x),optCombos{k})) % combined "combo -non-combo" menu
        
        isNestedCombo = find(cellfun(@(x) iscell(x),optCombos{k}));
%         optCombosNested = optCombos{k}(isNestedCombo);
        for k2=isNestedCombo
            if all(cellfun(@(x) isequal(class(x),class(optCombos{k}{k2}{1})),optCombos{k}{k2}(1:end-1))) && ...
                iscell(optCombos{k}{k2}{end}) && all(cellfun(@(x) isequal(class(x),class(optCombos{k}{k2}{end}{1})),optCombos{k}{k2}(1:end-1)))
                
%                 iscell(optCombos{k}{k2}{end}) && ischar(optCombos{k}{k2}{end}{1})
                optCombos{k}{k2} = optCombos{k}{k2}{end}{1};%optCombos{k}{optCombos{k}{end}{1}};
            else
                optCombos{k}{k2} = optCombos{k}{k2}{1};
            end
        end
        parameters(1).(optComboName{k}) = optCombos{k};
        
    % Not all elements are the same, but none is a variable with >1 dimensions: Keep the parameter
    % Example: filter_focal_circularROI = {[-Inf]    [Inf]    'BL'}
    elseif ~all(cellfun(@(x) isequal(class(x),class(optCombos{k}{1})),optCombos{k})) && ...
            all(cellfun(@(x) numel(x)<=1,optCombos{k}) | cellfun(@(x) ischar(x),optCombos{k})) 
        parameters(1).(optComboName{k}) = optCombos{k};
    else
        parameters(1).(optComboName{k}) = optCombos{k}{1};
    end
end
end