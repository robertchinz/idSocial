function trajectory = idSocial_recursiveCorrectTrajectoryStructure(trajectory,level)
if nargin < 2 || isempty(level)
    level=0;
end


no_nodes = numel(trajectory);

if isa(trajectory,'char')
    trajectory = cellstr(trajectory);
    return
end
level = level + 1;
if level < 3 
    for k = 1:no_nodes
        try
        trajectory{k} = idSocial_recursiveCorrectTrajectoryStructure(trajectory{k},level);
        catch
            keyboard
        end
    end
else
    return
end
end
