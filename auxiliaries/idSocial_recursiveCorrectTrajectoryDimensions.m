function trajectory = idSocial_recursiveCorrectTrajectoryDimensions(trajectory)
% if nargin < 2 || isempty(level)
%     level=1;
% end

if ~isempty(trajectory) && isa(trajectory,'cell') && ~all(cellfun(@(x) ischar(x),trajectory))
    no_nodes = numel(trajectory);
    if size(trajectory,1)>1 && size(trajectory,2)==1
        trajectory = trajectory';
    end
%     if level==2 && size(trajectory,1)>1 && size(trajectory,2)==1
%         trajectory = trajectory';
%     end
    for k = 1:no_nodes
%         level = level + 1;
        trajectory{k} = idSocial_recursiveCorrectTrajectoryDimensions(trajectory{k});
    end
else
    return
end
end
