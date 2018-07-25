function trajectoryLocation = idSocial_recursiveBuildTrajectoryLocationFromList(trajectoryLocation,pathList)



if size(pathList{2},2) > 1
    act_idx=pathList{2}(1);
%     fprintf('%d,',act_idx)
    if numel(trajectoryLocation)<act_idx
        trajectoryLocation{act_idx}=[];
    end
    pathList{2}=pathList{2}(2:end);
    try
    trajectoryLocation{act_idx} = ...
        idSocial_recursiveBuildTrajectoryLocationFromList(trajectoryLocation{act_idx},pathList);
    catch
        keyboard
    end
else
    trajectoryLocation{pathList{2}(1),1} = pathList{1};
    return
end
    