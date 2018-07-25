function trajectoryLocation=idSocial_buildTrajectoryLocationTreeFromList(pathList)


no_entries=size(pathList,1);
trajectoryLocation=cell(1,1);
for file=1:no_entries
    act_idx = pathList{file,2}(1);
    pathList{file,2} = pathList{file,2}(2:end);
    try
        if numel(trajectoryLocation)<act_idx
         trajectoryLocation{act_idx}=[];
        end
   trajectoryLocation{act_idx} = ...
       idSocial_recursiveBuildTrajectoryLocationFromList(trajectoryLocation{act_idx},pathList(file,:)); 
    catch
        keyboard
    end
%    fprintf('%s\r\n','')
end