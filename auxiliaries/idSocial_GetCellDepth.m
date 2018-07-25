function [elist_all, pathlist_all]= idSocial_GetCellDepth(trajectory_location)

elist_all=[];
pathlist_all=[];
while ~isempty(trajectory_location)
    
    [trajectory_location,elist,pathlist]= idSocial_recursiveGetCellDepth(trajectory_location);
    if ~isempty(elist)
      elist_all=vertcat(elist_all,{elist});
    end
    if ~isempty(pathlist)
      pathlist_all=vertcat(pathlist_all,{pathlist});
    end
end

