function [trajectory_location,elist,pathlist]= idSocial_recursiveGetCellDepth(trajectory_location,elist,pathlist)

if nargin<2 || isempty(elist)
    elist=[];
end
if nargin<2 || isempty(pathlist)
    pathlist=[];
end

try
if  isa(trajectory_location{1,1},'cell') && ~isempty(trajectory_location{1,1})...
        && isa(trajectory_location{1,1}{1,1},'cell') && ~isempty(trajectory_location{1,1}{1})
    
    try
        [trajectory_location{1,1},elist,pathlist]= ...
            idSocial_recursiveGetCellDepth(trajectory_location{1,1},[elist trajectory_location{1,2}],pathlist);
    catch
        keyboard
    end
  
    
elseif isa(trajectory_location{1,1},'cell')  && ~isempty(trajectory_location{1,1}) && ...
        (~isa(trajectory_location{1,1}{1,1},'cell')  && trajectory_location{1,2}==1)
    %             elist = vertcat(elist,new_list);
    %     elist=;
    disp(trajectory_location{1,1}{1});
    pathlist=trajectory_location{1,1}{1};
    trajectory_location(1,:)=[];
    return;
elseif isa(trajectory_location{1,1},'cell') && ~isempty(trajectory_location{1,1}) && ...
        (isempty(trajectory_location{1,1}{1,1}))% && trajectory_location{1,1}{1,2}==1)
    %             elist = vertcat(elist,new_list);
    %     elist=;
    trajectory_location{1,1}(1,:)=[];
    elist=[];
    pathlist=[];
    return;
elseif isa(trajectory_location{1,1},'cell') && isempty(trajectory_location{1,1})
    %             elist = vertcat(elist,new_list);
    %     elist=;
    trajectory_location(1,:)=[];
    elist=[];
    pathlist=[];
    return;
end
catch
    keyboard
end






