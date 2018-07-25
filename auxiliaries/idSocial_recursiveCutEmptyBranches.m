function [trajectory_location]= idSocial_recursiveCutEmptyBranches(trajectory_location)


count=1;
for k=1:size(trajectory_location,1)
    try
        if  isa(trajectory_location,'cell') && isa(trajectory_location{count,1},'cell') && trajectory_location{count,2}>0
            
            [trajectory_location{count,1}]= ...
                idSocial_recursiveCutEmptyBranches(trajectory_location{count,1});
            count=count+1;
        elseif isa(trajectory_location,'cell') && isa(trajectory_location{count,1},'cell') && trajectory_location{count,2}==0
            trajectory_location(count,:)=[];
        elseif isa(trajectory_location,'cell') && ~isa(trajectory_location{count,1},'cell')
            trajectory_location=trajectory_location{count,1};
             count=count+1;
        elseif ~isa(trajectory_location,'cell')
          
%             1
        end
    catch
        keyboard
    end
   
end
trajectory_location(:,2)=[];

