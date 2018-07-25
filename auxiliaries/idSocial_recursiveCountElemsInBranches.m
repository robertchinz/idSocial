function [trajectory_location]= idSocial_recursiveCountElemsInBranches(trajectory_location)


if isa(trajectory_location,'cell') 
   for k=1:size(trajectory_location,1)
         [trajectory_location{k}]= ...
           idSocial_recursiveCountElemsInBranches(trajectory_location{k});
%         trajectory_location{k,2}=
   end
   %    trajectory_location=trajectory_location(cellfun(@(x) x{1},trajectory_location));
   
   
   for k=1:size(trajectory_location,1)
       if isa(trajectory_location{k,1},'cell')
           trajectory_location{k,2} =  sum(vertcat(trajectory_location{k,1}{:,2}));
%            trajectory_location{k,1}{1};
           
       end
   end
   
elseif ~isa(trajectory_location,'cell') && ~isempty(trajectory_location)
    trajectory_location = {trajectory_location, true};
elseif ~isa(trajectory_location,'cell') && isempty(trajectory_location)
    trajectory_location = {[],false};
end