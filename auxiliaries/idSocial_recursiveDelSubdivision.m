function [trajectory_location, act_div] = idSocial_recursiveDelSubdivision(trajectory_location,div_id,act_div)


if nargin < 3 || isempty(act_div)
    act_div=1;
end

if div_id==act_div
    try
        trajectory_location=vertcat(trajectory_location{:});
    catch
        keyboard
    end
    return
%     act_div=act_div-1;
elseif div_id>act_div
    for k=1:size(trajectory_location,1)
        [trajectory_location{k},  ~] = ...
            idSocial_recursiveDelSubdivision(trajectory_location{k},div_id,act_div+1);
    end
end



