function a = idSocial_auxiliaries_vector2angle2pi(v1,v2)
        a = atan2(-v1(:,1).*v2(:,2)-v2(:,1).*(-v1(:,2)), ...
            (-v1(:,1)).*v2(:,1)+(-v1(:,2)).*v2(:,2));
end
    
% This might be wrong (i.e., wrong sign of one of the vectors), but I do not want to change it because it works for Onset3_simulation4modelCreateTrajectory 
% The actual formula is:
% a = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
% http://www.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
% " between the vectors as measured in a counterclockwise direction from v1
% to v2"