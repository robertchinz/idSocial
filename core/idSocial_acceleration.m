function [mean_acceleration, acceleration_per_frame]=idSocial_acceleration(trajectory,framerate,bodylength)
%   Acceleration from trajectory.
%
%   mean_acceleration = ...
%           IDSOCIAL_ACCELERATION(TR) 
%   calculates the average magnitude of acceleration over 
%   time for trajectory TR.
%
%   TR is either an N-D array with dimensions 
%   #Frames-by-#Individuals-by-#Dimensions, in which case 
%   mean_acceleration is a 1-by-#Individuals vector, or an 
%   N-D array with dimensions 
%   #Frames-by-#Individuals-by-#Individuals-by-#Dimensions,
%   in which case mean_acceleration is an 
%   #Individuals-by-#Individuals vector. 
%   The second case enables non-symmetric filtering, e.g., 
%   when values of the focal individual are filtered 
%   depending on position, speed, etc. of each neighbor 
%   individual (see idSocial_prepareTrajectories3D). 
%   #Dimensions can be 2 or 3.
%
%   mean_acceleration = ...
%           IDSOCIAL_ACCELERATION(TR,FRAMERATE,BODYLENGTH) 
%   scales the result according to the given framerate
%   (normally the framerate at which the video has been
%   recorded) and the body length in pixels of either each
%   individual in a 1-by-#Individuals vector or the same 
%   scalar value for all individuals. The default values 
%   are:
%   FRAMERATE=1, BODYLENGTH=1 
%   (resulting in units of 1 pixel/frame^2)
%
%   [mean_acceleration, acceleration_per_frame] = ...
%           IDSOCIAL_ACCELERATION(TR) 
%   in addition returns an array acceleration_per_frame of 
%   dimensions #Frames-by-#Individuals containing the
%   magnitude of acceleration for each time step and each
%   individual. 
%
%   See also idSocial_prepareTrajectories3D, 
%   idSocial_accelerationDistribution, idSocial_speed

%   2014 Robert C. Hinz, Gonzalo G. de Polavieja, Consejo Superior de Investigaciones Científicas

if nargin<2 || isempty(framerate)
    framerate=1;
end
if nargin<3 || isempty(bodylength)
    bodylength=1;
end

tr_size=size(trajectory);
no_arraydims=ndims(trajectory);


acc=cat(1,diff(trajectory,2,1),NaN([1 tr_size(2:end)]));
acc_magnitude=sqrt(sum(acc.^2,no_arraydims))./bodylength*framerate^2;

mean_acceleration=squeeze(nanmean(acc_magnitude,1));
acceleration_per_frame={permute(acc_magnitude,[2:no_arraydims-1 1])};



