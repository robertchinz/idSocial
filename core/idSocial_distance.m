function [distnce,nndist,fardist]=idSocial_distance(trajectory,bodylength)
%   Inter-indivual distance from trajectory.
%
%   distnce = ...
%           IDSOCIAL_DISTANCE(TR) 
%   calculates the distance between all possible pairs of
%   individuals in each frame. 
%
%   TR is either an N-D array with dimensions 
%   #Frames-by-#Individuals-by-#Dimensions, 
%   or an N-D array with dimensions 
%   #Frames-by-#Individuals-by-#Individuals-by-#Dimensions.
%   (The second case enables non-symmetric filtering, e.g., 
%   when values of the focal individual are filtered 
%   depending on position, speed, etc. of each neighbor 
%   individual (see idSocial_prepareTrajectories3D)).
%   In both cases, distnce is a 
%   #Frame-by-#Individuals-by-#Individuals array.
%   #Dimensions can be 2 or 3.
%
%   distnce = ...
%           IDSOCIAL_DISTANCE(TR,BODYLENGTH) 
%   scales the result according to the body length in pixels
%   of either each individual in a 1-by-#Individuals vector 
%   or the same scalar value for all individuals. 
%   The default is:
%   BODYLENGTH=1 
%   (resulting in units of 1 pixel)
%
%   [distnce, nndist, fardist] = ...
%           IDSOCIAL_DISTANCE(TR) 
%   in addition returns an array nndist (fardist) of 
%   dimensions #Frame-by-#Individuals-by-#Individuals array,
%   containing only the distance to the nearest (furthest)
%   neighbor for each time step and each individual. 
%   Elements for not-nearest (furthest) neighbors are set NaN.   
%
%   See also idSocial_prepareTrajectories3D, 
%   idSocial_distanceDistribution

%   2014 Robert C. Hinz, Gonzalo G. de Polavieja, Consejo Superior de Investigaciones Científicas


if nargin<2 || isempty(bodylength)
    bodylength=1;
end



[trajectory,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(trajectory);


% Calculations
foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,no_dim);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf
            foc_to_nb_vec(:,ff,nf,:)=squeeze(trajectory(:,nf,ff,:)-trajectory(:,ff,ff,:));
        end
    end
end
distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,4))./bodylength;

distnce=permute(distance_focal_neighbour,[2,3,1]);
distnce(logical(repmat(eye(no_fish,no_fish),[1 1 no_frames])))=NaN;

nndist=NaN(no_fish,no_fish,no_frames);
fardist=NaN(no_fish,no_fish,no_frames);
for ff=1:no_fish
    act_dist=squeeze(distance_focal_neighbour(:,ff,:));
    act_dist(:,ff)=NaN;
    mn=min(act_dist,[],2);
    act_dist_mn=act_dist;
    ind = bsxfun(@gt, act_dist_mn, mn);
    act_dist_mn(ind)=NaN;
    nndist(ff,:,:)=act_dist_mn';
       
    mx=max(act_dist,[],2);
    act_dist_mx=act_dist;
    ind = bsxfun(@lt, act_dist_mx, mx);
    act_dist_mx(ind)=NaN;
    fardist(ff,:,:)=act_dist_mx';
end





