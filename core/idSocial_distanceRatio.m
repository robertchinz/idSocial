function [distnceHist,inbins,funcinfo]=idSocial_distanceRatio(trajectory,social_radius,bodylength)
%   Inter-indivual distance from trajectory.
%
%   distnce = ...
%           IDSOCIAL_DISTANCERATIO(TR,SOCIAL_RADIUS) 
%   calculates the ratio of time, a focal spends within a distance of
%   1,2,..N-1 of its neighbors.
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
%   SOCIAL_RADIUS determines the maximum distance for which two individuals 
%   are considered to be interacting socially. the default is:
%   SOCIAL_RADIUS = 5
%   distnceRatio = ...
%           IDSOCIAL_DISTANCERATIO(TR,SOCIAL_RADIUS,BODYLENGTH) 
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

if nargin<2 || isempty(social_radius)
    social_radius=5;
end

if nargin<3 || isempty(bodylength)
    bodylength=1;
end


invertY = true;
[trajectory,~,~,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(trajectory,invertY);

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

distnceRatio = distnce < social_radius;

ctrs=0:no_fish-1;
no_bins=size(ctrs,2);
distnceHist=NaN(no_fish,no_fish,no_bins);
inbins=NaN(no_fish,no_fish,no_frames);
no_good_frames = NaN(no_fish,no_fish);
for ff=1:no_fish
    
    val=squeeze(nansum(distnceRatio(ff,:,:),2));
    
    no_good_frames = sum(~isnan(val));

    
    [hitemp,bins]=histc(val,[ctrs-.5 ctrs(end)+.5]);
    
    distnceHist(ff,ff,:) = hitemp(1:end-1);
    
    if ~all(bins==0)
        inbins(ff,ff,:)=val;%accumarray(bins(~isnan(val) & bins>0),val(~isnan(val)& bins>0),[no_bins 1],@(x) {x});
    end
end
distnceHist=mat2cell(distnceHist,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));

fstack=dbstack(3);
funcinfo.Function = mfilename;
funcinfo.callerFunction = fstack(1).name;
funcinfo.no_fish = no_fish;
funcinfo.XTick = ctrs;
funcinfo.XTickLabel = cellstr(num2str(ctrs'+1))';


