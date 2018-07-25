function [tr_smooth, tol] =  idSocial_adaptiveTrajectorySmoothing(trajectories,weights,noise_mean,max_deviation,smooth_method)

no_frames=size(trajectories,1);
no_fish=size(trajectories,2);
no_dim=size(trajectories,3);

if nargin >= 4 && ~isempty(max_deviation)
    if numel(max_deviation)==1
        max_deviation = [max_deviation max_deviation];
    end
end

if nargin < 2 || isempty(weights)
    weights = ones(no_frames,no_fish);
end
if nargin < 4 || isempty(max_deviation)
    max_deviation = [inf inf];
end
if nargin < 5 || isempty(smooth_method)
    smooth_method = 'smoothing_spline';
end

%%
trajectories=trajectories(:,:,:,1:2);
rand_ang=2*pi*(rand(no_frames,no_fish)-.5);
rand_ang = repmat(rand_ang,[1 1 no_fish]);
tr_noise=trajectories+repmat(noise_mean,[no_frames no_fish no_fish 2]).*repmat(abs(randn(no_frames,no_fish,no_fish)),[1 1 1 2]).* cat(4,cos(rand_ang), sin(rand_ang));
% cm=lines(128);

% tr_noise=NaN(size(trajectories));
% for ff=1:no_fish
%     tr_noise(:,ff,1)=trajectories(:,ff,1)+noise_mean(ff).*randn(no_frames,1);
%     tr_noise(:,ff,2)=trajectories(:,ff,2)+noise_mean(ff).*randn(no_frames,1);
% end
% The distribution of distances between tr_noise and trajectories does not have its maximum
% at zero, because there are more data points on a circle with a radius >0 (~0.5) than
% at cero.

% figure; hist(sqrt((tr_noise(:,1,1)-trajectories(:,1,1)).^2+(tr_noise(:,1,2)-trajectories(:,1,2)).^2),200)
% title('Distribution of distances between noisy and original trajectory')

% Calculating the resulting tolerance: Distance between
% noisy and smoothed trajectory divided by 2 since splines
% will be applied seperately to x and y coordinates.
% tol=nansum(sqrt((tr_noise(:,1,1)-trajectories(:,1,1)).^2+(tr_noise(:,1,2)-trajectories(:,1,2)).^2))/2;
tol = (nansum(sqrt((tr_noise(:,1,1,1)-trajectories(:,1,1,1)).^2)) + ...
    nansum(sqrt((tr_noise(:,1,1,2)-trajectories(:,1,1,2)).^2)))/2;

% smooth_degree = 31;

% There are two Matlab functions for smoothing splines:
% spaps, which lets you define the tolerance, and
% csaps, which asks for the 'smoothing parameter rho'.
% In my case, I am more interested in define the tolerance
% (which I can roughly determine from blob size, see above).
% Given the tolerance, spaps calculates the adequate rho.
% [NOTE: Tried the following, but with strange results:]
% Since I have to apply the smoothing to x- and y-coordinates
% seperately, and thus spaps calculates different rhos for
% both cases, I decided to use a combination of spaps and
% csaps: First, I use spaps to define the tolerance
% and calculate rho for each individual and dimension. Then
% for each individual I take the average rho over all dimensions,
% (here only x and y), and use the new, average rho to get
% the smoothed trajectory using csaps.
within_error = false;
count = 1;
hi=cell(1,1);
figure;
% disth=axes;

while ~within_error
    tr_smooth=NaN(size(trajectories));
    
    % Spline smoothing lets us chooses weights for each coordinate,
    % i.e., when calculating the tolerance ('sum over distances
    % between smoothed and original trajectories') it uses a weighted
    % sum.
    % I will use idTrackers 'probtrajectories', which reflects the
    % reliability of the assignment of a given coordinate, to define
    % weights.
    % In addition, probtrajectories contains values of -1 and -2,
    % referring to idTrackers internal interpolation procedure:
    % -1: There is a crossing of individuals, but blobs can be
    % recovered by erosion, and the center of mass coordinates
    % of those blobs are interpolated in order to get approximate
    % coordinates of the overlapping individuals during the crossing.
    % -2: Similar to -1, but erosion does not recover seperate
    % blobs. Interoplation is done under the only condition that
    % the interpolated trajectory coordinates lie within the
    % blob in a certain frame.
    % I will use probtrajectories as follows:
    % * 0 < Values <= 1 in probtrajectories -> weight = Value
    % * Value = 0 (even though I am not sure this happens) ->
    %               weight = min_weight.
    %   (When setting weight = 0, the coordinate simply disappears
    %   in the smoothed trajectories. By assigning it a small weight
    %   it will still somehow used, even though with a possibly
    %   unrealistic result. This can later be dealt with by
    %   filtering the trajectory using probtrajectory and a threshhold
    %   (as already implemented in idSocial)
    % * In order to treat the effective gaps between 'reliable'
    %   coordinates and gaps where interpolation hs been done
    %   (probtrajectory=-1 or -2), I increase the weight of the
    %   reliable coordinates at the start and the end of the gaps,
    %   forcing the smoothed trajectory to pass very closely
    %   by those reliable coordinates.
    %   In order to do so, I will multiply the weights/values
    %   of probtrajecry in those coordinates by
    %   'increase_border_weight' (this way, at least somehow the
    %   original weight is preseerved, and weights for values
    %   of probatrajectory==1 remain more important than those
    %   0 < values < 1.
    
    min_weight = 1e-1;
    increase_border_weight = 2;
    spline_weights=weights;
    
    probtr_lt_zero = spline_weights<0;
    gap_start_border = vertcat(diff(probtr_lt_zero,1,1), zeros(1,no_fish)) == +1;
    gap_end_border = vertcat(zeros(1,no_fish),diff(probtr_lt_zero,1,1)) == -1;
    gap_border = gap_start_border | gap_end_border;
    
    spline_weights(gap_border) = spline_weights(gap_border) * increase_border_weight;
    
    spline_weights(probtr_lt_zero) = min_weight;
    
    
    switch smooth_method
        
            
        case 'smoothing_spline_adaptive'
            rho=NaN(no_fish,1);
            for ff=1:no_fish
                
                
                [~,values,rho(ff)] = spaps(1:no_frames,squeeze(trajectories(:,ff,ff,:))',tol,...
                    squeeze(spline_weights(:,ff))');
                %         [values,rho(ff)] = csaps(1:no_frames,squeeze(trajectories(:,ff,:))',[],1:no_frames,...
                %             squeeze(spline_weights(:,ff))');
                tr_smooth(1:size(values,2),ff,ff,:)=values';
            end
            
            
            % Distance between original and smoothed trajectory
            dist_orig_smooth=NaN(no_frames,no_fish);
            for ff=1:no_fish
                dist_orig_smooth(:,ff)=squeeze(spline_weights(:,ff)).*...
                    sqrt((tr_smooth(:,ff,1)-trajectories(:,ff,1)).^2+(tr_smooth(:,ff,2)-trajectories(:,ff,2)).^2);
                nanidx=isnan(dist_orig_smooth(:,ff));
                if all(dist_orig_smooth(~nanidx,ff)<max_deviation(ff))
                    within_error=true;
                else
                    within_error=false;
                end
            end
            
            disp(['Tolerance = ' num2str(tol) ', Smoothing parameter p=' num2str(1-1/(1+rho(1)))])
            
            if within_error==false
                tol = tol*.9;
            end
            
            
        case {'moving'}
            for ff=1:no_fish
                for d=1:no_dim
                    tr_smooth(:,ff,d)=smooth(trajectories(:,ff,d),smooth_degree,smooth_method);
                end
            end
            
            % Distance between original and smoothed trajectory
            dist_orig_smooth=NaN(no_frames,no_fish);
            for ff=1:no_fish
                dist_orig_smooth(:,ff)=squeeze(spline_weights(:,ff)).*...
                    sqrt((tr_smooth(:,ff,1)-trajectories(:,ff,1)).^2+(tr_smooth(:,ff,2)-trajectories(:,ff,2)).^2);
                nanidx=isnan(dist_orig_smooth(:,ff));
                if all(dist_orig_smooth(~nanidx,ff)<max_deviation(ff))
                    within_error=true;
                else
                    within_error=false;
                end
            end
            disp(['Smoothing span = ' num2str(smooth_degree) ', max distance smooth to original = ' num2str(max(dist_orig_smooth(:)))])
            
            if within_error==false
                smooth_degree = smooth_degree - 2;
            end
            
            
    end
    
    % Filtering using weights. In principle this will
    % be done later if using idSocial, but should be done because
    % splines tend to act very strangely at sites with small weights/
    % small reliability.
    % probtr_filter = weights<.5;
    % for ff=1:no_fish
    %        tr_smooth(probtr_filter(:,ff),ff,:)=NaN;
    % end
    
    dist_orig_smooth_noweights=NaN(no_frames,no_fish);
    for ff=1:no_fish
        dist_orig_smooth_noweights(:,ff)=...
            sqrt((tr_smooth(:,ff,1)-trajectories(:,ff,1)).^2+(tr_smooth(:,ff,2)-trajectories(:,ff,2)).^2);
    end
    hi{count}=histc(dist_orig_smooth_noweights(:),0:.02:3);
    
    %     axes(disth);plot(0:.02:3,hi{count},'Color',cm(count,:))
    %     hold on
    %     drawnow;
    count=count+1;
end
% title(disth,'Distribution of distances between smoothed and original trajectory')

% % Filtering frames in which distance between original and
% % smoothed trajectory is greater than maxdeviation
% dist_orig_smooth_filter=dist_orig_smooth;
% dist_filter=false(size(dist_orig_smooth));
% for ff=1:no_fish
%
%     dist_filter(dist_orig_smooth_filter(:,ff)>max_deviation(ff),ff)=true;
%
%     tr_smooth(dist_filter(:,ff),ff,:)=NaN;
% end
