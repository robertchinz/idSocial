function [tr_temp, tol] = idSocial_auxiliaries_smoothTrajectories(tr,smooth_method,smooth_degree,max_deviation,noise_mean)

if nargin < 2 || isempty(smooth_method)
    smooth_method = 'moving';
end
if nargin < 3 || isempty(smooth_degree)
    smooth_degree = 30;
end
if nargin < 4 || isempty(max_deviation)
    max_deviation = NaN;
end
if nargin < 5 || isempty(noise_mean)
    noise_mean = .1;
end

no_frames = size(tr,1);
no_fish = size(tr,2);
if ndims(tr)==4
    no_dim = size(tr,4);

    trtemp = NaN(no_frames,no_fish,no_dim);
    for ff = 1:no_fish;
        trtemp(:,ff,:) = tr(:,ff,ff,:);
    end
else
    no_dim = size(tr,3);
end

% if no_fish==1 && nadims(trtemp)
%     tr_temp=reshape(tr,[no_frames no_fish no_dim]);
% else
%     tr_temp=squeeze(tr(:,:,1,:));
% end

tol = NaN;
switch smooth_method
    case {'moving', 'lowess','loess','sgolay','rlowess','rloess'}
        if strcmpi('smooth_method','sgolay')
            keyboard
        end
        tr_temp = NaN(size(tr));
        for ff=1:no_fish
            for d=1:no_dim
                tr_temp(:,ff,d)=smooth(tr(:,ff,d),smooth_degree,smooth_method);
            end
        end
    case 'adaptive_moving'
        [tr_temp, mavgwidth] = idSocial_auxiliaries_smoothTrajectoryAdaptiveMovAvg(tr,smooth_degree);
    case 'adaptive_moving_acc'
        %                 keyboard
        
        [tr_temp, mavgwidth] = idSocial_auxiliaries_smoothAcc2Vel2TrajectoryAdaptiveMovAvg(tr,smooth_degree);
    case 'smoothing_spline'
        %                 keyboard
        tr_temp = idSocial_auxiliaries_splineSmoothing(tr,max_deviation,[],smooth_degree);
    case 'smoothing_spline_acc'
        %                 keyboard
        tr_temp = idSocial_auxiliaries_smoothAcc2Vel2TrajectoryAdaptiveSpline(tr,max_deviation);
    case 'smoothing_spline_adaptive'
        [tr_temp, tol] = ...
            idSocial_adaptiveTrajectorySmoothing(tr,[],noise_mean,max_deviation,smooth_method);
        
    case {'DouglasPeucker','iterative_end_point','split_and_merge'}
        tr_temp = idSocial_auxiliaries_smoothTrajectoryDouglasPeucker(tr,smooth_degree);
    case 'moving_separate'
        tr_temp = idSocial_auxiliaries_smoothTrVelAccSep(tr,smooth_degree,'moving');
end