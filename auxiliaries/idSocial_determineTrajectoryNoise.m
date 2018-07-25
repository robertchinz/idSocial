function [mean_std mean_wgt_std maxdeviation]=idSocial_determineTrajectoryNoise(datosegm,area)

%% Apply smoothing spline
no_fish = datosegm.n_peces;

% I want to calculate the tolerance for the smoothing spline
% (tol>=E(f)=\sum_{j=1}_{n}(w(j)*|y(:,j)-f(x(j))|^2)), which measures the OVERALL
% distance between the original trajectory and the smoothed
% one. Therefore I want to convert the tolerance I would like
% to apply to every single data point (i.e.,the distance between
% trajectory(frame,focal,:) and tr_smooth(frame,focal,:) should
% be smaller than tol_SinglePoint) to the overall tolerance
% used for the smoothing spline.
% To do so, I assume that each measured coordinate is subject to
% noise (changes in light, noise in photo sensor), and that
% the noise is gaussian, centered at the measured coordinate
% and standard deviation given by the changes in blob size
% between frames (see below to see the exact measure).
% I apply this noise to the actual trajectory, and calculate
% the resulting sum of distances between original and noisy
% trajectory:

% Firstly, in order to get an idea of the 'amount of noise'
% I have a look at how a change of blob size influences the
% coordinate of center of mass
segm_dir=datosegm.directorio;
try
    segm=load_encrypt([segm_dir 'segm_1.mat'],1);
catch
    segm=load_encrypt([segm_dir 'segm_1.mat'],0);
end

no_frames_thresh=size(segm,2);

threshh = 0:.01:1;
use_frames = 1:10:no_frames_thresh;

areadens=NaN(no_fish,100);
areaedges=NaN(no_fish,100);
for ff = 1:no_fish
    [areadens(ff,:), areaedges(ff,:)]= ksdensity(area(:,ff));
end

minarea = min(area);
maxarea = max(area);

mean_coord=NaN(numel(use_frames),no_fish,2);
wgt_mean_coord=NaN(numel(use_frames),no_fish,2);
std_coord=NaN(numel(use_frames),no_fish,2);
no_frames4wgt_mean=NaN(numel(use_frames),no_fish);
wgt_std_coord=NaN(numel(use_frames),no_fish,2);
maxdev=NaN(numel(use_frames),no_fish);
minfr_array=cell(numel(use_frames),1);
fr_count=1;
for fr_act=use_frames
    
    miniframe=segm(fr_act).miniframes;
    coord = NaN(numel(threshh),no_fish,2);
    thr_count=ones(1,no_fish);
    area_check=NaN(numel(threshh),no_fish);
    for ff=1:no_fish
        if ff<=numel(miniframe) && ~isempty(miniframe{ff})
            try
                minfr_norm=double(miniframe{ff})./double(max(miniframe{ff}(:)));
            catch
                keyboard
            end
            
            for thr=threshh
                minfr_bw=minfr_norm<thr;
                minfr_array{thr_count(ff)}=minfr_bw;
                CC = bwconncomp(minfr_bw);
                RP = regionprops(CC,'Centroid','Area');
                
                if  ~isempty(RP) && CC.NumObjects==1
                    area_act=RP.Area;
                    if area_act > minarea(ff) && ...
                            area_act < maxarea(ff)
                        area_act=RP.Area;
                        coord(thr_count(ff),ff,:)=RP.Centroid;
                        area_check(thr_count(ff),ff) = area_act;
                        thr_count(ff)=thr_count(ff)+1;
                    end
                end
                %             if ~isempty(RP) && CC.NumObjects==1 % Only take new coordinate if blob has not split up into various blobs
                %                 area_act=RP.Area;
                %                 coord(thr_count(ff),ff,:)=RP.Centroid;
                %                 area_check(thr_count) = area_act;
                %                 thr_count(ff)=thr_count(ff)+1;
                %             end
            end
        end
        
    end
    mean_coord(fr_count,:,:)=nanmean(coord,1);
    try
    maxdev(fr_count,:)=mean(max(abs(coord-squeeze(repmat(mean_coord(fr_count,:,:),[size(coord,1) 1 1])))),3);
    catch
            maxdev(fr_count,:)=mean(max(abs(coord-repmat(mean_coord(fr_count,:,:),[size(coord,1) 1 1]))),3);

    end
        std_coord(fr_count,:,:)=nanstd(coord,[],1);
    
    
    %Calculate weights for weighted mean and std:
    for ff=1:no_fish
        % Look how often the area in question occurs 'in
        % real life':
        [~,inbin]=histc(area_check(:,ff),areaedges(ff,:));
        % For later use save number of 'cuts' which enter
        % into each mean value and standard deviation:
        no_frames4wgt_mean(fr_count,ff)=nansum(inbin>0);
        % Calculate the weighted mean and standard deviation from each cut:
        ttt=areadens(inbin(inbin>0));
        if ndims(ttt)==2 && size(ttt,1)==1;ttt=ttt'; end;
        wgt = repmat(ttt,[1 1 2]);
        
        wgt_mean_coord(fr_count,ff,:)=sum(squeeze(coord(inbin>0,ff,:).*wgt))/sum(inbin>0);
        for dim=1:size(coord,3)
            wgt_std_coord(fr_count,ff,dim)=sqrt(var(squeeze(coord(inbin>0,ff,dim)),wgt(:,:,dim),1));
        end
    end
    fr_count=fr_count+1;
end
mean_std = nanmean(nanmean(std_coord,3),1);
maxdeviation = max(maxdev,[],1); % Approximation of the max. deviation from the mean coordinate

try
mean_wgt_std = nanmean(squeeze(nansum(wgt_std_coord.*repmat(no_frames4wgt_mean,[1 1 2])))./repmat(nansum(no_frames4wgt_mean,1)',[1 2]),2)';
catch
    mean_wgt_std = nanmean(squeeze(nansum(wgt_std_coord.*repmat(no_frames4wgt_mean,[1 1 2])))./repmat(nansum(no_frames4wgt_mean,1)',[1 2])',2)';

end
