function [tr_smooth, tol_out] = idSocial_auxiliaries_smoothAcc2Vel2TrajectoryAdaptiveSpline(trajectories,tol_frame,filter_outliers)
if nargin<2 || isempty(tol_frame)
    tol_frame = 3;
end

if nargin<3 || isempty(filter_outliers)
    filter_outliers = true;
end

display_plot = false;
err = .05*tol_frame;
max_count = 50;
spl_deg = 2; % 2=cubic, 3=quintic
no_frames = size(trajectories,1);
no_fish = size(trajectories,2);


if ndims(trajectories)==4
    no_dims = size(trajectories,4);
    trtemp = NaN(no_frames,no_fish,no_dims);
    for ff=1:no_fish
        trtemp(:,ff,:)=trajectories(:,ff,ff,:);
    end
    trajectories = trtemp;
else
    no_dims = size(trajectories,3);
end


scrsz = get(0,'ScreenSize');

% We don't want to have NaNs when calculating acc. They mess with cumsum!
trNaN = isnan(trajectories);
% Get rid of 'fake frames' from interpolation, which through acc
% (diff(_,2,_)) have infected the 2 following frames as well!
trNaNext = trNaN | cat(1,trNaN(2:end,:,:),false(1,no_fish,no_dims)) | cat(1,trNaN(3:end,:,:),false(2,no_fish,no_dims));

trajectories=idSocial_interpolateTrajectories(trajectories,'linear');

vel = vertcat(diff(trajectories,1,1),NaN(1,no_fish,size(trajectories,3)));
acc = vertcat(diff(trajectories,2,1),NaN(2,no_fish,size(trajectories,3)));
tr_smooth = NaN(size(trajectories));



for ff = 1:no_fish
    if display_plot
        fh=figure('Name',['Smoothing Spline Information, Individual ' num2str(ff)],'Position',[scrsz(3)*.2 scrsz(4)*.2 scrsz(3)*.5 scrsz(4)*.5 ]);
    end
    tol = 1-8.2864e-07;%75000;
    tolPrev = 0;
    tolPrev2 = NaN;
    deltaPrev = NaN;
    delta = inf;
    tolDeltaArray = NaN(max_count+2,3);
    tolDeltaArray(1,:) = [tolPrev2 deltaPrev 1];
    tolDeltaArray(2,:) = [tolPrev delta 2];
    ttt=0;
    count=1;
    prctile99 = inf;
    prctile99Prev = inf;
    converged = false;
    dist2o = inf(no_frames,no_dims);
    
    while (all(isinf(dist2o(:))) || prctile(dist2o(:),99)-tol_frame > 0 || prctile(dist2o(:),99)-tol_frame< -err)  &&  ...
            count<max_count &&  ~converged
        % while (all(isinf(dist2o(:))) || max(dist2o(:))-tol_frame > 0 || max(dist2o(:))-tol_frame< -err)  &&  ...
        %         count<max_count &&  ~converged
        acc2 = NaN(no_frames,no_dims);
        
        
        for dm = 1:no_dims
            
            accnan = isnan(squeeze(acc(:,ff,dm)));
            idx = 1:size(trajectories,1); idx = idx(~accnan);
            %             ppaccx = spaps(idx,squeeze(acc(~accnan,ff,dm)),tol,spl_deg);
            %             acc2(~accnan,dm) = fnval(idx,ppaccx);
            ppaccx = csaps(idx,squeeze(acc(~accnan,ff,dm)),1-tol);
            acc2(~accnan,dm) = fnval(idx,ppaccx);
            %             acc2(~accnan,dm) = ppaccx;
            
            
        end
        
        
        
        % Integrate acc to vel
        acc4vel = NaN(no_frames,no_dims);
        acc4vel(2:end,:) = acc2(1:end-1,:); acc4vel(1,:) = zeros(1,no_dims);
        acc2vel = cumsum(acc4vel,1);
        v0 = -nanmean(acc2vel-squeeze(vel(:,ff,:)));
        acc2vel = acc2vel + repmat(v0,[no_frames 1]);
        
        % Integrate vel to tr
        vel4tr = NaN(no_frames,no_dims);
        vel4tr(2:end,:) = acc2vel(1:end-1,:); vel4tr(1,:) = zeros(1,no_dims);
        vel2tr = cumsum(vel4tr,1);
        tr0 = -nanmean(vel2tr-squeeze(trajectories(:,ff,:)));
        vel2tr = vel2tr + repmat(tr0,[no_frames 1]);
        
        % Get rid of 'fake frames' from interpolation, which through acc
        % (diff(_,2,_)) have infected the 2 following frames as well!
        acc2vel(squeeze(trNaNext(:,ff,:)))=NaN;
        vel2tr(squeeze(trNaNext(:,ff,:)))=NaN;
        
        % Distance measures to originial trajectory
        delta2o = sqrt(nansum((acc2vel-squeeze(acc(:,ff,:))).^2,2));
        dist2o = sqrt(nansum((vel2tr-squeeze(trajectories(:,ff,:))).^2,2));
        
        % Check if smoothed trajectory meets conditions
        deltaPrev2 = deltaPrev;
        deltaPrev = delta;
        prctile99Prev2 = prctile99Prev;
        prctile99Prev = prctile99;
        prctile99 = prctile(dist2o(:),99);% max(dist2o(:));%
        delta = prctile99-tol_frame;
        
        try
            tolDeltaArray(count+1,:) = [tol delta count];
        catch
            keyboard
        end
        
        if display_plot
            % Control plots
            figure(fh);
            subplot(4,1,1)
            ctrs = 0:(max([tol_frame*2 prctile99Prev2 prctile99Prev prctile99]))/100:max([tol_frame*2 prctile99Prev2 prctile99Prev prctile99]);
            hist(dist2o(:),ctrs); hold on;
            plot([tol_frame tol_frame],get(gca,'YLim'),'k');
            plot([prctile99 prctile99],get(gca,'YLim'),'r');
            plot([prctile99Prev prctile99Prev],get(gca,'YLim'),'--r');
            plot([prctile99Prev2 prctile99Prev2],get(gca,'YLim'),':r');
            plot([tol_frame-err tol_frame-err],get(gca,'YLim'),'g');
            hold off;
            set(gca,'XLim',[0 max([tol_frame*2 prctile99Prev2 prctile99Prev prctile99])])
            title({['tolfr ' num2str(tol_frame) ', Repetition #' num2str(count) ', ' num2str(ttt)]; ...
                ['Prctle = ' num2str(prctile99) ', Delta=' num2str(max(dist2o(:)))]});
            
            idces = 1:size(trajectories,1);
            subplot(4,1,2); plot(squeeze(acc(idces,ff,1)),'--'); hold on; %plot(squeeze(vel(idces,ff,2)),':b');
            plot(acc2(idces,1),'--g','LineWidth',2);hold off;
            subplot(4,1,3); plot(squeeze(vel(idces,ff,1)),'--'); hold on; %plot(squeeze(vel(idces,ff,2)),':b');
            plot(acc2vel(idces,1),'--g','LineWidth',2);hold off;
            subplot(4,1,4); plot(squeeze(trajectories(idces,ff,1)),'k','MarkerSize',3,'LineWidth',1); hold on; %plot(squeeze(vel(idces,ff,2)),':b');
            plot(squeeze(trajectories(idces,ff,2)),'k','MarkerSize',3,'LineWidth',1);
            plot(vel2tr(idces,1),'--b','LineWidth',1.1);
            plot(vel2tr(idces,2),'--g','LineWidth',1.1);hold off;
            
            %%%
        end
        
        
        tolPrev3=tolPrev2;
        tolPrev2=tolPrev;
        
        tolPrev=tol;
        
        
        %     disp(['Count ' num2str(count)])
        
        sarray = tolDeltaArray;
        sarray(:,3) = max(sarray(:,3))-sarray(:,3)+1;
        [~, sidx] = sortrows(sarray,[2 1 3]);
        tolDeltaArray = tolDeltaArray(sidx,:);
        
        lowLims = double(tolDeltaArray(:,2) < -err) ;%& ~isnan(tolDeltaArray(:,2));
        lowLims(isnan(tolDeltaArray(:,2))) = NaN;
        highLims = double(tolDeltaArray(:,2) > 0) ;%& ~isnan(tolDeltaArray(:,2));
        highLims(isnan(tolDeltaArray(:,2))) = NaN;
        
        low2highID = diff(lowLims) == -1;
        lowID = find(low2highID);
        highID = lowID + 1;
        if ~isempty(lowID)
            %         if count ==5   ; keyboard; end
            %         if isequal(lowhighpair,[tolDeltaArray(lowID,3) tolDeltaArray(highID,3)])
            %             lowID = lowID - 1;
            %             highID = highID + 1;
            %         end
            lowhighpair = [tolDeltaArray(lowID,3) tolDeltaArray(highID,3)];
            
            %         [tolDeltaArray(lowID,1) tolDeltaArray(highID,1)]
            %         [tolDeltaArray(lowID,2) tolDeltaArray(highID,2)]
            tol = (tolDeltaArray(lowID,1) + tolDeltaArray(highID,1))/2;
            
            
            ttt=1;
            
            
            while any(tol == tolDeltaArray(:,1))
                if abs(tolDeltaArray(lowID,2)) <= tolDeltaArray(lowID+1,2)  && ...
                        tol ~= tolDeltaArray(lowID,1)
                    tol = (tol+tolDeltaArray(lowID,1))/2;
                    %                 fprintf('%.10f\n',tol)
                    ttt=2;
                elseif abs(tolDeltaArray(lowID,2)) > tolDeltaArray(lowID+1,2)  && ...
                        tol ~= tolDeltaArray(lowID+1,1)
                    tol = (tol+tolDeltaArray(lowID+1,1))/2;
                    ttt=3;
                else
                    tol = tol*.99;
                    ttt=4;
                end
            end
%             fprintf('%.30f\n',tol)
            % If indices of high-low-pair have not changed
            
        elseif ~any(lowLims==1) && any(highLims==1)
            tol = tolDeltaArray(tolDeltaArray(:,1) == min(tolDeltaArray(:,1)),1)/2;
            tol = tol(1); % In case there is more than one candidate.
        elseif any(lowLims==1) && ~any(highLims==1)
            tol = tolDeltaArray(tolDeltaArray(:,1) == max(tolDeltaArray(:,1)),1)*2+1;
            tol = tol(1); % In case there is more than one candidate.
        end
        
        count = count + 1;
        
        
        
    end
    if display_plot
        close(fh)
    end
    tr_smooth(:,ff,:) = vel2tr;
    tol_out = tolPrev;
    if count==max_count
        disp('Max. number of repetitions reached.')
    end
    disp(['Final smoothing parameter is ' num2str(tolPrev)])
    disp([num2str(sum(dist2o(:)>tol_frame)) ' Frames exceed max. distance. Max. error: ' num2str(max(dist2o(:))) ' pxl.'])
    
    if filter_outliers
        tr_smooth(dist2o(:)>tol_frame,ff,:)=NaN;
        
        disp([num2str(sum(dist2o(:)>tol_frame)) ' exceeding frames filtered.'])
    end
end