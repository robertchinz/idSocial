function tr_out = idSocial_auxiliaries_ProjectTrajectoryOnSpline(trajectories,tr_spline)

% nntol = 5*tol;

no_frames=size(trajectories,1);
no_fish=size(trajectories,2);

delta = sqrt(nansum(diff(tr_spline,1,1).^2,4));
dst2o = NaN(no_frames,no_fish);
for ff=1:no_fish
    dst2o(:,ff) = sqrt(sum((squeeze(trajectories(:,ff,:))-squeeze(tr_spline(:,ff,1:2))).^2,2));
end



idces=1:no_frames;
%%
tr_out = NaN(size(trajectories));
select_idx = NaN(size(trajectories,1),size(trajectories,2));
maxcount = no_frames;
tic
ipf = ceil(max(delta(:)));
disp(['Interpolation factor: ' num2str(ipf)]);
for ff = 1:no_fish
    
    
    fprintf('Focal %d: ',ff)
    tr0X = squeeze(trajectories(idces,ff,1));
    tr0Y = squeeze(trajectories(idces,ff,2));
    trX = squeeze(tr_spline(idces,ff,1));
    trY = squeeze(tr_spline(idces,ff,2));
    
    trXip = NaN(size(trajectories,1)*ipf-(ipf-1),1);
    trYip = NaN(size(trajectories,1)*ipf-(ipf-1),1);
    trXip((idces(1)-1)*ipf+1:(idces(end)-1)*ipf+1) = interp1(idces,trX,idces(1):1/ipf:idces(end));
    trYip((idces(1)-1)*ipf+1:(idces(end)-1)*ipf+1) = interp1(idces,trY,idces(1):1/ipf:idces(end));
    
    for fr = idces
        nntol = dst2o(fr,ff)*1.1;
        if mod(fr,100)==0; fprintf('%d,',fr); end
        if mod(fr,2000)==0; fprintf('\n'); end
        pt = [tr0X(fr-idces(1)+1) tr0Y(fr-idces(1)+1)];
        % List of possible candidates (prev):
        dst = -inf;
        cnt = 1;
        while (dst < nntol || isnan(dst)) && ((fr-1)*ipf)+1 - cnt+1 > 0 && cnt<maxcount
            if ~isnan(dst)
                
                act_trXip =[trXip(((fr-1)*ipf)+1 - cnt+1) trYip(((fr-1)*ipf)+1-cnt+1)]; dst = sqrt(   sum((pt - act_trXip).^2));
                
                cnt = cnt+1;
            else
                NotNaNidx = find(~isnan(trXip(((fr-1)*ipf)+1 - cnt+1:-1:1)) & ...
                    ~isnan(trYip(((fr-1)*ipf)+1 - cnt+1:-1:1)),1,'first');
                if ~isempty(NotNaNidx)
                    dst = sqrt(   sum((pt - [trXip(((fr-1)*ipf)+1 - cnt + 1 - NotNaNidx + 1) trYip(((fr-1)*ipf)+1-cnt+1 - NotNaNidx + 1)]).^2));
                    cnt = cnt + NotNaNidx;
                else % There is no more not-NaN coming. Exit.
                    dst = inf;
                end
            end
        end
        if cnt>=maxcount; disp('MAX!'); end
        prev_cnd = cnt-2;
        % List of possible candidates (after):
        dst = -inf;
        cnt = 1;
        while dst < nntol && ((fr-1)*ipf)+1 + cnt < numel(trXip)  && cnt<maxcount
            if ~isnan(dst)
                dst = sqrt(sum((pt - [trXip(((fr-1)*ipf)+1 + cnt) trYip(((fr-1)*ipf)+1 + cnt)]).^2));
                cnt = cnt+1;
            else
                NotNaNidx = find(~isnan(trXip(((fr-1)*ipf)+1 + cnt - cnt+1:1:end)) & ...
                    ~isnan(trYip(((fr-1)*ipf)+1 + cnt - cnt+1:1:end)),1,'first');
                if ~isempty(NotNaNidx)
                    dst = sqrt(   sum((pt - [trXip(((fr-1)*ipf)+1 + cnt + NotNaNidx) trYip(((fr-1)*ipf)+1 + cnt + NotNaNidx)]).^2));
                    cnt = cnt + NotNaNidx;
                else % There is no more not-NaN coming. Exit.
                    dst = inf;
                end
            end
        end
        if cnt>=maxcount; disp('MAX!'); end
        after_cnd = cnt-2;
        cnd_idces = ((fr-1)*ipf)+1 - prev_cnd+1 : ((fr-1)*ipf)+1 + after_cnd;
        
        if ~isempty(cnd_idces)
            dst = sqrt(sum((repmat(pt,[size(cnd_idces,2),1]) - [trXip(cnd_idces) trYip(cnd_idces)]).^2,2));
            [~, nn_idx] = min(dst);
            nn_idx = cnd_idces(nn_idx);
            % Force minimum movement
            %             select_idx(fr,ff) = max(nn_idx,select_idx(max(fr-1,1),ff)+1);
            select_idx(fr,ff) = nn_idx;
            %         tr1(fr,1,:) = [trXip(nn_idx) trYip(nn_idx)];
        else
            if fr>1 && ~any(isnan(pt)); disp('No candidates found!'); end;
            %             keyboard
            select_idx(fr,ff) = ((fr-1)*ipf)+1;
            %         tr1(fr,1,:) = [trXip(((fr-1)*ipf)+1) trYip(((fr-1)*ipf)+1)];
        end
    end
    
    % Sort idces: Allow no backward movements, guarantee that
    % direction of movement always in the same direction even when there is hardly any movement at all (noise!)
    select_idx(idces,ff) = sort(select_idx(idces,ff),1);
   
    % Smooth the smooth
%     pp=csaps(1:size(select_idx(idces,ff),1),select_idx(idces,ff),.2); 
%     select_idx(idces,ff) = round(fnval(pp,1:size(select_idx(idces,ff),1)));
    select_idx(idces,ff) = round(smooth(select_idx(idces,ff),5));
    tr_out(idces,ff,:) = [trXip(select_idx(idces,ff)) trYip(select_idx(idces,ff))];
    fprintf('\n')
end
toc