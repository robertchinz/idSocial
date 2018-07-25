function [smooth_tr] = idSocial_auxiliaries_smoothTrajectoryDouglasPeucker(trajectories,tol)

if nargin < 2 || isempty(tol)
    tol=10;
end

no_frames = size(trajectories,1);
no_fish = size(trajectories,2);
smooth_tr = NaN(size(trajectories));

vel = NaN(size(smooth_tr));
vel(1:no_frames-1,:,:)=diff(trajectories,1,1);
dir_vec_all = NaN(size(trajectories));
for ff=1:no_fish
    [pt,idx] = ...
        dpsimplify(squeeze(trajectories(:,ff,:)),tol);
    smooth_tr(idx,ff,:) = pt(~isnan(pt(:,1)),:);
    
    dirvec = NaN(size(idx,1),2);
    dv=squeeze(smooth_tr(~isnan(smooth_tr(:,ff,1)),ff,:));
    dirvec(1:end-1,:) = diff(dv,1,1);
    %     dir_vec(1:size(dv,1)-1,ff,:)=diff(dv,1,1);
    dirvec_magn=sqrt(sum(dirvec.^2,2));
    dirvec_norm=bsxfun(@rdivide,dirvec,dirvec_magn);
    
    
    vel_proj = NaN(size(trajectories));
    for id=1:size(idx,1)-1
%         if id==size(idx,1)-1 && ff ==2
%             keyboard
%         end
        is=idx(id);
        ie=idx(id+1);
        
        if ~all(isnan(trajectories(is+1:ie-1,ff,1))) % Leave out pieces with only Nans
            dir_vec_all(is:ie-1,ff,:) = squeeze(repmat(dirvec_norm(id,:),[ie-is,1]));
            proj=sum(dir_vec_all(is:ie-1,ff,:) .* vel(is:ie-1,ff,:),3);
            proj0=proj;
            proj0(proj0<0)=0;
            
            cproj =cumsum(proj);
            
            cproj0 = NaN(1,numel(proj));
            cproj0(1) = max(proj(1),0);
            for k=2:numel(proj)
                % Project only forward movements
                if proj(k) > 0 && cproj(k-1) + proj(k) >= max(cproj0)
                    pract = min(proj(k),max(cproj(1:k))-max(cproj0));
                elseif ~isnan(proj(k))
                    pract = 0;
                elseif isnan(proj(k))
                    pract=NaN;
                end
                
                cproj0(k) = cproj0(k-1) + pract;
            end
            cproj0(cproj0>cproj(end))=repmat(cproj(end),[sum(cproj0>cproj(end)),1]);
            
            
            new_trace = squeeze(repmat(smooth_tr(is,ff,:),[ie-is,1,1])) +  repmat(cproj0(1:end),[2 1])' .*squeeze(dir_vec_all(is:ie-1,ff,:));
            
            if  sqrt((new_trace(end,1)-smooth_tr(idx(id+1),ff,1)).^2+(new_trace(end,2)-smooth_tr(idx(id+1),ff,2)).^2) > 1e-5 && ...
                    ~all( isnan(new_trace(end,:)) & isnan(squeeze(smooth_tr(idx(id+1),ff,:))')) && ...
                    ~all(isnan(squeeze(smooth_tr(idx(id+1),ff,:))'))
                %
                warning(['Inaccuracy: New trace= ' num2str(new_trace(end,:)) ', End point= ' num2str(squeeze(smooth_tr(idx(id+1),ff,:))')   ]);
                
            end
            
            smooth_tr(is+1:ie-1,ff,:) =new_trace(1:end-1,:);
        end
    end
    
end
%%
% idces=30000:30100;
% figure; plot(trajectories(idces,ff,1),trajectories(idces,ff,2),'-+'); hold on; plot(smooth_tr(idces,ff,1),smooth_tr(idces,ff,2),'-*g');
% idx1=64944;
% plot(trajectories(idx1,ff,1),trajectories(idx1,ff,2),'-or');
% plot(smooth_tr(idx1,ff,1),smooth_tr(idx1,ff,2),'-ok');
% plot(smooth_tr(idx1-10:idx1-1,ff,1),smooth_tr(idx1-10:idx1-1,ff,2),'-*k');
% plot(smooth_tr(idx1+1:idx1+10,ff,1),smooth_tr(idx1+1:idx1+10,ff,2),'-+k');
%
% plot(trajectories(idx1-1,ff,1),trajectories(idx1-1,ff,2),'-*k');
% plot(trajectories(idx1+1,ff,1),trajectories(idx1+1,ff,2),'-+k');
% figure; plot(new_trace(:,1),new_trace(:,2),'*'); hold on; plot(smooth_tr(is:ie,ff,1),smooth_tr(is:ie,ff,2),'+g');
% figure; plot(cproj); hold on; plot(cproj0,'g'); %plot(cproj01,'r')
% ['Numerical inaccuracy? New trace= ' num2str(new_trace(end,:)) ', End point= ' num2str(squeeze(smooth_tr(idx(id+1),ff,:))')   ]
% % dir in each projected coordinate:
% ttt=squeeze(repmat(smooth_tr(ie,ff,:),[ie-is+1,1,1])-new_trace)./repmat(sqrt((squeeze(smooth_tr(ie,ff,1))'-new_trace(:,1)).^2+(squeeze(smooth_tr(ie,ff,2))'-new_trace(:,2)).^2),[1 2])