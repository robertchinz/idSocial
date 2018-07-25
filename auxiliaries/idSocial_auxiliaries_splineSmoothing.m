function [tr_smooth, sp, vel, acc, tol_out, max_dist, d2oAll, no_frames_off] = idSocial_auxiliaries_splineSmoothing(trajectories,tol_frame,separate_smoothing_vel_acc,spl_deg)
filter_outliers =true;
display_plot = false;
if nargin<2 || isempty(tol_frame)
    tol_frame = 3;
end
% if nargin<3 || isempty(bodylength)
%     bodylength = inf;
% end
if nargin<3 || isempty(separate_smoothing_vel_acc)
    separate_smoothing_vel_acc = false;
end
if nargin<4 || isempty(spl_deg)
spl_deg = 2; % 2=cubic, 3=quintic
end

weight_increment = .2;
err = .1*tol_frame;
max_count = 30;

len = size(trajectories,1);
no_fish = size(trajectories,2);
dim = size(trajectories,4);

if ndims(trajectories)==4
    trtemp = NaN(len,no_fish,dim);
    for ff=1:no_fish
        trtemp(:,ff,:)=trajectories(:,ff,ff,:);
    end
    trajectories = trtemp;
end



% tol = len*tol_frame;% * tol_frame;%.^2;

dst_vel = NaN(len,size(trajectories,2));
dst_velx = NaN(len,size(trajectories,2));
dst_vely = NaN(len,size(trajectories,2));
nan_velx = false(len,size(trajectories,2));
nan_vely = false(len,size(trajectories,2));
dst_velz = NaN(len,size(trajectories,2));
d2oAll = NaN(len,size(trajectories,2),3);
dst2orig = NaN(len,size(trajectories,2),3);
zero_idces = NaN(len,size(trajectories,2),1);
zero_idcesx = false(len,size(trajectories,2),1); % obsolete
zero_idcesy = false(len,size(trajectories,2),1); % obsolete    
zero_idcesz = NaN(len,size(trajectories,2),1);
tr_smooth = NaN(size(trajectories));
vel = NaN(size(trajectories));
acc = NaN(size(trajectories));
sp=cell(no_fish,1);
scrsz = get(0,'ScreenSize');
if display_plot
    fh=figure('Name','Smoothing Spline Information','Position',[scrsz(3)*.2 scrsz(4)*.2 scrsz(3)*.5 scrsz(4)*.3 ]);
end

deltaX = inf;
deltaY = inf;
deltaZ = inf;
tol_out=NaN(no_fish,3);
max_dist=NaN(no_fish,3);
no_frames_off = NaN(no_fish,3); 
for ff=1:no_fish
    
%     dst_vel = cat(1,trajectories(2:end,:,:) - trajectories(1:end-1,:,:),NaN(1,no_fish,2));
%     dst_vel_norm =sqrt(nansum(dst_vel.^2,3));
%     dst_vel_norm(isnan(dst_vel_norm(:,ff)),ff)=0;
%     dst_vel_norm(:,ff) = vertcat(0,dst_vel_norm(2:end,ff));
%     zero_idces(:,ff) = dst_vel_norm(:,ff)==0;
    
    w = ones(1,len);
    wx = ones(1,len);
    wy = ones(1,len);
    wz = ones(1,len);
    d2ox = inf(len,1);
    d2oy = inf(len,1);
    d2oz = inf(len,1);
    tol = sum((tol_frame*(rand(1, len)-.5)).^2);
    tol = sum((tol_frame*(rand(1, len))).^2);
    tol = (0.00051*len*tol_frame).^2;
    
    tolx_prev = 0;
    toly_prev = 0;
    tolx=tol;
    toly=tol;
    count=1;
    %     while (~all(d2ox(~isnan(d2ox))<tol_frame) || prctile(d2ox,99)-tol_frame > 0 || prctile(d2ox,99)-tol_frame< -err)  && count<max_count
    while (all(isinf(d2ox)) || prctile(d2ox,99)-tol_frame > 0 || prctile(d2ox,99)-tol_frame< -err)  && count<max_count
        
        dst_velx(2:end,ff) = abs(trajectories(1:end-1,ff,1)-trajectories(2:end,ff,1));
        nan_velx(:,ff) = isnan(dst_velx(:,ff));
        dst_velx(nan_velx(:,ff),ff)=0;
        dst_velx(:,ff) = vertcat(0,dst_velx(2:end,ff));
        zero_idcesx(:,ff) = dst_velx(:,ff)==0;
        
%         keyboard
        xx=1:len; xx=xx(~zero_idcesx(:,ff));
        %         xx=cumsum(dst_velx(:,ff)); xx=xx(~zero_idcesx(:,ff));
%                 xx=cumsum(dst_vel_norm(:,ff)); xx=xx(~zero_idces(:,ff));
        
                [spX, trX]= ...
                    spaps(xx,trajectories(~zero_idcesx(:,ff),ff,1),tolx,wx(~zero_idcesx(:,ff)),spl_deg);
        
        tr_smooth(~zero_idcesx(:,ff),ff,1) = trX;
%         [spX, trX]= ...
%             spaps(xx,trajectories(~zero_idces(:,ff),ff,1),tolx,wx(~zero_idces(:,ff)),spl_deg);
%         tr_smooth(~zero_idces(:,ff),ff,1) = trX;
        deltaXprev = deltaX;
        d2ox = abs(squeeze(trajectories(:,ff,1))-squeeze(tr_smooth(:,ff,1)));
        deltaX = max(d2ox)-tol_frame;
        deltaX = prctile(d2ox,99)-tol_frame;
        
        %         wx(d2ox' > prctile(d2ox,99)) = 1000000;
        %         wx(d2ox' < prctile(d2ox,99)) = 1;
        %         wx = wx + ...
        %             sign((d2ox' - tol_frame))*weight_increment;
        %         wx = max(wx,.1 * ones(1,len));
        %         wx = wx/sum(wx)*len;
        
        tolx_prev2=tolx_prev;
        tolx_prev=tolx;
        
        if deltaX > 0 || deltaX < - err
            if deltaX < - err && deltaXprev > 0
                tolx = (tolx_prev2+tolx_prev)/2;
            elseif deltaX > 0 && deltaXprev < - err
                tolx = (tolx_prev2+tolx_prev)/2;
            elseif deltaX < - err && deltaXprev < - err
                tolx = tolx_prev*2;
            elseif deltaX > 0  && deltaXprev > 0
                tolx = tolx_prev/2;%nansum(d2ox(d2ox<tol_frame)) + (nansum(d2ox>tol_frame)*tol_frame/2) ;%nanmedian(d2o)*len;
                
            end
        end
        count = count + 1;
        
        if display_plot
            figure(fh);
%             subplot(1,2,1);
            hist(d2ox,200); hold on;
%             plot([tol_frame tol_frame],get(gca,'YLim'),'k','LineWidth',2);
            plot([prctile(d2ox,99) prctile(d2ox,99)],get(gca,'YLim'),'r','LineWidth',2);
            plot([max(d2ox) max(d2ox)],get(gca,'YLim'),'k','LineWidth',2);
            yl = get(gca,'YLim');
            x1 = [tol_frame-err tol_frame tol_frame tol_frame-err];
            y1 = [0 0 yl(2) yl(2)];
            patch('XData',x1, ...
                'YData',y1, ...
                'FaceColor',[0 1 0],'FaceAlpha',.5,'EdgeColor','none');
%             plot([tol_frame-err tol_frame-err],get(gca,'YLim'),'g','LineWidth',2);
            hold off;
            set(gca,'XLim',[0 tol_frame*2])
%             set(gca,'XLim',[0 max(d2ox)+5])
            
            set(gcf,'Name',['X: Individual ' num2str(ff) ' (tolfr ' num2str(tol_frame) '), Repetition #' num2str(count) '. Delta=' num2str(max(d2ox)) 'Rel. Error: ' num2str(sum((d2ox'.*wx).^2)/(tol_frame*len))]);
%             title({['X: Individual ' num2str(ff) ' (tolfr ' num2str(tol_frame) '), Repetition #' num2str(count)];['Delta=' num2str(max(d2ox)) 'Rel. Error: ' num2str(sum((d2ox'.*wx).^2)/(tol_frame*len))]});
            xlabel('Distance [pxl]')
            ylabel('No. of Frames')
            
            
            text(tol_frame-err/2,yl(2)/2,'Goal','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',20)
            text(prctile(d2ox,99),yl(2)/2,'99th Percentile','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20)
            text(max(d2ox),yl(2)/2,'Max. Distance','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20)
            %             subplot(1,2,2); plot(wx); title('Weights x')
            %         if max(d2ox)<tol_frame; disp(num2str(max(d2ox))); end
            idSocial_improveGraphics(fh)
%             export_fig('E:\usuarios\Robert\Documents\Syncplicity Folders\Onset3\figures\supplement\Smoothing\SplineAdaptiveParameter.png')
%             keyboard
        end
    end
    for fr = 2:size(trajectories,1)
        if zero_idcesx(fr,ff) && ~nan_velx(fr,ff)
            tr_smooth(fr,ff,1) = tr_smooth(fr-1,ff,1);
        end
    end
    if count==max_count
        disp('Max. number of repetitions reached for X.')
    end
    disp(['X: ' num2str(sum(d2ox>tol_frame)) ' Frames exceed max. distance. Max. distance: ' num2str(max(d2ox))])
    count=1;
    while (all(isinf(d2oy)) || prctile(d2oy,99)-tol_frame > 0 || prctile(d2oy,99)-tol_frame< -err)  && count<max_count
        
        %     while (~all(d2oy(~isnan(d2oy))<tol_frame) || prctile(d2oy,99)-tol_frame > 0 || prctile(d2oy,99)-tol_frame< -err)  && count<max_count
        dst_vely(2:end,ff) = abs(trajectories(1:end-1,ff,2)-trajectories(2:end,ff,2));
        nan_vely(:,ff) = isnan(dst_vely(:,ff));
        dst_vely(nan_vely(:,ff),ff)=0;
        dst_vely(:,ff) = vertcat(0,dst_vely(2:end,ff));
        zero_idcesy(:,ff) = dst_vely(:,ff)==0;
        yy=1:len;yy=yy(~zero_idcesy(:,ff));
        %         yy=cumsum(dst_vely(:,ff)); yy=yy(~zero_idcesy(:,ff));
%         yy=cumsum(dst_vel_norm(:,ff)); yy=yy(~zero_idces(:,ff));
        try
                [spY, trY]=...
                    spaps(yy,trajectories(~zero_idcesy(:,ff),ff,2),toly,wy(~zero_idcesy(:,ff)),spl_deg);
        catch
            keyboard
        end
        tr_smooth(~zero_idcesy(:,ff),ff,2) = trY;
%         [spY, trY]=...
%             spaps(yy,trajectories(~zero_idces(:,ff),ff,2),toly,wy(~zero_idces(:,ff)),spl_deg);
%         tr_smooth(~zero_idces(:,ff),ff,2) = trY;
        
        deltaYprev = deltaY;
        d2oy = abs(squeeze(trajectories(:,ff,2))-squeeze(tr_smooth(:,ff,2)));
        deltaY = max(d2oy)-tol_frame;
        deltaY = prctile(d2oy,99)-tol_frame;
        
        %         wy = wy + ...
        %             sign((d2oy' - tol_frame))*weight_increment;
        %         wy = max(wy,.1 * ones(1,len));
        %         wy = wy/sum(wy)*len;
        %         toly = nansum(d2oy(d2oy<tol_frame)) + (nansum(d2oy>tol_frame)*tol_frame/2) ;%nanmedian(d2o)*len;
        
        toly_prev2=toly_prev;
        toly_prev=toly;
        if deltaY > 0 || deltaY < - err
            if deltaY < - err && deltaYprev > 0
                toly = (toly_prev2+toly_prev)/2;
            elseif deltaY > 0 && deltaYprev < - err
                toly = (toly_prev2+toly_prev)/2;
            elseif deltaY < - err && deltaYprev < - err
                toly = toly_prev*2;
            elseif deltaY > 0  && deltaYprev > 0
                toly = toly_prev/2;%nansum(d2ox(d2ox<tol_frame)) + (nansum(d2ox>tol_frame)*tol_frame/2) ;%nanmedian(d2o)*len;
                
            end
        end
        
        count = count + 1;
        if display_plot
            figure(fh);
            subplot(1,2,1);hist(d2oy,200); hold on; plot([tol_frame tol_frame],get(gca,'YLim'),'k');
            plot([prctile(d2oy,99) prctile(d2oy,99)],get(gca,'YLim'),'r');
            plot([tol_frame-err tol_frame-err],get(gca,'YLim'),'g');
            hold off;
            set(gca,'XLim',[0 tol_frame*1.3])
            title({['Y: Individual ' num2str(ff) '(tolfr ' num2str(tol_frame) '), Repetition #' num2str(count)];['Delta=' num2str(max(d2oy)) 'Rel. Error: ' num2str(sum((d2oy'.*wy).^2)/(tol_frame*len))]});
            subplot(1,2,2); plot(wy); title('Weights y')
        end
    end
    for fr = 2:size(trajectories,1)
        if zero_idcesy(fr,ff) && ~nan_vely(fr,ff)
            tr_smooth(fr,ff,2) = tr_smooth(fr-1,ff,2);
        end
    end
    if count==max_count
        disp('Max. number of repetitions reached for Y.')
    end
    disp(['Y: ' num2str(sum(d2oy>tol_frame)) ' Frames exceed max. distance. Max. distance: ' num2str(max(d2oy))])
    
%         dst_velz(2:end,ff) = abs(trajectories(1:end-1,ff,2)-trajectories(2:end,ff,2));
%         dst_velz(isnan(dst_velz(:,ff)),ff)=0;
%         dst_velz(:,ff) = vertcat(0,dst_velz(2:end,ff));
%         zero_idcesz(:,ff) = dst_velz(:,ff)==0;
%     count=1;
%     tolz = NaN;
%     if dim == 3 && ~all(trajectories(~zero_idcesz(:,ff),ff,3)==0)
%         while ~all(d2oz(~isnan(d2oz))<tol_frame) && count<max_count
%             %         zz=cumsum(dst_vel(~zero_idcesz(:,ff),ff));
%             zz=1:len; zz=zz(~zero_idcesz(:,ff));
%             [~, trZ]=...
%                 spaps(zz,trajectories(~zero_idcesz(:,ff),ff,3),tolz,wz(~zero_idcesz(:,ff)),spl_deg);
%             tr_smooth(~zero_idcesz(:,ff),ff,3) = trZ;
%             d2oz = abs(squeeze(trajectories(:,ff,3))-squeeze(tr_smooth(:,ff,3)));
%             wz = wz + ...
%                 sign((d2ox' - tol_frame))*weight_increment;
%             wz = max(wz,.1 * ones(1,len));
%             wz = wz/sum(wz)*len;
%             
%             tolz = nansum(d2oz(d2oz<tol_frame)) + (nansum(d2oz>tol_frame)*tol_frame/2) ;%nanmedian(d2o)*len;
%             count = count + 1;
%             figure(fh);
%             subplot(1,2,1);hist(d2oz,200); hold on; plot([tol_frame tol_frame],get(gca,'YLim'),'k'); plot([max(d2oz) max(d2oz)],get(gca,'YLim'),'r'); hold off;
%             set(gca,'XLim',[0 tol_frame*1.3])
%             title({['Z: Individual ' num2str(ff) ', Repetition #' num2str(count)];['Delta=' num2str(max(d2oz)) 'Rel. Error: ' num2str(sum((d2oz'.*w).^2)/(tol_frame*len))]});
%             subplot(1,2,2); plot(wz); title('Weights x')
%         end
%     end
    tolz = NaN;
    tol_out(ff,:) = [tolx toly tolz];
    max_dist(ff,:) = [max(d2ox) max(d2oy) max(d2oz)];
    no_frames_off(ff,:) = [sum(d2ox>tol_frame) sum(d2oy>tol_frame) sum(d2oz>tol_frame)];
    d2oAll(:,ff,1) = d2ox;
    d2oAll(:,ff,2) = d2oy;
    d2oAll(:,ff,3) = d2oz;
    
    
    sp{ff}.X=spX;
    sp{ff}.Y=spY;
    sp{ff}.indices_xx = ~zero_idcesx(:,ff);
    sp{ff}.indices_yy = ~zero_idcesy(:,ff);
%     sp{ff}.indices_xx = ~zero_idces(:,ff);
%     sp{ff}.indices_yy = ~zero_idces(:,ff);
    sp{ff}.xx = xx;%cumsum(dst_velx(~zero_idcesx(:,ff),ff));
    sp{ff}.yy = yy;%cumsum(dst_vely(~zero_idcesy(:,ff),ff));
    
   
    if ~separate_smoothing_vel_acc
        
        spvX = fnder(spX,1);
        spvY = fnder(spY,1);
        vxtemp = fnval(spvX,sp{ff}.xx);
        vel(sp{ff}.indices_xx,ff,1) = vxtemp;
        vytemp = fnval(spvY,sp{ff}.yy);
        vel(sp{ff}.indices_yy,ff,2) = vytemp;
        
        spaX = fnder(spX,2);
        spaY = fnder(spY,2);
        axtemp = fnval(spaX,sp{ff}.xx);
        acc(sp{ff}.indices_xx,ff,1) = axtemp;
        aytemp = fnval(spaY,sp{ff}.yy);
        acc(sp{ff}.indices_yy,ff,2) = aytemp;
    else
        [~, vel, acc] = idSocial_auxiliaries_smoothVelAndAcc(trajectories);
%         vel = vertcat(diff(trajectories,1,1),NaN(1,no_fish,size(trajectories,3)));
%         acc = vertcat(diff(trajectories,2,1),NaN(2,no_fish,size(trajectories,3)));
%         
%         disp('Smoothing parameter is so random!')
%         
%         velx = NaN(size(trajectories,1),1);
%         vely = NaN(size(trajectories,1),1);
%         velnanx = isnan(vel(:,ff,1));
%         velnany = isnan(vel(:,ff,2));
%         idxx = 1:size(trajectories,1); idxx = idxx(~velnanx);
%         idxy = 1:size(trajectories,1); idxy = idxy(~velnany);
%         ppvelx = csaps(idxx,squeeze(vel(~velnanx,ff,1)),.5);
%         ppvely = csaps(idxy,squeeze(vel(~velnany,ff,2)),.5);
%         velx(idxx) = fnval(idxx,ppvelx);
%         vely(idxy) = fnval(idxy,ppvely);
%         
%         accx = NaN(size(trajectories,1),1);
%         accy = NaN(size(trajectories,1),1);
%         accnanx = isnan(acc(:,ff,1));
%         accnany = isnan(acc(:,ff,2));
%         idxx = 1:size(trajectories,1); idxx = idxx(~accnanx);
%         idxy = 1:size(trajectories,1); idxy = idxy(~accnany);
%         ppaccx = csaps(idxx,squeeze(acc(~accnanx,ff,1)),.5);
%         ppaccy = csaps(idxy,squeeze(acc(~accnany,ff,2)),.5);
%         accx(idxx) = fnval(idxx,ppaccx);
%         accy(idxy) = fnval(idxy,ppaccy);
%         
% %         ppaccx = csaps(1:size(trajectories,1),squeeze(acc(:,ff,1)),.5);
% %         ppaccy = csaps(1:size(trajectories,1),squeeze(acc(:,ff,2)),.5);
% %         accx = fnval(1:size(trajectories,1),ppaccx);
% %         accy = fnval(1:size(trajectories,1),ppaccy);
%         
%         vel(1:size(trajectories,1),ff,1) = velx;
%         vel(1:size(trajectories,1),ff,2) = vely;
%         acc(1:size(trajectories,1),ff,1) = accx;
%         acc(1:size(trajectories,1),ff,2) = accy;
        
        sp{ff}.velx = squeeze(vel(:,ff,1));
        sp{ff}.vely = squeeze(vel(:,ff,2));
        sp{ff}.accx = squeeze(acc(:,ff,1));
        sp{ff}.accy = squeeze(acc(:,ff,2));
        %% acc
%                 idces = 200:300;
%                 acctol = 100000;
%                 ff = 1;
%         figure; plot(squeeze(acc(idces,1,1)),'--'); hold on; plot(squeeze(acc(idces,1,2)),':b');
%         
%         [~,accx] = spaps(1:size(trajectories,1),squeeze(acc(:,ff,1))',acctol);
%         [~,accy] = spaps(1:size(trajectories,1),squeeze(acc(:,ff,2))',acctol);
%         
% moving
%                 accx = smooth(squeeze(acc(:,ff,1))',10,'moving');
%                 accy = smooth(squeeze(acc(:,ff,2))',10,'moving');
%         
%         plot(accx(idces),'--g','LineWidth',2)
%         plot(accy(idces),':g','LineWidth',2)
%         plot(sqrt(nansum(acc(idces,ff,:).^2,3)),'k')
%          plot(sqrt(accx(idces).^2 + accy(idces).^2),'g')
         
         %% vel
%         idces = 200:300;
%         veltol = 5;%100000;
%         ff = 1;
%         figure; plot(squeeze(vel(idces,1,1)),'--'); hold on; plot(squeeze(vel(idces,1,2)),':b');
%         
% %         [~,velx] = spaps(1:size(trajectories,1),squeeze(vel(:,ff,1))',veltol);
% %         [~,vely] = spaps(1:size(trajectories,1),squeeze(vel(:,ff,2))',veltol);
%         
%        
% %         velx = smooth(squeeze(vel(:,ff,1))',veltol,'moving');
% %         vely = smooth(squeeze(vel(:,ff,2))',veltol,'moving');
%         
%         plot(velx(idces),'--g','LineWidth',2)
%         plot(vely(idces),':g','LineWidth',2)
%         plot(sqrt(nansum(vel(idces,ff,:).^2,3)),'k')
%         plot(sqrt(accx(idces).^2 + accy(idces).^2),'g')
        %% Try to reconstruct trajectory from smooth acc
%         idces = 1:size(trajectories,1);
%         accx4velx = NaN(size(trajectories,1),1);
%         accx4velx(2:end) = accx(2:end); accx4velx(1) = 6.36999;%nanmean(cumsum(accx));
%         accx2velx = cumsum(accx4velx);
%         figure; plot(squeeze(vel(idces,1,1)),'--'); hold on; plot(squeeze(vel(idces,1,2)),':b');
%          plot(accx2velx(idces),'--g','LineWidth',2)
%         
%         accy4vely = NaN(size(trajectories,1),1);
%         accy4vely(2:end) = accy(2:end); accy4vely(1) = -3;%nanmean(cumsum(accx));
%         accy2vely = cumsum(accy4vely);
%          figure; plot(squeeze(vel(idces,1,2)),'--'); hold on; plot(squeeze(vel(idces,1,2)),':b');
%          plot(accy2vely(idces),'--g','LineWidth',2)
%          
%          velx4trx = NaN(size(trajectories,1),1);
%          velx4trx(2:end) = accx2velx(2:end); velx4trx(1) = trajectories(1,1,1);%nanmean(cumsum(accx));
%          velx2trx = cumsum(velx4trx);
%          figure; plot(squeeze(trajectories(idces,1,1)),'--'); hold on; plot(squeeze(trajectories(idces,1,2)),':b');
%          plot(velx2trx(idces),'--g','LineWidth',2)
%          
%          vely4try = NaN(size(trajectories,1),1);
%          vely4try(2:end) = accy2vely(2:end); vely4try(1) = trajectories(1,1,1);%nanmean(cumsum(accy));
%          vely2try = cumsum(vely4try);
%          figure; plot(squeeze(trajectories(idces,1,1)),'--'); hold on; plot(squeeze(trajectories(idces,1,2)),':b');
%          plot(vely2try(idces),'--g','LineWidth',2)
        
    end
    
    dst2orig(:,ff,:) = [d2ox d2oy d2oz];
    
end
if display_plot
    close(fh)
end
% find(diff(find(exceeds_tol))>1);
% exceeds_tol = dst2orig > tol_frame;
% max_deviation = max(dst2orig);
% if filter_outliers
%     exceeds_tol2 = dst2orig(:,ff) > tol_frame*2;
%     
%     for ff=1:no_fish
%         tr_smooth(exceeds_tol2,ff,:)=NaN;
%     end
%     disp([mfilename ': Filtered ' num2str(sum(exceeds_tol2)) ' Frames after smoothing which exceeded tolerance.'])
% else
%     if any(sum(exceeds_tol)>0)
%         warning(sprintf('%s\n%s',[mfilename ': Deviation from original trajectory exceeds tolerance!'], ...
%             ['Max. deviation: ' num2str(max_deviation) ' pxl.']))
%     end
%     if any(sum(exceeds_tol)>max_count) && any(prctile(dst2orig,99.9) > bodylength )
%         keyboard
%         error(sprintf('%s\n%s',[mfilename ': Deviation from original trajectory exceeds tolerance!'], ...
%             ['Max. deviation: ' num2str(max_deviation) ' pxl.']))
%     end
% end
%%
% ff=1;
% fidx=34385;
% idces = fidx-400:fidx+400;%1:len;%14000:16000;
% figure
% plot(trajectories(idces,ff,1),trajectories(idces,ff,2),'-')
% hold on
% plot(tr_smooth(idces,ff,1),tr_smooth(idces,ff,2),'g-')
% plot(trajectories(fidx,ff,1),trajectories(fidx,ff,2),'o')
% plot(tr_smooth(fidx,ff,1),tr_smooth(fidx,ff,2),'g*')
