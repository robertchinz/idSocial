function [tr_smooth, vel, acc, tol_out, max_dist, d2oAll, no_frames_off] = idSocial_auxiliaries_smoothTrajectoryAdaptiveMovAvg(trajectories,tol_frame)
filter_outliers =true;
display_plot =true;% false;
if nargin<2 || isempty(tol_frame)
    tol_frame = 3;
end
% if nargin<3 || isempty(bodylength)
%     bodylength = inf;
% end



weight_increment = .2;
err = .2*tol_frame;
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
d2oAll = NaN(len,size(trajectories,2),1);
dst2orig = NaN(len,size(trajectories,2),3);
zero_idces = NaN(len,size(trajectories,2),1);
zero_idcesx = false(len,size(trajectories,2),1); % obsolete
zero_idcesy = false(len,size(trajectories,2),1); % obsolete    
zero_idcesz = NaN(len,size(trajectories,2),1);
tr_smooth = NaN(size(trajectories));
vel = NaN(size(trajectories));
acc = NaN(size(trajectories));
scrsz = get(0,'ScreenSize');
if display_plot
    fh=figure('Name','Adaptive Moving Information','Position',[scrsz(3)*.2 scrsz(4)*.2 scrsz(3)*.5 scrsz(4)*.3 ]);
end

deltaX = inf;
deltaY = inf;
deltaZ = inf;
tol_out=NaN(no_fish,1);
max_dist=NaN(no_fish,1);
no_frames_off = NaN(no_fish,1); 
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
    tol = 30;%(0.00051*len*tol_frame).^2;
    
    tolx_prev = 0;
    toly_prev = 0;
    tolx=tol;
%     toly=tol;
    count=1;
    deltaX = prctile(d2ox,99)-tol_frame;
    %     while (~all(d2ox(~isnan(d2ox))<tol_frame) || prctile(d2ox,99)-tol_frame > 0 || prctile(d2ox,99)-tol_frame< -err)  && count<max_count
    while (all(isinf(d2ox)) || deltaX > 0 || deltaX< -err)  && count<max_count
        
        dst_velx(2:end,ff) = abs(trajectories(1:end-1,ff,1)-trajectories(2:end,ff,1));
        nan_velx(:,ff) = isnan(dst_velx(:,ff));
        dst_velx(nan_velx(:,ff),ff)=0;
        dst_velx(:,ff) = vertcat(0,dst_velx(2:end,ff));
        zero_idcesx(:,ff) = dst_velx(:,ff)==0;
        
        %         keyboard
        xx=1:len; xx=xx(~zero_idcesx(:,ff));
        %         xx=cumsum(dst_velx(:,ff)); xx=xx(~zero_idcesx(:,ff));
        %                 xx=cumsum(dst_vel_norm(:,ff)); xx=xx(~zero_idces(:,ff));
 
      
            trX=smooth(trajectories(:,ff,1),tolx,'moving');
            trY=smooth(trajectories(:,ff,2),tolx,'moving');
      
        
        tr_smooth(:,ff,1) = trX;
        tr_smooth(:,ff,2) = trY;
        
       deltaXprev = deltaX;
        d2ox = sqrt(sum((trajectories(:,ff,1:2)-tr_smooth(:,ff,1:2)).^2,3));%abs(squeeze(trajectories(:,ff,1))-squeeze(tr_smooth(:,ff,1)));
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
%         [ tolx deltaX]
        if deltaX > 0 || deltaX < - err
            if deltaX < - err && deltaXprev > 0
                tolx = (tolx_prev2+tolx_prev)/2;
            elseif deltaX > 0 && deltaXprev < - err
                tolx = (tolx_prev2+tolx_prev)/2;
            elseif deltaX < - err && deltaXprev < - err
                tolx = tolx_prev*2;
            elseif deltaX > 0  && deltaXprev > 0
                tolx = tolx_prev/2;%nansum(d2ox(d2ox<tol_frame)) + (nansum(d2ox>tol_frame)*tol_frame/2) ;%nanmedian(d2o)*len;
            else
                keyboard
            end
            tolx = round(tolx);

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
%             drawnow
        end
    end
    for fr = 2:size(trajectories,1)
        if zero_idcesx(fr,ff) && ~nan_velx(fr,ff)
            tr_smooth(fr,ff,1) = tr_smooth(fr-1,ff,1);
        end
    end
    if count==max_count
        disp('Max. number of repetitions reached.')
    end
    disp([num2str(sum(d2ox>tol_frame)) ' Frames exceed max. distance. Max. distance: ' num2str(max(d2ox))])
%     count=1;
%     while (all(isinf(d2oy)) || prctile(d2oy,99)-tol_frame > 0 || prctile(d2oy,99)-tol_frame< -err)  && count<max_count
%         
%         %     while (~all(d2oy(~isnan(d2oy))<tol_frame) || prctile(d2oy,99)-tol_frame > 0 || prctile(d2oy,99)-tol_frame< -err)  && count<max_count
%         dst_vely(2:end,ff) = abs(trajectories(1:end-1,ff,2)-trajectories(2:end,ff,2));
%         nan_vely(:,ff) = isnan(dst_vely(:,ff));
%         dst_vely(nan_vely(:,ff),ff)=0;
%         dst_vely(:,ff) = vertcat(0,dst_vely(2:end,ff));
%         zero_idcesy(:,ff) = dst_vely(:,ff)==0;
%         yy=1:len;yy=yy(~zero_idcesy(:,ff));
%         %         yy=cumsum(dst_vely(:,ff)); yy=yy(~zero_idcesy(:,ff));
% %         yy=cumsum(dst_vel_norm(:,ff)); yy=yy(~zero_idces(:,ff));
%         try
% %                 [spY, trY]=...
% %                     spaps(yy,trajectories(~zero_idcesy(:,ff),ff,2),toly,wy(~zero_idcesy(:,ff)),spl_deg);
%                 trY=smooth(trajectories(:,ff,2),toly,'moving');
% 
%         catch
%             keyboard
%         end
%         tr_smooth(:,ff,2) = trY;
% %         [spY, trY]=...
% %             spaps(yy,trajectories(~zero_idces(:,ff),ff,2),toly,wy(~zero_idces(:,ff)),spl_deg);
% %         tr_smooth(~zero_idces(:,ff),ff,2) = trY;
%         
%         deltaYprev = deltaY;
%         d2oy = abs(squeeze(trajectories(:,ff,2))-squeeze(tr_smooth(:,ff,2)));
%         deltaY = max(d2oy)-tol_frame;
%         deltaY = prctile(d2oy,99)-tol_frame;
%         
%         %         wy = wy + ...
%         %             sign((d2oy' - tol_frame))*weight_increment;
%         %         wy = max(wy,.1 * ones(1,len));
%         %         wy = wy/sum(wy)*len;
%         %         toly = nansum(d2oy(d2oy<tol_frame)) + (nansum(d2oy>tol_frame)*tol_frame/2) ;%nanmedian(d2o)*len;
%         
%         toly_prev2=toly_prev;
%         toly_prev=toly;
%         if deltaY > 0 || deltaY < - err
%             if deltaY < - err && deltaYprev > 0
%                 toly = (toly_prev2+toly_prev)/2;
%             elseif deltaY > 0 && deltaYprev < - err
%                 toly = (toly_prev2+toly_prev)/2;
%             elseif deltaY < - err && deltaYprev < - err
%                 toly = toly_prev*2;
%             elseif deltaY > 0  && deltaYprev > 0
%                 toly = toly_prev/2;%nansum(d2ox(d2ox<tol_frame)) + (nansum(d2ox>tol_frame)*tol_frame/2) ;%nanmedian(d2o)*len;
%                 
%             end
%             toly = round(toly);
%         end
%         
%         count = count + 1;
%         if display_plot
%             figure(fh);
%             subplot(1,2,1);hist(d2oy,200); hold on; plot([tol_frame tol_frame],get(gca,'YLim'),'k');
%             plot([prctile(d2oy,99) prctile(d2oy,99)],get(gca,'YLim'),'r');
%             plot([tol_frame-err tol_frame-err],get(gca,'YLim'),'g');
%             hold off;
%             set(gca,'XLim',[0 tol_frame*1.3])
%             title({['Y: Individual ' num2str(ff) '(tolfr ' num2str(tol_frame) '), Repetition #' num2str(count)];['Delta=' num2str(max(d2oy)) 'Rel. Error: ' num2str(sum((d2oy'.*wy).^2)/(tol_frame*len))]});
%             subplot(1,2,2); plot(wy); title('Weights y')
%             drawnow
%         end
%     end
%     for fr = 2:size(trajectories,1)
%         if zero_idcesy(fr,ff) && ~nan_vely(fr,ff)
%             tr_smooth(fr,ff,2) = tr_smooth(fr-1,ff,2);
%         end
%     end
%     if count==max_count
%         disp('Max. number of repetitions reached for Y.')
%     end
%     disp(['Y: ' num2str(sum(d2oy>tol_frame)) ' Frames exceed max. distance. Max. distance: ' num2str(max(d2oy))])
%     

    tolz = NaN;
    tol_out(ff,1) = tolx;
    max_dist(ff,1) = max(d2ox);
    no_frames_off(ff,1) = sum(d2ox>tol_frame) ;
    d2oAll(:,ff,1) = d2ox;
   
    
    
   
 dst2orig(:,ff,1) = d2ox ;
    
end
if display_plot
%     close(fh)
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
