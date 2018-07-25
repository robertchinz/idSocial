 function [h, timemap_fitparameters,sf] =idSocial_displayApproachMap(edges,tmap,interval_in_seconds) 

 
 hold on
 
 [~, mth]=contourf(tmap,0:.5:interval_in_seconds,'-');
 z=tmap;
 x=edges{1}(1):(edges{1}(end)-edges{1}(1))/(size(tmap,1)-1):edges{1}(end);
 y=edges{2}(1):(edges{2}(end)-edges{2}(1))/(size(tmap,2)-1):edges{2}(end);
 [X,Y] = meshgrid(x,y);
 R = sqrt(X.^2+Y.^2);
 Th = abs(atan2(Y,X)');
 
 

%  myfun = @(params) krause_poly4fitting(X(~isnan(z(:))),Y(~isnan(z(:))),params) - z(~isnan(z(:)));
%  params0 = [.1,.1,1,10,1];
%  opts = optimset('Display','Iter','TolFun',1e-10);

 %Fit
%  try
%      timemap_fitparameters = lsqnonlin(myfun,params0,[],[],opts);
%      
%  catch
%      timemap_fitparameters=NaN(1,5);
%      disp([mfilename ': Objective function is returning undefined values at initial point. lsqnonlin cannot continue.'])
%  end
%  fitf=krause_poly4fitting(X,Y,timemap_fitparameters);
%  
%  [~, cth2]=contour(fitf,0:.5:interval_in_seconds,'-k');
%  [~, cth]=contour(fitf,0:.5:interval_in_seconds,'-');

 sf = fit([R(:), Th(:)],z(:),'poly11');
 timemap_fitparameters = [sf.p10 sf.p01 sf.p00];
 fitf=krause_poly4fittingPolar(R,Th,timemap_fitparameters);
 
 [~, cth2]=contour(fitf,0:.5:interval_in_seconds,'-k');
 [~, cth]=contour(fitf,0:.5:interval_in_seconds,'-');
 
 uistack(mth,'bottom');
 set(mth,'edgecolor','none');
caxis([0 interval_in_seconds])
 
 axis equal tight
 set(gca,'YDir','normal')

 
 colorbar
 colormap(jet(ceil(interval_in_seconds)*2))
 
 set(cth,'LineWidth',1.5)
 set(cth2,'LineWidth',2.0)
 hold off
       
 h=[cth mth];
 idSocial_auxiliaries_setImageTicks(gca,size(tmap),edges)
  
 
 end
 
function T=krause_poly4fitting(x,y,params)
% see Krause et. al., Anim. Behav., 1994, 48, 353-359.
d=(x.^2+y.^2).^.5;
a=abs(atan2(x,y));
T=params(1)*d+params(2)*a+params(3);

end

function T=krause_poly4fittingPolar(d,a,params)
% see Krause et. al., Anim. Behav., 1994, 48, 353-359.
T=params(1)*d+params(2)*a+params(3);

end