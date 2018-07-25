function [H p CI bootstrpStat sampStat]=idSocial_bootstrap(vec1,vec2,nrep,alpha,statfuncstring,tail)
% Idea from http://courses.washington.edu/matlab1/Bootstrap_examples.html

if nargin<5 || isempty(statfuncstring)
    statfuncstring='Mean';
end
if nargin<3 || isempty(nrep)
    nrep=100;
end

if nargin<4 || isempty(alpha)
    alpha=.05;
end

if nargin<6 || isempty(tail)
    tail='both';
end
if iscell(vec1); vec1=vec1{1}; end
if iscell(vec2); vec2=vec2{1}; end

if size(vec1,2)==1; vec1=vec1'; end
if size(vec2,2)==1; vec2=vec2'; end



len1=size(vec1,2);
len2=size(vec2,2);
minlen=min(len1,len2);

if strcmpi(statfuncstring,'mean')
    statfunc = @(x,y) nanmean(x)-nanmean(y);
elseif strcmpi(statfuncstring,'median')
    statfunc = @(x,y) nanmedian(x)-nanmedian(y);
elseif strcmpi(statfuncstring,'circ_var')
    statfunc = @(x,y) idSocial_circ_var(x)-idSocial_circ_var(y);
end
sampStat=statfunc(vec1,vec2);

if minlen*nrep<100
    rvec1=vec1(ceil(rand(minlen,nrep)*len1));
    rvec2=vec2(ceil(rand(minlen,nrep)*len2));
    bootstrpStat=statfunc(rvec1,rvec2);
else
    
    bootstrpStat=NaN(nrep,1);
    
    stp=min(nrep,20);
    for k=1:nrep/stp
        rd=rand(minlen,stp);
%         rvec1=vec1(ceil(rd*len1));
%         rvec2=vec2(ceil(rd*len2));
        rvec1=vec1(ceil(rd*len1));
        rvec2=vec2(ceil(rd*len2));
%         disp(['x : ' num2str(nanmean(rvec1(:))) ', ' num2str(nanstd(rvec1(:)))])
%         disp(['y : ' num2str(nanmean(rvec2(:))) ', ' num2str(nanstd(rvec2(:))) ])
        bootstrpStat((k-1)*stp+1:k*stp)=statfunc(rvec1,rvec2);
    end

%     bstat1 =bootstrp(nrep,@mean,vec1);
%     bstat2 =bootstrp(nrep,@mean,vec2);
%     bootstrpStat=bstat1-bstat2;
end


CI=prctile(bootstrpStat,[100*alpha/2, 100*(1-alpha/2)]);
H=CI(1)>0 | CI(2)<0;
p=1;
if strcmpi(tail,'both')
    if CI(1)>0
        p=sum(bootstrpStat<=0)/nrep;
    elseif CI(2)<0
        p=sum(bootstrpStat>=0)/nrep;
    end
end
if strcmpi(tail,'right') % alternative hypothesis states mean/median of x is greater than the mean/median of y
   p=sum(bootstrpStat<=0)/nrep;
end
if strcmpi(tail,'left') % alternative hypothesis states mean/median of x is less than the mean/median of y
        p=sum(bootstrpStat>=0)/nrep;
end
% keyboard
%%
% figure
% xx = min(bootstrpStat):.01:max(bootstrpStat);
% hist(bootstrpStat,xx);
% hold on
% ylim = get(gca,'YLim');
% h1=plot(sampStat*[1,1],ylim,'y-','LineWidth',2);
% h2=plot(CI(1)*[1,1],ylim,'r-','LineWidth',2);
% plot(CI(2)*[1,1],ylim,'r-','LineWidth',2);
% h3=plot([0,0],ylim,'b-','LineWidth',2);
% xlabel('Difference between means');
% 
% decision = {'Fail to reject H0','Reject H0'};
% title(decision(H+1));
% legend([h1,h2,h3],{'Sample mean',sprintf('%2.0f%% CI',100*alpha),'H0 mean'},'Location','NorthWest');