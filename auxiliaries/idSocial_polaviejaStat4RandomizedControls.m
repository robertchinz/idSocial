function [H p CI stat_vec sampStat]=idSocial_polaviejaStat4RandomizedControls(vec1,vec2,no_rep,alpha,statfuncstring,tail,mode_edges,mode_bandwidth,mode_method)

if nargin<5 || isempty(statfuncstring)
    statfuncstring='Mean';
end
if nargin<3 || isempty(no_rep)
    no_rep=100;
end

if nargin<4 || isempty(alpha)
    alpha=.05;
end

if nargin<6 || isempty(tail)
    tail='both';
end

if nargin<7 || isempty(mode_edges)
    mode_edges=[];
end

if nargin<8 || isempty(mode_bandwidth)
    mode_bandwidth=[];
end

if nargin<9 || isempty(mode_method)
    mode_method=[];
end


if iscell(statfuncstring) && size(statfuncstring,2) == 2
    statfunc_param = statfuncstring{1,2};
    statfuncstring = statfuncstring{1,1};
else
    statfunc_param = [];
end

if strcmpi(statfuncstring,'mean')
    statfunc = @(x) nanmean(x);
elseif strcmpi(statfuncstring,'median')
    statfunc = @(x) nanmedian(x);
elseif strcmpi(statfuncstring,'var')
    statfunc = @(x) nanvar(x);
elseif strcmpi(statfuncstring,'circ_var')
    statfunc = @(x) idSocial_circ_var(x);
elseif strcmpi(statfuncstring,'angle_ratio')
    if isempty(statfunc_param)
        statfunc_param = pi/2;
    end
    statfunc = @(x) idSocial_angle_ratio(x,[],statfunc_param);
elseif strcmpi(statfuncstring,'positive_ratio')
    if isempty(statfunc_param)
        statfunc_param = 0;
    end
    statfunc = @(x) idSocial_positive_ratio(x,[],statfunc_param);
elseif strcmpi(statfuncstring,'mode') && ~isempty(mode_edges) && ~isempty(mode_bandwidth) && ~isempty(mode_method)
    statfunc = @(x) idSocial_distribution_mode(x,mode_edges,mode_bandwidth,mode_method);  
elseif strcmpi(statfuncstring,'1stmode') && ~isempty(mode_edges) && ~isempty(mode_bandwidth) && ~isempty(mode_method)
    statfunc = @(x) idSocial_distribution_mode(x,mode_edges,mode_bandwidth,mode_method,'first'); 
elseif strcmpi(statfuncstring,'1stmode_height') && ~isempty(mode_edges) && ~isempty(mode_bandwidth) && ~isempty(mode_method)
    statfunc = @(x) idSocial_distribution_modeHeight(x,mode_edges,mode_bandwidth,mode_method,'first');  
else
    error([mfilename ': Unknown test statistics.'])
end

% vec1=squeeze(val(group,:,:));
% vec2=squeeze(val(group2,:,:));

% Eliminate NaN frames
vec1=vec1(~isnan(vec1));
vec2=vec2(~isnan(vec2));

no_frames1=numel(vec1);
no_frames2=numel(vec2);

try
if size(vec1,2)==1 && size(vec1,1) >= 1
    vec1 = vec1';
end
if size(vec2,2)==1 && size(vec2,1) >= 1
    vec2 = vec2';
end
catch
    keyboard
end
% no_rep = 1000;

stat_vec=NaN(1,no_rep);
% keyboard

try
sampStat=statfunc(vec1);
% keyboard
catch
    keyboard
end
if ~(no_frames1==0 || no_frames2==0) && no_frames1<no_frames2
    for rep = 1:no_rep
        no_framesSample = no_frames1;
      
        idces = randsample(no_frames2,no_framesSample);
        vec1rand = vec2(idces);
        try
        stat_vec(rep)=statfunc(vec1rand);
        catch
            keyboard
        end
    end
end
% figure; hist(stat_vec,200)
% hold on; plot([sampStat sampStat],get(gca,'YLim'),'r')
CI=prctile(stat_vec,[100*alpha/2, 100*(1-alpha/2)]);

p=1;
if strcmpi(tail,'both')
%     Is this correct?
    H=CI(1)>sampStat | CI(2)<sampStat;
    if CI(1)>sampStat
        p=(sum(stat_vec<sampStat))/(no_rep);
    elseif CI(2)<sampStat
        p=(sum(stat_vec>sampStat))/(no_rep);
    end
     H = p<alpha;
end
if strcmpi(tail,'right') % alternative hypothesis states mean/median/... of x is greater than the mean/median of y
%     H= CI(2)<sampStat;
    if ~isnan(sum(stat_vec)) && ~isnan(sampStat)
        p=(sum(stat_vec>sampStat))/(no_rep);
    else
        p=1;
    end
       H = p<alpha;
end
if strcmpi(tail,'left') % alternative hypothesis states mean/median/... of x is less than the mean/median of y
%     H= CI(1)>sampStat;
    %     p=(sum(stat_vec<=sampStat)+1)/(no_rep+1);
    if ~isnan(sum(stat_vec)) && ~isnan(sampStat)
        p=(sum(stat_vec<sampStat))/(no_rep);
    else
        p=1;
    end
     H = p<alpha;
end
H = H & p<alpha;
% p = max(p,1/no_rep); % Even if the sample statistics lies outside 
% the confidence interval, there is still the probability of
% 1/no_rep to draw the actual sample distribution from the
% joined distribution!
% http://faculty.washington.edu/kenrice/sisg/SISG-08-06.pdf