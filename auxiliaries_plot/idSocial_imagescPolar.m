function [h, ax_lines, Xctrs, Yctrs]=idSocial_imagescPolar(edgesR,edgesTh,C,clims,interpolation_factor,area_normalization,totalNumber_normalization)


if nargin<5 || isempty(interpolation_factor)
    interpolation_factor=[1 1];
end

if size(edgesR,2)>1 && size(edgesR,1)==1
    edgesR = edgesR';
end
if size(edgesTh,1)>1 && size(edgesTh,2)==1
    edgesTh = edgesTh';
end
if nargin<6 || isempty(area_normalization)
    area_normalization = false;
end
if nargin<6 || isempty(totalNumber_normalization)
    totalNumber_normalization = false;
end


if area_normalization
    ar = pi*edgesR(2:end).^2 - pi*edgesR(1:end-1).^2;
    ft =     (edgesTh(2)-edgesTh(1))/(2*pi);
    atotal = repmat(ar .* ft,[1,size(edgesTh,2)-1]);
    C = C./atotal;
end
if totalNumber_normalization
    C = C/nansum(C(:));
end

if nargin<4 || isempty(clims)
    clims=[min(C(:)) max(C(:))];
end

dr = edgesR(2)-edgesR(1);
dtheta = edgesTh(2)-edgesTh(1);

r=(edgesR(1):dr/interpolation_factor(1):edgesR(end-1))';
theta = edgesTh(1):dtheta/interpolation_factor(2):edgesTh(end-1);
X = -r*cos(theta+pi/2);
Y = r*sin(theta+pi/2);
Xctrs = -(edgesR+dr/2)*cos((edgesTh+dtheta/2)+pi/2);
Yctrs = (edgesR+dr/2)*sin((edgesTh+dtheta/2)+pi/2);
Cinterp=NaN(numel(r),numel(theta));
cnt = [1 1];
% Xctrs= NaN(size(C));
% Yctrs= NaN(size(C));
for rw=1:interpolation_factor(1):size(Cinterp,1)-1
    for cl=1:interpolation_factor(2):size(Cinterp,2)-1
        Cinterp(rw:rw+interpolation_factor(1)-1,cl:cl+interpolation_factor(2)-1) = C(cnt(1),cnt(2));
        cnt(2) = cnt(2)+1;
    end
    cnt(2)=1;
    cnt(1) = cnt(1)+1;
end

 
 h=pcolor(X,Y,Cinterp);

% gridl_h = axes('Position',get('Position',gca));
hold on
set(h, 'EdgeColor', 'none');
ax_r=NaN(1,numel(edgesR)-1);
for k = 1:numel(edgesR)-1
    ax_r(k)=plot(edgesR(k)*cos(0:2*pi/500:2*pi),... 
    edgesR(k)*sin(0:2*pi/500:2*pi),...
    'k','LineWidth',get(gca,'LineWidth'));
end
ax_Th=NaN(1,numel(edgesTh)-1);
for k = 1:numel(edgesTh)
    ax_Th(k)=plot([0 edgesR(end-1)*cos(edgesTh(k))],... 
    [0 edgesR(end-1)*sin(edgesTh(k))],...
    'k','LineWidth',get(gca,'LineWidth'));
end
ax_lines = [ax_r ax_Th];

% warning([mfilename ': Fixe LineWidth for thesis!!!!'])
% polar_linewidth =1;
% set(ax_lines,'LineWidth',polar_linewidth);
if all(isnan(clims))
    clims = [-1 1];
end
caxis(clims)
set(gca,'XLim',[-edgesR(end-1) edgesR(end-1)],'YLim',[-edgesR(end-1) edgesR(end-1)])
