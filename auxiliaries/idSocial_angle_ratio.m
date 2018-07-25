function aratio = idSocial_angle_ratio(vec,dim,ang)

if nargin<2 || isempty(dim)
    siz=size(vec);
    firstidx = find(siz>1,1,'first');
    if isempty(firstidx)
        dim = 1;
    else
        dim = firstidx;
    end
end

if nargin<3 || isempty(ang)
    ang = pi/2;
end

if numel(ang)==2
    aratio = nansum(vec<ang(1) & ~isnan(vec) | vec>ang(2) & ~isnan(vec),dim)./sum(~isnan(vec),dim);
else
    % vec=vec(~isnan(vec));
    aratio = nansum(vec<ang & ~isnan(vec),dim)./sum(~isnan(vec),dim);
    % if any(isnan(squeeze(aratio))); keyboard; end
end