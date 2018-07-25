function aratio = idSocial_positive_ratio(vec,dim,thresh)

if nargin<2 || isempty(dim)
    siz=size(vec);
    firstidx = find(siz>1,1,'first');
    if isempty(firstidx)
        dim = 1;
    else
        dim = firstidx;
    end
end

if nargin<3 || isempty(thresh)
    thresh = 0;
end

aratio = nansum(vec>thresh & ~isnan(vec),dim)./sum(~isnan(vec)&vec~=thresh,dim);