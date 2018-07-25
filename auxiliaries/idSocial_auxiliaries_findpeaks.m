function [sortvals,sortlocs] = idSocial_auxiliaries_findpeaks(vec,sortstr,no_peaks)

if nargin < 2 || isempty(sortstr)
        sortstr = [];
end
if nargin <3 || isempty(no_peaks)
    no_peaks = [];
end

dx = diff(vec);
locs=find(dx(1:end-1)>0 & dx(2:end)<0) +1;
vals = vec(locs);

if ~isempty(sortstr) && (strcmpi(sortstr,'ascend') || strcmpi(sortstr,'descend'))
    [sortvals, sortlocs] = sort(vals,sortstr);
    sortlocs = locs(sortlocs);
else 
    sortvals = vals;
    sortlocs = locs;
end

if ~isempty(no_peaks) && ~isempty(sortvals)
    sortvals = sortvals(1:no_peaks);
    sortlocs = sortlocs(1:no_peaks);
end

% figure; plot(hitemp); hold on; plot(idx,hitemp(locs),'+')