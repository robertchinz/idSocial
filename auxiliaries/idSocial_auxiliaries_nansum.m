function nsum = idSocial_auxiliaries_nansum(vec,dim,weights)

if nargin <2 || isempty(dim)
    dim = [];
end

if ~isempty(vec)
    nsum = nansum(vec,dim);
    % end
    % If == 0, check if this is due to "all nans"
    
    if isempty(dim)
        dim = max([1 find(size(vec)>1,1,'first')]);
    end
    try
        
        nanidx = all(isnan(vec),dim);
        if any(nanidx(:))
            
            nsum(nanidx) = NaN;
        end
        
        
    catch
        keyboard
    end
else
    nsum = NaN;
end