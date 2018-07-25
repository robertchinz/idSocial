function [rand_check, no_fishORIG] = idSocial_auxiliaries_trRandCheck(trajectory)

if ndims(trajectory) ==4
    no_focals = size(trajectory,2);
    no_neighbors = size(trajectory,3);
    rand_check = NaN(no_focals,no_neighbors);
if ndims(trajectory)==4
    for ff = 1:no_focals
        for nf = 1:no_neighbors
            rand_check(ff,nf) = all(isnan(trajectory(:,ff,nf,1)) & isinf(trajectory(:,ff,nf,2)));
        end
    end
    if any(rand_check(:)==1)
        no_fishORIG = find(any(rand_check==1,1),1,'first')-1;
    else
        no_fishORIG = no_focals;
    end
else
    rand_check = false(no_focals,no_neighbors);
end
elseif ndims(trajectory) ==3
    no_focals = size(trajectory,2);
    rand_check = all(isnan(trajectory(:,:,1)),1) & all(isinf(trajectory(:,:,1)),1);
    if any(rand_check(:)==1)
        no_fishORIG = find(any(rand_check==1,1),1,'first')-1;
    else
        no_fishORIG = no_focals;
    end
end
