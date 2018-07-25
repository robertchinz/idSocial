function act_ndims = idSocial_recursiveCellSize(inpcell)

act_ndims=0;

rec_cell=inpcell;
% max_size=NaN(1,2);

while isa(rec_cell,'cell')
    act_ndims = act_ndims + 1;
%     max_size(act_ndims,:) = size(rec_cell);
    rec_cell = rec_cell{1};
end

end