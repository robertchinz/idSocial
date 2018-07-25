% function pooled_cell=idSocial_function_wrapper_poolCellDims(orig_cell,pool_idx,dim_names_rearranged)

no_legendentries=size(pool_idx,1);

if size(orig_cell,1)==1 && no_legendentries>1
    rep_cell=cell(no_legendentries,1);
    dim_names_rearranged=cell(no_legendentries,size(dim_names_rearranged,2));
    for k=1:no_legendentries
        rep_cell{k}=orig_cell;
        dim_names_rearranged{k,:}=dim_names_rearranged;
    end
else
    rep_cell{1}=orig_cell;
%     dim_names_rearranged{1,:}=dim_names_rearranged;
end

S.type='()';
for k=1:no_legendentries
     ndims_cell=ndims(rep_cell{k});
     idxarray=[setxor(1:ndims_cell,pool_idx{k}) pool_idx{k}];
     temp_cell=permute(rep_cell{k},idxarray);
     dim_names_rearranged(k,:)=dim_names_rearranged([idxarray setxor(idxarray,1:length(dim_names))])';
     S.subs=repmat({':'},1,size(setxor(1:ndims_cell,pool_idx{k}),2)+1);
     pooled_cell=subsref(temp_cell,S);
end