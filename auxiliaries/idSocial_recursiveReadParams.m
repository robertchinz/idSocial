function rec_cell = idSocial_recursiveReadParams(rec_cell,input_data,def_options,act_method)

if isa(rec_cell,'cell')
    for k = 1:numel(rec_cell)
 
        rec_cell{k} = idSocial_recursiveReadParams(rec_cell{k},input_data,def_options,act_method);
    end
else
    rec_cell.test = 2;  
    [~, rec_cell]=idSocial_readparams(input_data,rec_cell,def_options,act_method);
    return
end
end
