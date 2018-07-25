function combi=vectorsize2indexcombinations(vec_size)

sets=cell(length(vec_size),1);
for k=1:size(sets,1)
    sets{k}=1:vec_size(k);
end
c = cell(1, numel(sets));
[c{:}] = ndgrid( sets{:} );
combi = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );