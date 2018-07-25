function input_data = idSocial_changeOutputPathInData(input_data,old_path,new_path)

fnames = fieldnames(input_data(1,1,1));


for fn = 1:numel(fnames)
    if isfield(input_data(1,1,1).(fnames{fn}),'output') && ...
            ~isempty(input_data(1,1,1).(fnames{fn}).output) && ...
            ischar(input_data(1,1,1).(fnames{fn}).output)
        input_data(1,1,1).(fnames{fn}).output = ...
            strrep(input_data(1,1,1).(fnames{fn}).output,old_path,new_path);
    end
end
