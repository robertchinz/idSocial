function orig_cell = idSocial_auxiliaries_loadResults(orig_cell,files)

if strcmp(orig_cell(end-4),'_')
    
    if nargin>1 && ~isempty(files)
        first_file = files(1);
    else
        first_file = 1;
    end
    try
        disp(['Loading ' orig_cell(1:end-4) '_' num2str(first_file) '.mat ...'])
        loaddat=load([orig_cell(1:end-4) '_' num2str(first_file) '.mat']);
%         no_files = loaddat.output_new.no_files;
        if nargin>1 && ~isempty(files)
            file_idces = files;
        else
            file_idces=1:loaddat.output_new.no_files;
        end
        size_file = loaddat.output_new.size;
        act_size = size(loaddat.output_new.output,2);
        orig_cell_data = cell(1,prod(size_file));
        if isfield(loaddat,'output_new')
            orig_cell_data(1:act_size) = loaddat.output_new.output;
        elseif isfield(loaddat,'output_rand_new')
            orig_cell_data(1:act_size) = loaddat.output_rand_new.output;
        end
        %             orig_cell_data(1:act_size) = loaddat.output_new.output;
        disp('Done.')

        act_idx = act_size + 1;
        for k=file_idces(2:end)
            disp([orig_cell(1:end-4) '_' num2str(k) '.mat ...'])

            loaddat=load([orig_cell(1:end-4) '_' num2str(k) '.mat']);
            
            %% Filter out all-nan dims, mainly to save some space
            if ~any(cellfun(@(x) iscell(x), loaddat.output_new.output))
                allnan=cellfun(@(x) ~isempty(x)&&all(isnan(x)), loaddat.output_new.output);
                loaddat.output_new.output(allnan) = {[]};
            end
            
            act_size = size(loaddat.output_new.output,2);
            orig_cell_data(act_idx:act_idx + act_size -1) = ...
                loaddat.output_new.output;
            act_idx = act_idx + act_size;
            disp('Done.')

        end
        orig_cell = reshape(orig_cell_data,size_file);

        clear orig_cell_data;
    catch
        warning(mfilename)
        keyboard
    end
else
    try
        loaddat=load(orig_cell);
    catch
        warning(mfilename)
        keyboard
    end
    if isfield(loaddat,'output_new')
        orig_cell=loaddat.output_new;
    elseif isfield(loaddat,'output_rand_new')
        orig_cell=loaddat.output_rand_new;
    end
end