function input_data = idSocial_changeTrajectoryPathInData(input_data,old_path,new_path)


if numel(input_data)>1
    data_siz = size(input_data);
    if ndims(data_siz)==2; data_siz=[data_siz 1]; end
    for gr = 1 : data_siz(1)
        for day = 1 : data_siz(2)
            for tr = 1 : data_siz(3)
                if isfield(input_data(gr,day,tr),'movementdata')
                    if ~isempty(input_data(gr,day,tr).movementdata)
                        input_data(gr,day,tr).movementdata = ...
                            strrep(input_data(gr,day,tr).movementdata,old_path,new_path);
                    end
                end
                if isfield(input_data(gr,day,tr),'rand_movementdata')
                    if ~isempty(input_data(gr,day,tr).rand_movementdata)
                        input_data(gr,day,tr).rand_movementdata = ...
                            strrep(input_data(gr,day,tr).rand_movementdata,old_path,new_path);
                    end
                end
            end
        end
    end
else
    data_siz = size(input_data.movementdata);
    if ndims(data_siz)==2; data_siz=[data_siz 1]; end
    for gr = 1 : data_siz(1)
        for day = 1 : data_siz(2)
            for tr = 1 : data_siz(3)
                if isfield(input_data,'movementdata')
                    if ~isempty(input_data.movementdata{gr,day,tr})
                        input_data.movementdata{gr,day,tr} = ...
                            strrep(input_data.movementdata{gr,day,tr},old_path,new_path);
                    end
                end
                if isfield(input_data,'rand_movementdata')
                    if ~isempty(input_data.rand_movementdata{gr,day,tr})
                        input_data.rand_movementdata{gr,day,tr} = ...
                            strrep(input_data.rand_movementdata{gr,day,tr},old_path,new_path);
                    end
                end
            end
        end
    end
    
    
end