function [filter_cell, input_data] = idSocial_auxiliaries_output2filter(input_data,act_method)
if isfield(input_data,act_method)
    load(input_data.(act_method).output)
    
    
    
    no_groups = size(output_new,1);
    no_days = size(output_new,2);
    no_trials = size(output_new,3);
    no_time = size(output_new,4);
    
    filter_cell = cell(no_groups,no_days,no_trials,no_time);
    
    for gr = 1:no_groups
        for dd = 1:no_days
            for tr = 1:no_trials
                for tm = 1:no_time
                    
                    act_data=squeeze(output_new(gr,dd,tr,tm,:,:));
                    
                    no_frames=max(cellfun(@(x) numel(x),act_data(:)));
                    no_fish = size(act_data,1);
                    
                    data_array = NaN(no_frames, no_fish, no_fish);
                    for ff = 1:no_fish
                        for nf = 1:no_fish
                            
                            if ~isempty(act_data{ff,nf})
                                if ~(islogical(act_data{ff,nf}) || all(act_data{ff,nf}==1 | act_data{ff,nf} ==0)) 
                                    error([mfilename ': ' act_method ' does not contain logicals.'])
                                end
                                data_array(:,ff,nf) = act_data{ff,nf};
                            end
                        end
                    end
                    filter_cell{gr,dd,tr,tm} = data_array;
                    
                end
            end
        end
    end
else
    error([mfilename ': ' act_method ' not found.'])
end



