function [outTableCell,outTable] = idSocial_auxiliaries_getOutputTable(output_plot,output_plotRAND,method_name)

if nargin<2 || isempty(output_plotRAND)
    output_plotRAND=[];
end
if nargin<3 || isempty(method_name)
    method_name='';
end



% outp = method.output_plot;
outp =output_plot;
outpRand =output_plotRAND;

XTicks = [];
XTickLabel = [];
YTicks = [];
YTickLabel = [];
% if isfield(outp,'ticks')
%     XTicks = outp.XTick;
% end
if isfield(outp,'XTick')
    XTicks = outp.XTick;
end
if isfield(outp,'XTickLabel')
    XTickLabel = outp.XTickLabel;
end
if isfield(outp,'YTick')
    YTicks = outp.YTick;
end
if isfield(outp,'YTickLabel')
    YTickLabel = outp.YTickLabel;
end
if isfield(outp,'xlabel')
    xlabel = outp.xlabel;
end
if isfield(outp,'ylabel')
    ylabel = outp.ylabel;
end
    
no_filter_entries = numel(outp.data);

if isempty(output_plotRAND)
    outTableCell = cell(8,no_filter_entries);
else
    outTableCell = cell(12,no_filter_entries);
end

% header_string = {'Description' 'Results' 'Results before statistics' 'Results (Distribution)' 'Results before statistics (Distribution)'};
row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; ...
    'Results (RAND)'; 'Results Raw (RAND)'; 'Results Distribution (RAND)'; 'Results Raw Distribution (RAND)'; 'Axis'};
header_stringTable = outp.legendstring;
    
for fe = 1:no_filter_entries
    
    outTableCell{1,fe} = squeeze(outp.data_string{fe});
    
    data = squeeze(outp.data{fe});
    if size(data,1)>1 && size(data,2)==1 % If a vector, make it row vector ("x-direction")
        data = data';
    end
    outTableCell{2,fe} = data;
    
    data_sign = squeeze(outp.data_sign{fe});
    if size(data_sign,2) ~= size(squeeze(outp.data{fe}),1) && ...
            size(data_sign,1) == size(squeeze(outp.data{fe}),1)
        data_sign = data_sign';
    end
    
    try
        outTableCell{3,fe} = data_sign(~all(isnan(data_sign),2),:);
    catch
        outTableCell{3,fe} = data_sign;
    end
    
    data = squeeze(outp.data_Median{fe});
    if size(data,1)>1 && size(data,2)==1 % If a vector, make it row vector ("x-direction")
        data = data';
    end
    outTableCell{4,fe} = data;
    
    data_sign = squeeze(outp.data_signMedian{fe});
    if size(data_sign,2) ~= size(squeeze(outp.data_Median{fe}),1) && ...
            size(data_sign,1) == size(squeeze(outp.data_Median{fe}),1)
        data_sign = data_sign';
    end
    outTableCell{5,fe} = data_sign;
    
    outTableCell{6,fe} = {xlabel; ylabel};
    outTableCell{7,fe} = {XTicks; YTicks};
    outTableCell{8,fe} = {XTickLabel; YTickLabel};

    row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; 'Axis'; 'Ticks'; 'TickLabels'};
    
end

if ~isempty(output_plotRAND)
    for fe = 1:no_filter_entries
        data = squeeze(outpRand.data{fe});
        if size(data,1)>1 && size(data,2)==1 % If a vector, make it row vector ("x-direction")
            data = data';
        end
        outTableCell{9,fe} = data;
        
        data_sign = squeeze(outpRand.data_sign{fe});
        if size(data_sign,2) ~= size(squeeze(outpRand.data{fe}),1) && ...
                size(data_sign,1) == size(squeeze(outpRand.data{fe}),1)
            data_sign = data_sign';
        end
        try
            outTableCell{10,fe} = data_sign(~all(isnan(data_sign),2),:);
        catch
            outTableCell{10,fe} = data_sign;
        end
        
        data = squeeze(outpRand.data_Median{fe});
        if size(data,1)>1 && size(data,2)==1 % If a vector, make it row vector ("x-direction")
            data = data';
        end
        outTableCell{11,fe} = data;
        
        data_sign = squeeze(outpRand.data_signMedian{fe});
        if size(data_sign,2) ~= size(squeeze(outpRand.data{fe}),1) && ...
                size(data_sign,1) == size(squeeze(outpRand.data{fe}),1)
            data_sign = data_sign';
        end
        outTableCell{12,fe} = data_sign;
        
        
%         outTableCell{13,fe} = {xlabel; ylabel};
% %         row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; ...
% %             'Results (RAND)'; 'Results Raw (RAND)'; 'Results Distribution (RAND)'; 'Results Raw Distribution (RAND)'; 'Axis';};
%         
%         outTableCell{14,fe} = {XTicks; YTicks};
% %         row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; ...
% %             'Results (RAND)'; 'Results Raw (RAND)'; 'Results Distribution (RAND)'; 'Results Raw Distribution (RAND)'; 'Axis'; 'Ticks';};
%         
%         outTableCell{14,fe} = {XTickLabel; YTickLabel};
        row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; 'Axis'; 'Ticks'; 'TickLabels'; ...
            'Results (RAND)'; 'Results Raw (RAND)'; 'Results Distribution (RAND)'; 'Results Raw Distribution (RAND)'; };
        
        
        %         if  isfield(outp,'xaxis') && ~isempty(outp.) % Ticks = Xticks
        %             outTableCell{11,fe} = XTickLabel;
%             row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; ...
%                 'Results (RAND)'; 'Results Raw (RAND)'; 'Results Distribution (RAND)'; 'Results Raw Distribution (RAND)'; 'Axis'; 'Axis Label'};
%         end
       
    end
end
% if isempty(output_plotRAND) 
%     outTableCell{6,fe} = {xlabel; ylabel};
%     outTableCell{7,fe} = {XTicks; YTicks};
%     outTableCell{8,fe} = {XTickLabel; YTickLabel};
% 
%     row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; 'Axis'; 'Ticks'; 'TickLabels'};
% 
% end
% else
%     outTableCell(6,:) = [];
% end
try
    %     outTable = cell2table(outTableCell,...
    %     'RowNames',row_stringTable(1:size(outTableCell,1)),...
    %     'VariableNames',cellfun(@(str) str(ismember(str,['A':'Z' 'a':'z' '0':'9'])),header_stringTable(:),'UniformOutput',false));
    if ~verLessThan('matlab', '14.0.1')
        outTable = cell2table(outTableCell,...
            'RowNames',row_stringTable(1:size(outTableCell,1)),...
            'VariableNames',header_stringTable(:));
        if ~isempty(method_name)
            outTable.Properties.Description = ['Results for ' method_name];
        end
        
    
    else
        outTable = [row_stringTable(1:size(outTableCell,1)) outTableCell];
        outTable = vertcat([cell(1,1) header_stringTable(:)'],outTable);
    end
catch
    keyboard
end
% function [outTableCell,outTable] = idSocial_auxiliaries_getOutputTable(output_plot,output_plotRAND,method_name)
% 
% if nargin<2 || isempty(output_plotRAND)
%     output_plotRAND=[];
% end
% if nargin<3 || isempty(method_name)
%     method_name='';
% end
% 
% % outp = method.output_plot;
% outp =output_plot;
% outpRand =output_plotRAND;
% 
% no_filter_entries = numel(outp.data);
% 
% if isempty(output_plotRAND)
%     outTableCell = cell(6,no_filter_entries);
% else
%     outTableCell = cell(10,no_filter_entries);
% end
% 
% % header_string = {'Description' 'Results' 'Results before statistics' 'Results (Distribution)' 'Results before statistics (Distribution)'};
% row_stringTable = {'Description'; 'Results'; 'Results Raw'; 'Results Distribution'; 'Results Raw Distribution'; ...
%     'Results (RAND)'; 'Results Raw (RAND)'; 'Results Distribution (RAND)'; 'Results Raw Distribution (RAND)'; 'Axis'};
% header_stringTable = outp.legendstring;
%     
% for fe = 1:no_filter_entries
%     
%     outTableCell{1,fe} = squeeze(outp.data_string{fe});
%     outTableCell{2,fe} = squeeze(outp.data{fe});
%     outTableCell{3,fe} = squeeze(outp.data_sign{fe});
%     outTableCell{4,fe} = squeeze(outp.data_Median{fe});
%     outTableCell{5,fe} = squeeze(outp.data_signMedian{fe});
%     
% end
% 
% if ~isempty(output_plotRAND)
%     for fe = 1:no_filter_entries
%         outTableCell{6,fe} = squeeze(outpRand.data{fe});
%         outTableCell{7,fe} = squeeze(outpRand.data_sign{fe});
%         outTableCell{8,fe} = squeeze(outpRand.data_Median{fe});
%         outTableCell{9,fe} = squeeze(outpRand.data_signMedian{fe});
%         if  isfield(outp,'ticks')
%             outTableCell{10,fe} = outpRand.XTick;
%         end
%     end
% end
% if isempty(output_plotRAND) && isfield(outp,'ticks')
%     outTableCell{6,fe} = outp.XTick;
% else
%     outTableCell(6,:) = [];
% end
% try
%     %     outTable = cell2table(outTableCell,...
%     %     'RowNames',row_stringTable(1:size(outTableCell,1)),...
%     %     'VariableNames',cellfun(@(str) str(ismember(str,['A':'Z' 'a':'z' '0':'9'])),header_stringTable(:),'UniformOutput',false));
%     if ~verLessThan('matlab', '14.0.1')
%         outTable = cell2table(outTableCell,...
%             'RowNames',row_stringTable(1:size(outTableCell,1)),...
%             'VariableNames',header_stringTable(:));
%         if ~isempty(method_name)
%             outTable.Properties.Description = ['Results for ' method_name];
%         end
%         
%     
%     else
%         outTable = [row_stringTable(1:size(outTableCell,1)) outTableCell];
%         outTable = vertcat([cell(1,1) header_stringTable(:)'],outTable);
%     end
% catch
%     keyboard
% end