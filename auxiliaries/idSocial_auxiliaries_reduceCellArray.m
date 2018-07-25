function output=idSocial_auxiliaries_reduceCellArray(orig_cell,plot_mode,orig_cell4statistics,extra_dim_names,options)
% Each idSocial function creates a cell array ORIG_CELL,
% which is saved in a file at the location given in
% data.('name of the actual function').output.
%
% ORIG_CELL contains the results of all the calculations for
% each group, subset (e.g., day), trial (single trajectory),
% time step, focal, neighbor, and additional fields
% depending on the actual function. Thus, OUTPUT is a cell
% array with dimensions
%
% (#Groups -by- #Subsets -by- #Trials -by- #Individuals
%   -by- #Individuals -by- ...).
%
% NOTE: In each group there may be a different number of
% subsets, and for each subset a different number of trials
% etc. Therefore, dimension sizes here refer to the maximum
% number of subsets in all groups, and analogously for
% trials, individuals etc.
%
% In general, we want to plot and compare different pieces
% of all the data. This is done by the setting PLOT_MODE.
% PLOT_MODE uses string specifiers to distinguish subsets of
% the data:
% 'Group', 'Day', 'Trial', 'Time', 'Focal', 'Neighbor',
%                                           'EdgeX', 'EdgeY'
%
% PLOT_MODE is a structure with the following fields:
%       XAXIS:  Specifies the x-axis:
%                   plot_mode.xaxis = {'Time'}
%               will plot the data versus time.
%               Some functions do not permit changing the
%               axes, for example in case of distributions/
%               histograms, the x-axis is set to 'EdgeX'
%               automatically, corresponding to the bin
%               edges set before running the function (and,
%               analogously, the y-axis is set to 'EdgeY'
%               in case of bivariate histograms/maps).
%
%       YAXIS:  see XAXIS. Only applicable to bivariate
%               histograms/maps.
%
%       FILTER: Specifies/filters data for calculations.
%               Filter is a cell array with as many
%               rows as the number of data subsets you want
%               to compare (e.g., the number of legend
%               entries in a plot). Each row has the form
%
%                    {'Specifier 1' [indices] ...
%                          'Specifier 2' [indices] ...}
%               where 'Specifier' is a string specifier
%               from the list above, and 'indices' determine
%               the data indices which will be used for the
%               calculations. E.g.,
%
%                    plot_mode.filter =
%                           {'Group' 1 'Trial' 1; ...
%                           'Group' 2 'Trial' 3 }
%
%               corresponds to data from Group 1,
%               Trial 1 in the first plot, and Group 2,
%               Trial 3 in the second.
%
%               If the field filter does not exist or is
%               empty, all the data will be used.
%
%       DATA:   Specifies which statistical operation will
%               be applied to which dimension of OUTPUT.
%               Like FILTER, it is a cell array with as many
%               rows as the number of data subsets you want
%               to compare.
%               Each row has the form
%                   {'Stat. op.' 'Specifier 1' ...
%                          'Stat. op.' 'Specifier 2' ...}
%               Statistical operation may be 'Mean',
%               'Median', 'Min', 'Max', 'Hist' or 'Pool'.
%               E.g.,
%
%                    plot_mode.data =
%                           {'Mean' 'Trial'; ...
%                           'Median' 'Trial'}
%
%               will calculate the average over all trials
%               (which are left after filtering) for the
%               first plot, and the median of all trials for
%               the second.
%
%               If the field data does not exist or is
%               empty, all the data will be pooled together
%               and the mean is calculated.
%
%       AUTOCOMBS: Instead of selecting each plot entry
%               seperately by defining various rows in
%               FILTER, AUTOCOMBS can be used to generate
%               automatically all the possible index
%               combinations for a list of string
%               specifiers.
%               E.g., if #Groups = 2 and #Trials = 2
%
%                    plot_mode.autocombs =
%                           {'Group Trial'}
%
%               will generate combinations
%               (Group 1, Trial 1), (Group 1, Trial 2),
%               (Group 2, Trial 1), (Group 2, Trial 2)
%               and plot the corresponding data separately
%               (in the same figure window).
%
%               Note that if FILTER is defined, combinations
%               are only generated from the indices present
%               in FILTER. In this case, FILTER and DATA can
%               only have one row, and the statistical
%               operation(s) in data will be applied to all
%               combinations.
%
%       LEGENDSTRING: Can be used to manually set the legend
%               string. For each row in FILTER/each
%               combination generated by autocombs, there
%               must be one row in LEGENDSTRING:
%                   plot_mode.legendstring =
%                           {'Mean over trials'; ...
%                           'Median of trials'}
%
%       STATISTICS does not apply to all functions, since
%               for many functions STATISTICS is fixed. It
%               determines the statistical operation which
%               is applied directly to the output of the
%               core function. This only makes sense if
%               the output is a numeric vector or a
%               n -by- m cell of numeric vectors (e.g., a 1
%               -by- m cell where each element of the cell
%               corresponds to the bin of a histogram and
%               contains all the values falling into this
%               bin).
%               The statistical operation can be 'Mean',
%               'Median' and 'Hist' (the latter simply
%               counts the number of values).
%               E.g.,
%
%                    plot_mode.statistics =
%                           'Mean';
%
%        SUBPLOTS

%% TODO:
% - If we want to skip every other index of dimension 'time', we have to take into
% account not only that time is not the x-axis, but neither appears in the filter. Which
% makes it more complicated because we have to generate a 'orig_cell' for each filter
% entry.

% plot_mode.display_mode='extern';

prctl=5;
output=[];
% A huge try-catch in order to make sure that calculations are not lost
% only because of a stupid problem with reduce_cell
% try

% If orig_cell is a string giving the path to the data

if ischar(orig_cell)
    orig_cell = idSocial_auxiliaries_loadResults(orig_cell);
end
% keyboard
if nargin<3 || isempty(orig_cell4statistics)
    orig_cell4statistics=orig_cell;
end
if ischar(orig_cell4statistics)
    orig_cell4statistics = idSocial_auxiliaries_loadResults(orig_cell4statistics);
end
disp('Preparing data for presentation...')

if nargin<5 || isempty(options)
    options = [];
end

if isfield(plot_mode,'statistics') && ischar(plot_mode.statistics)
    plot_mode.statistics={plot_mode.statistics};
end

%% Get some information about the input cell etc.
% if ~isempty(orig_cell4statistics)
%     if any(cellfun(@(x) isa(x,'cell'),orig_cell4statistics(:)))
%         cell4statistics_max_elems = max(cellfun(@(x) numel(x),[orig_cell4statistics{:}]));
%     elseif any(cellfun(@(x) isnumeric(x) & ~isempty(x),orig_cell4statistics(:)))
%         temp = cellfun(@(x) numel(x),orig_cell4statistics);
%         cell4statistics_max_elems = max(temp(:));
%     end
% end

if isfield(plot_mode,'parts_and_slices')
    parts_and_slices = plot_mode.parts_and_slices;
else
    parts_and_slices = 0;
end


datacell_size=cellfun(@(x) size(x),orig_cell,'UniformOutput',false);
cell_ndims=ndims(orig_cell);
cell_size=size(orig_cell);

if cell_ndims < 5 % Only one animal in trajectories
    cell_ndims = 6;
    cell_size = [cell_size 1 1];
end

cell_ndims_total=cell_ndims+length(max(vertcat(datacell_size{cellfun(@(x) any(x>1),datacell_size)})));
if all(cellfun(@(x) isa(x,'double'),orig_cell(cellfun(@(x) ~isempty(x),orig_cell)))) && ~strcmp(plot_mode.display_mode,'hist')
    data_size=1;%max(vertcat(datacell_size{:}));
elseif all(cellfun(@(x) isa(x,'double'),orig_cell(cellfun(@(x) ~isempty(x),orig_cell)))) && strcmp(plot_mode.display_mode,'hist')
    data_size=max(vertcat(datacell_size{:}));
elseif all(cellfun(@(x) isa(x,'cell'),orig_cell(cellfun(@(x) ~isempty(x),orig_cell))))
    data_size=max(vertcat(datacell_size{:}));
end
% % For 'rank' an additional dimension will be created
% if strcmpi(plot_mode.statistics{1},'rank') 
%     data_size = [cell4statistics_max_elems 1];
%     cell_ndims_total = cell_ndims_total + 1;
%     plot_mode.xaxis = {'EdgeX'};
% end
% data_type=class(orig_cell{1});
% if iscell(orig_cell{1})
%     data_type_incell=class(orig_cell{1}{1});
%     data_maxsize_incell=max(max(cellfun(@(x) length(x),vertcat(orig_cell{:}))));
% else
%     data_type_incell=[];
%     data_maxsize_incell=1;
% end
cell_size_total=[cell_size data_size];

dim_names_orig={'Group'; 'Subset'; 'Trial'; 'Time'; 'Focal'; 'Neighbor';  'EdgeX'; 'EdgeY'; 'EdgeZ'};
dim_names_short={'Gr'; 'Subset'; 'Tr'; 'Tm'; 'Foc'; 'Nb';  'EdX'; 'EdY'; 'EdZ'};
legend_identifiers=cellstr(['A':'Z']');
% legend_identifiers=arrayfun(@(x) {x},strcat('A':'Z',': '));
for k=1:400
    legend_identifiers=vertcat(legend_identifiers,cellstr(strcat(['A':'Z']',[num2str(k) ' '])));
end
legend_identifiers=legend_identifiers';
dim_names=dim_names_orig(1:cell_ndims_total);
dim_names_short=dim_names_short(1:cell_ndims_total);
timeidx=find(strcmpi(dim_names,'time'));
groupidx=find(strcmpi(dim_names,'group'));
trialidx=find(strcmpi(dim_names,'trial'));
focidx=find(strcmpi(dim_names,'focal'));
nbidx=find(strcmpi(dim_names,'neighbor'));
dayidx=find(strcmpi(dim_names,'subset'));


if isfield(plot_mode,'reorder') && ~isempty(plot_mode.reorder)
    act_dim = find(strcmpi(dim_names_orig,plot_mode.reorder{1}));
    if ~isempty(act_dim)
        
        bl = plot_mode.reorder{2}(:);
        bl = bl(~isnan(bl));
        
        
        bledges = plot_mode.reorder{3};
        [ct id]=histc(bl,bledges);
        
        new_size_dim = numel(plot_mode.reorder{3})-1;
        
        new_cell_size = cell_size;
        new_cell_size(act_dim) = new_size_dim;
        new_cell_size(trialidx) = max(ct);
        
        new_cell_size_total = cell_size_total;
        new_cell_size_total(act_dim) = new_size_dim;
        new_cell_size_total(trialidx) = max(ct);
        new_cell = cell(new_cell_size);
        new_cell4statistics = cell(new_cell_size);
        
        for gr=1:cell_size(1)
            for dy=1:cell_size(2)
                for tr=1:cell_size(3)
                    bl_act = plot_mode.reorder{2}(gr,dy,tr);
                    [~,act_id] = histc(bl_act,bledges);
                    if act_id>0 && ~isempty(act_id)
                        for tm=1:cell_size(4)
                            for foc=1:cell_size(5)
                                for nb=1:cell_size(6)
                                    tr_cnt = find(cellfun(@(x) isempty(x),new_cell(gr,act_id,:,tm,foc,nb)),1,'first');
                                    new_cell(gr,act_id,tr_cnt,tm,foc,nb) = orig_cell(gr,dy,tr,tm,foc,nb);
                                    new_cell4statistics(gr,act_id,tr_cnt,tm,foc,nb) = orig_cell4statistics(gr,dy,tr,tm,foc,nb);
        
                                end
                            end
                        end
                    end
                end
            end
        end
        cell_size = new_cell_size;
        cell_size_total = new_cell_size_total;
        orig_cell = new_cell;
        orig_cell4statistics=new_cell4statistics;
    end
end


if nargin<4 || isempty(extra_dim_names)
    extra_dim_names=dim_names;
end
if nargin<2 || isempty(plot_mode)
    plot_mode.newfigure={[]};
    plot_mode.subplots={[]};
    plot_mode.data={[]};
    plot_mode.statistics={'MEAN'};
    plot_mode.filter=dim_names(max([find(cell_size>1,1) 1]));
    warning([mfilename ': No plot mode specified'],[mfilename ':  No plot mode specified, set automatically.']);
end


if ~isfield(plot_mode,'newfigure') || isempty(plot_mode.newfigure)
    plot_mode.newfigure={[]};
end
if ~isfield(plot_mode,'autocombs') || isempty(plot_mode.autocombs)
    plot_mode.autocombs={[]};
end
if ~isfield(plot_mode,'subplots') || isempty(plot_mode.subplots)
    plot_mode.subplots={[]};
end
if ~isfield(plot_mode,'pool') || isempty(plot_mode.pool)
    plot_mode.pool=[];
end
if ~isfield(plot_mode,'filter') || isempty(plot_mode.filter)
    plot_mode.filter=[];
end
% if ~isfield(plot_mode,'data') || isempty(plot_mode.data)
%     plot_mode.data=[];
% end
if ~isfield(plot_mode,'display_mode') || isempty(plot_mode.display_mode)
    plot_mode.display_mode='plot2d';
end

if ~isfield(plot_mode,'data') || isempty(plot_mode.data)
    %     plot_mode.filter=dim_names(max([find(cell_size>1,1) 1]));
    disp([mfilename ': No data specified, set automatically.']);
    plot_mode.data=[];
end
% Set plot_mode.statistics='POOL', which in case of data_maxsize_incell==1 means: Do nothing.
if (~isfield(plot_mode,'statistics') || isempty(plot_mode.statistics) || isempty([plot_mode.statistics{:}]))
    plot_mode.statistics={'MEAN'};
end
% The following is to change the old format: If
% plot_mode.data is empty and plot_mode.pool_data==true,
% data will be pooled before 'statistics' is applied.
if isfield(plot_mode,'pool_data') && plot_mode.pool_data
    if isempty(plot_mode.data) && ~isempty(plot_mode.statistics)
        plot_mode.data={plot_mode.statistics{1}};
        plot_mode.statistics = {'Pool'};
    end
end


% If filter is given as a single string out of the possible 'dim_names', add a second
% column (indices) with element '[]', which later will be translated to all possible
% indices of that particular dimension:
if all(size(plot_mode.filter)==1) && ischar(plot_mode.filter{1}) && any(strcmp(dim_names,plot_mode.filter{1}))
    plot_mode.filter={plot_mode.filter{1} []};
end
% Similarly, if the data has only two columns, 'MEAN'/'POOL' etc. and the fieldname,
% e.g., 'Trial', add a third column (indices) with element '[]', which later will be
% translated to all possible indices of that particular dimension:
if isfield(plot_mode,'data') && isequal(size(plot_mode.data),[1 2]) && all(cellfun(@(x) ischar(x),plot_mode.data)) && any(strcmp(dim_names,plot_mode.data{1,2}))
    plot_mode.data={plot_mode.data{1} plot_mode.data{2} []};
end
% Similarly for subplots:
if all(size(plot_mode.subplots)==1) && ischar(plot_mode.subplots{1}) && any(strcmp(dim_names,plot_mode.subplots{1}))
    plot_mode.subplots={plot_mode.subplots{1} []};
end
% If display_mode=='hist', xaxis should always be 'EdgeX'
if strcmpi(plot_mode.display_mode,'hist')
    plot_mode.xaxis={'EdgeX'};
end
% If data_size>1, i.e., if the output of the core function is a vector (histogram) or map,
% the x-axis should be 'EdgeX' (in theory, you could collapse the additional dimensions again by
% for example calculate the mean over the whole vector representing the histogram, but what
% sense does it make, after calculating it?); same for y-axis.
if length(data_size)>=1 && data_size(1)>1
    plot_mode.xaxis={'EdgeX'};
end
if length(data_size)>1 && data_size(2)>1
    plot_mode.yaxis={'EdgeY'};
end
% If display_mode is 'plot2d', and therefore no definition of y-axis is needed, set
% yaxis=empty
if (~isfield(plot_mode,'yaxis') || isempty(plot_mode.yaxis)) && ...
        (strcmpi(plot_mode.display_mode,'plot2d') || strcmpi(plot_mode.display_mode,'hist') ||...
        strcmpi(plot_mode.display_mode,'density') ||...
        strcmpi(plot_mode.display_mode,'prob') || strcmpi(plot_mode.display_mode,'extern') )
    plot_mode.yaxis=[];
end

if (~isfield(plot_mode,'data_filter') || isempty(plot_mode.data_filter))
    data_filter = [];
else
    data_filter = plot_mode.data_filter;
end


%% Checking input data

% 2. Pool fields come in couples or more elements

% 3. Are all fields used for each entry (group once, trial once, etc.)
% ----will be checked after the next cell-------

% 4. Indices given for filters do not exceed corresponding array dimension
% ----will be checked after the next cell-------

% 5. Number of indices given in plot_mode.filter equals the number of given
% fields (if indices are not empty)
% ----will be checked after the next cell-------
% 6. Is the x-axis defined?
if ~isfield(plot_mode,'xaxis') || isempty(plot_mode.xaxis)
    error([mfilename ': NoXAxis'],[mfilename ': No x-axis value defined.']);
end
% 7. Is the y-axis defined, if necessary (i.e., only for maps)?
if (~isfield(plot_mode,'yaxis') || isempty(plot_mode.yaxis)) && ...
        ~(strcmpi(plot_mode.display_mode,'plot2d') || strcmpi(plot_mode.display_mode,'hist') ||...
        strcmpi(plot_mode.display_mode,'density') ||...
        strcmpi(plot_mode.display_mode,'prob') || strcmpi(plot_mode.display_mode,'extern'))
    error([mfilename ': NoYAxis'],[mfilename ': No y-axis value defined.']);
end
%% filter: Filter certain indices from given dimensions.

% Dimension names can reappear in data, therefore they are neglected in
% dim_names_control
filter_idx=idSocial_function_wrapper_cellstr2idx(plot_mode.filter,1:2:2*size(dim_names,1),dim_names);

[data_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.data,2:2:3*size(dim_names,1),dim_names);
no_filterentries=max([size(filter_idx,1) size(data_idx,1) 1]);
if size(dim_names_control,1)<no_filterentries
    dim_names_control=repmat(dim_names_control,[no_filterentries 1]);
end
if size(data_idx,1)>no_filterentries
    data_idx=data_idx(1:no_filterentries,:);
end

[autocomb_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.autocombs,1:size(dim_names,1),dim_names, dim_names_control);


[pool_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.pool,1:size(dim_names,1),dim_names, dim_names_control);
% The first dimension in 'pool' is still available, containing all the
% pooled data; therefore return it to dim_names_control as available for further operations:
if ~isempty(plot_mode.pool)
    for lg=1:no_filterentries
        dim_names_control{lg,pool_idx{lg}(1)}=dim_names{pool_idx{lg}(1)};
    end
end


[figure_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.newfigure,1,dim_names,dim_names_control);
figure_idx=figure_idx(1,:);
[subplot_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.subplots,1,dim_names,dim_names_control);
subplot_idx=subplot_idx(1,:);

% if cell_ndims_total>5 || ~isempty(plot_mode.xaxis)
pm_xaxis=repmat(plot_mode.xaxis,[no_filterentries,1]);
[xaxis_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(pm_xaxis,1,dim_names,dim_names_control);
if isempty(xaxis_idx) || isempty(xaxis_idx{1})
    error([mfilename ': WrongXAxis'],[mfilename ': The given xaxis dimension does not exist for this function.']);
end
% end
% if cell_ndims_total>6 || ~isempty(plot_mode.yaxis)
pm_yaxis=repmat(plot_mode.yaxis,[no_filterentries,1]);
[yaxis_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(pm_yaxis,1,dim_names,dim_names_control);
% end



% In case plot_mode.data is empty, the idea is to pool all the superflous
% dimensions (see 'pool_rest' below) and calculate the MEAN over the pooled
% data. To do so, we look for the first not-yet-used field and apply 'Mean'
% to it. 'rest_pool' below will do the rest.
% Fist check if there are fields left:
[~,fst]=max(cellfun(@(x) ~isempty(x),dim_names_control),[],2);

% Firstly, in case plot_mode.data exists, but has only one
% entry specifying the statistical operation for all legend entries
if isfield(plot_mode,'data') && ~isempty(plot_mode.data) && size(plot_mode.data,2)==1 && ~isempty(plot_mode.data{1,1})
    plot_mode.data{1,1}= plot_mode.data{1,1};
    plot_mode.data{1,2}=dim_names_control{1,fst(1)};
    for lg=2:no_filterentries
        plot_mode.data{lg,1}= plot_mode.data{1,1};
        plot_mode.data{lg,2}= plot_mode.data{1,2};
        %             plot_mode.data{lg,3}=[];
    end
    [data_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.data,2:2:3*size(dim_names,1),dim_names,dim_names_control);
    
end


% Secondly, in case plot_mode.data exists, but is empty for only some filter
% entries.
if (isfield(plot_mode,'data') && ~isempty(plot_mode.data) && ~isempty([plot_mode.data{:}]))
    for lg=1:no_filterentries
        if isempty(plot_mode.data{lg,1})
            plot_mode.data{lg,1}='MEAN';
            plot_mode.data{lg,2}=dim_names_control{lg,fst(lg)};
            %             plot_mode.data{lg,3}=[];
            [data_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.data,2:2:3*size(dim_names,1),dim_names,dim_names_control);
            
        end
        
    end
end
% Thirdly, in case plot_mode.data does not exists or is empty for all filter
% entries.
if (~isfield(plot_mode,'data') || isempty(plot_mode.data) || isempty([plot_mode.data{:}]))
    plot_mode.data=cell(no_filterentries,2);
    for lg=1:no_filterentries
        plot_mode.data{lg,1}='MEAN';
        plot_mode.data{lg,2}=dim_names_control{lg,fst(lg)};
        %         plot_mode.data{lg,3}=[];
    end
    [data_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(plot_mode.data,2:2:3*size(dim_names,1),dim_names,dim_names_control);
end

% Finally, if there are dimensions not assigned to any operation
% (data,pool etc.) consider them superflous and merge ('pool') them with
% the first dimension a statistical operation (mean, median, etc.) is
% applied to.
first_statop_dim=vertcat(data_idx{:,2});
% [~,last_statop_idx]=unique(cellfun(@(x) ~isempty(x),data_idx(:,2:3:end)));
last_statop_idx=NaN(no_filterentries,1);
for lg=1:no_filterentries
    nempt=cellfun(@(x) ~isempty(x),data_idx(lg,2:2:end));
    fidx=find(nempt==0);
    if isempty(fidx)
        last_statop_idx(lg)=size(nempt,2);
    else
        last_statop_idx(lg)=fidx(end)-1;
    end
end

pool_rest=cell(no_filterentries,1);
for lg=1:no_filterentries
    %     pool_rest{lg}=[dim_names{data_idx{lg,2}} ' '];
    for k=1:size(dim_names_control(lg,:),2)
        if ~isempty(dim_names_control{lg,k})
            pool_rest{lg}=[pool_rest{lg} dim_names_control{lg,k} ' '];
        end
    end
end
[pool_rest_idx, dim_names_control]=idSocial_function_wrapper_cellstr2idx(pool_rest,1:size(dim_names,1),dim_names, dim_names_control);
for lg=1:no_filterentries
    pool_rest_idx{lg}=[first_statop_dim(lg) pool_rest_idx{lg}];
end

% If there is "normalization" for one legend/filter-entry extend it to all
% of them
if isfield(plot_mode,'normalization') && ~isempty(plot_mode.normalization) 
    if ischar(plot_mode.normalization)
        norm_entry = plot_mode.normalization;
    elseif iscell(plot_mode.normalization)
        norm_entry = plot_mode.normalization{1};
    end
    plot_mode.normalization = cell(no_filterentries,1);
    for lg=1:no_filterentries
        plot_mode.normalization{lg} = norm_entry;
    end
end

% Those basic dimensions (all those which do not depend on size of the function
% output) which have length==1 do not need to be taken care of especially, i.e., it is not a
% problem if there is no corresponding plot_mode entry. Therefore, check them out:
dim_names_control(cell_size==1)=repmat({''},[sum(cell_size==1),1]);

% Are all fields (with dim_size>1) used for each entry (group once, trial once, etc.)
if ~isempty(vertcat(dim_names_control{repmat(cell_size_total>1,[no_filterentries,1])}))
    error([mfilename ': LeftOvers'],[mfilename ': Some items have not been used']);
end
% Check if indices given for filters do exceed corresponding array dimension
if size(filter_idx,2)>1
    for lg=1:no_filterentries
        if  ~isempty(filter_idx) && ~isempty(filter_idx{lg,2})
            for k=1:size(filter_idx{lg,1},2)
                if cell_size_total(filter_idx{lg,1}(k))<filter_idx{lg,2}(k)
                    warning([mfilename ': IndexRange'],[mfilename ': filter indices exceeds available range.'])
                    error([mfilename ': IndexRange'],[mfilename ': filter indices exceeds available range.']);
                end
            end
        end
    end
end
if size(filter_idx,2)>1 && ~all(cellfun(@(x,y) size(x,1)==size(y,1) || isempty(y),filter_idx(:,1),filter_idx(:,2)))
    error([mfilename ': IndexRange'],[mfilename ': Number of filter indices does not equal number of field names.']);
end


%% Shift additional dimensions "from cell to cell" and apply statistics corresponding to plot_mode.statistics


ndms=unique(cellfun(@(x) ndims(x),orig_cell));
if length(ndms)>1
    error([mfilename ': NonHomogenousSize'],[mfilename ': Function outputs differ in size.']);
end


if size(plot_mode.statistics,1)==1 && no_filterentries>1
    plot_mode.statistics=repmat(plot_mode.statistics,[no_filterentries,1]);
end
extended_cell=cell(no_filterentries,1);
extended_cell_numberDataPoints=cell(no_filterentries,1);
% The following array will hold the median in case the core function output is a histogram.
% Later the median will be used to test for differences in distributons
extended_cell_histstatistics=cell(no_filterentries,1);


funcOutType=unique(cellfun(@(x) class(x),orig_cell(cellfun(@(x) ~isempty(x),orig_cell)),'UniformOutput',false)); % class of output of core function
if isa(orig_cell,'cell') && ~strcmp(plot_mode.display_mode,'hist') %%&& ...
    %         ~strfind(plot_mode.statistics{1},'hist')
    funcOutType='cell';
end
notempt=cellfun(@(x) ~isempty(x),orig_cell);
orig_cell_elem=orig_cell(notempt);

if isa(orig_cell_elem{1},'cell')
    funcOutType='cellofcells';
end

%% Complete filter strings
filter_complete=cell(no_filterentries,cell_ndims*2);
for row=1:no_filterentries
    filter_complete(row,1:2:end)=num2cell((1:cell_ndims)');
    for col=1:2:cell_ndims*2
        if ~isempty(filter_idx) && ~isempty(vertcat(filter_idx{row,1:2:end}))
            idx=find([filter_idx{row,1:2:end}]==filter_complete{row,col});
        else
            idx=[];
        end
        
        if ~isempty(idx) %&& col~=2*timeidx-1 % special treatment for dimension 'time'
            %             filter_complete{row,col+1}=filter_idx{row,idx*2-1+1};
            
            if  length(filter_idx{row,idx*2-1+1})==max(filter_idx{row,idx*2-1+1})-min(filter_idx{row,idx*2-1+1})+1 && ...
                    length(filter_idx{row,idx*2-1+1})>1
                filter_complete{row,col+1}=[num2str(min(filter_idx{row,idx*2-1+1})) '-' num2str(max(filter_idx{row,idx*2-1+1}))];
            else
                filter_complete{row,col+1}=strrep(num2str(filter_idx{row,idx*2-1+1}),'  ',',');
                
            end
            
        elseif isempty(idx) && xaxis_idx{row}~=filter_complete{row,col} && col~=2*timeidx-1 % speeial tratment for dimension 'time'
            
            if  cell_size(filter_complete{row,col})>2
                filter_complete{row,col+1}=['1-' num2str(cell_size(filter_complete{row,col}))];
            else
                ttt2=textscan(num2str(1:cell_size(filter_complete{row,col})),'%s','Delimiter',' ');
                ttt2=ttt2{1}(cellfun(@(x) ~isempty(x),ttt2{1}))';
                ttt2=strcat(ttt2,',');
                ttt2=[ttt2{:}];
                ttt2=ttt2(1:end-1);
                filter_complete{row,col+1}=ttt2;
                
            end
        elseif isempty(idx) && xaxis_idx{row}~=filter_complete{row,col} && col==2*timeidx-1 % special treatment for dimension 'time'
            filter_complete{row,col+1}=strrep(num2str(2:2:cell_size(filter_complete{row,col})-1),'  ',',');
            
            if  length(2:2:cell_size(filter_complete{row,col})-1)>2
                filter_complete{row,col+1}=['2-' num2str(cell_size(filter_complete{row,col})-1)];
            else
                ttt2=textscan(num2str(2:2:cell_size(filter_complete{row,col})),'%s','Delimiter',' ');
                ttt2=ttt2{1}(cellfun(@(x) ~isempty(x),ttt2{1}))';
                ttt2=strcat(ttt2,',');
                ttt2=[ttt2{:}];
                ttt2=ttt2(1:end-1);
                filter_complete{row,col+1}=ttt2;
                
            end
            
        end
        
        %         end
        filter_complete{row,col}=dim_names_short{ceil(col/2)};
        if isempty(filter_complete{row,col+1})
            filter_complete{row,col}=[];
        end
    end
end
% filter_complete(:,find(all(cellfun(@(x) isempty(x), filter_complete),1)))=[];


%%

% Default (if no filter is given): Take all indices.
gr_idces=1:size(orig_cell,groupidx);
day_idces=1:size(orig_cell,dayidx);
tr_idces=1:size(orig_cell,trialidx);
tm_idces=1:size(orig_cell,timeidx);
foc_idces=1:size(orig_cell,focidx);
nb_idces=1:size(orig_cell,nbidx);

if ~isempty(options) 
    if isfield(options,'filter_focal_list') && ~isempty(options.filter_focal_list) && all(isfinite(options.filter_focal_list))
        foc_idces = options.filter_focal_list;
    end
    if isfield(options,'filter_neighbor_list') && ~isempty(options.filter_neighbor_list) && all(isfinite(options.filter_neighbor_list))
        nb_idces = options.filter_neighbor_list;
    end
end

% if ~isempty(data_filter) && iscell(data_filter) && size(data_filter,1)==1
%     % data_filter contains id-string - index-vector pairs, e.g.,
%     % data_filter = {'Trial' [1 2 5] 'Focal' 1:2}; from all possible
%     % combinations of groups, trials etc, filter only trials 1 2 5 and
%     % focals 1 and 2
%     select_combs = idSocial_auxiliaries_allcombs(gr_idces,day_idces,tr_idces,tm_idces,foc_idces,nb_idces);
%     for k = 1:2:size(data_filter,2)
%         dimidx = find(strcmpi(dim_names,data_filter{1,k}));
%         if ~isempty(dimidx)
%             select_combs_idx = ...
%                 any(repmat(select_combs(dimidx,:)',[1 numel(data_filter{1,k+1})])==repmat(data_filter{1,k+1},[size(select_combs,2) 1]),2);
%             select_combs = select_combs(:,select_combs_idx);
%         end
%     end
%     data_filter = select_combs';
if ~isempty(data_filter) && iscell(data_filter) %&& size(data_filter,1)>1
    % Each row contains seperate filters. Add possible combinations of
    % 'free' indices
    all_combs = idSocial_auxiliaries_allcombs(gr_idces,day_idces,tr_idces,tm_idces,foc_idces,nb_idces)';
    final_combs = NaN(size(all_combs));
    act_idx = 1;
    for k1 = 1:size(data_filter,1)
%         idx_in_filter = ismember(dim_names,data_filter(k,1:2:end));
%         data_filter_idces = data_filter(k,2:2:end);
%         data_filter_idces
%         
%         
        select_combs = idSocial_auxiliaries_allcombs(gr_idces,day_idces,tr_idces,tm_idces,foc_idces,nb_idces)';
        for k = 1:2:size(data_filter,2)
            dimidx = find(strcmpi(dim_names,data_filter{k1,k}));
            if ~isempty(dimidx)
                select_combs_idx = ...
                    any(repmat(select_combs(dimidx,:)',[1 numel(data_filter{k1,k+1})])==repmat(data_filter{k1,k+1},[size(select_combs,2) 1]),2);
                select_combs = select_combs(:,select_combs_idx);
            end
        end
        if ~isempty(select_combs)
            final_combs(:,act_idx : act_idx + size(select_combs,2)-1) = select_combs;
            act_idx = act_idx + size(select_combs,2) ;
        end
        
        
    end
    data_filter = final_combs(:, ~all(isnan(final_combs),1))';

end

for pidx=1:no_filterentries
    
    %Extract array dimensions according to indices in data_idx:
    if ~isempty(filter_idx) && size(filter_idx,1)>=pidx && ~isempty(filter_idx(pidx,:))
        for lc=1:2:size(filter_idx(pidx,:),2)
            if ~isempty(filter_idx{pidx,lc}) && ~isempty(filter_idx{pidx,lc+1})
                switch filter_idx{pidx,lc}
                    %Group
                    case groupidx
                        gr_idces=filter_idx{pidx,lc+1};
                        if max(gr_idces)>size(orig_cell,groupidx)
                            error([mfilename ': WrongIndex'],[mfilename ': Given group index exceeds array size.']);
                        end
                        cell_size_total(groupidx)=length(filter_idx{pidx,lc+1});
                        cell_size(groupidx)=length(filter_idx{pidx,lc+1});
                        %Subset/day
                    case dayidx
                        day_idces=filter_idx{pidx,lc+1};
                        if max(day_idces)>size(orig_cell,dayidx)
                            warning([mfilename ': WrongIndex'],[mfilename ': Given subset index exceeds array size.']);
                            day_idces=day_idces(day_idces<=size(orig_cell,dayidx));
                            
                        end
                        cell_size_total(dayidx)=length(filter_idx{pidx,lc+1});
                        cell_size(dayidx)=length(filter_idx{pidx,lc+1});
                        %Trial
                    case trialidx
                        tr_idces=filter_idx{pidx,lc+1};
                        if max(tr_idces)>size(orig_cell,trialidx)
                            error([mfilename ': WrongIndex'],[mfilename ': Given trial index exceeds array size.']);
                        end
                        cell_size_total(trialidx)=length(filter_idx{pidx,lc+1});
                        cell_size(trialidx)=length(filter_idx{pidx,lc+1});
                        %Time
                    case timeidx
                        tm_idces=filter_idx{pidx,lc+1};
                        if max(tm_idces)>size(orig_cell,timeidx)
                            error([mfilename ': WrongIndex'],[mfilename ': Given time index exceeds array size.']);
                        end
                        cell_size_total(timeidx)=length(filter_idx{pidx,lc+1});
                        cell_size(timeidx)=length(filter_idx{pidx,lc+1});
                        %Focal
                    case focidx
                        foc_idces=filter_idx{pidx,lc+1};
                        if max(foc_idces)>size(orig_cell,focidx)
                            error([mfilename ': WrongIndex'],[mfilename ': Given focal index exceeds array size.']);
                        end
                        cell_size_total(focidx)=length(filter_idx{pidx,lc+1});
                        cell_size(focidx)=length(filter_idx{pidx,lc+1});
                        %Neighbor
                    case nbidx
                        nb_idces=filter_idx{pidx,lc+1};
                        if max(nb_idces)>size(orig_cell,nbidx)
                            error([mfilename ': WrongIndex'],[mfilename ': Given neighbor index exceeds array size.']);
                        end
                        cell_size_total(nbidx)=length(filter_idx{pidx,lc+1});
                        cell_size(nbidx)=length(filter_idx{pidx,lc+1});
                        
                end
            end
        end
    end
    extended_cell{pidx}=cell(cell_size_total);
    extended_cell_numberDataPoints{pidx}=NaN(cell_size_total);
    extended_cell_histstatistics{pidx}=cell(cell_size);
    
    fprintf('%s %s %s \n','Apply', plot_mode.statistics{pidx}, 'to')
    
    % Apply data_filter (filter: one legend entry per entry in plot_mode.filter. data_filter: simply omit data which is not in data_filter)
    for gr=gr_idces
        for dy=day_idces
            for tr=tr_idces
                for tm=tm_idces
                    for foc=foc_idces
                        for nb=nb_idces
                            S.type='()';
                            S.subs=vertcat({find(gr==gr_idces);...
                                find(dy==day_idces);...
                                find(tr==tr_idces);...
                                find(tm==tm_idces);...
                                find(foc==foc_idces);...
                                find(nb==nb_idces);},repmat({':'},ndms,1));
                            if mod(sum([find(gr==gr_idces) find(dy==day_idces) find(tr==tr_idces) find(tm==tm_idces) find(foc==foc_idces) find(nb==nb_idces)]),50) == 0
                                fprintf('Passing through  (%d,%d,%d,%d,%d,%d) \n ',gr,dy,tr,tm,foc,nb)
                            end
                            
                            try
%                                 isempty(data_filter) || ismember([gr dy tr tm foc nb],data_filter,'rows');
                            
                            if isempty(data_filter) || ismember([gr dy tr tm foc nb],data_filter,'rows')
%                                 if dy==13 && tr==2 && tm ==2; keyboard; end
                                if isa(orig_cell4statistics{gr,dy,tr,tm,foc,nb},'cell')
                                    cell4statistics = vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb}{:});
                                elseif isa(orig_cell4statistics{gr,dy,tr,tm,foc,nb},'double')
                                    cell4statistics = orig_cell4statistics{gr,dy,tr,tm,foc,nb};
                                end
                                
                                if strcmp(funcOutType,'cellofcells')
                                    
                                    try
                                        if ~isempty(orig_cell{gr,dy,tr,tm,foc,nb}) &&  ~isempty(vertcat(orig_cell{gr,dy,tr,tm,foc,nb}{:}))%~all(all(cellfun(@(x) isempty(x),orig_cell{gr,dy,tr,tm,foc,nb})))
                                            
                                            
                                            
                                            switch plot_mode.statistics{pidx}
                                                case {'Mean','mean','MEAN'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) nanmean(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(cell4statistics));
                                                case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) nanmean(abs(x)),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(abs(cell4statistics)));
                                                    
                                                case {'Circ_Mean','circ_mean','CIRC_MEAN'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_circ_mean(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_mean(cell4statistics));
                                                case {'Circ_Var','circ_var','CIRC_VAR'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_circ_var(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_var(cell4statistics));
                                                    
                                                case {'Std','std','STD'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) nanstd(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanstd(cell4statistics));
                                                case {'Var','var','VAR'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) nanvar(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanvar(cell4statistics));
                                                    
                                                case {'Median','median','MEDIAN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) nanmedian(x),orig_cell4statistics{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
%                                                 case {'Rank', 'rank', 'RANK'}
%                                                     if length(plot_mode.statistics(pidx,:)) == 1
%                                                         no_datapoints = 1;
%                                                     else
%                                                         no_datapoints = plot_mode.statistics{pidx,2};
%                                                     end
%                                                     sortvals = sort(cell4statistics,'descend');
%                                                     sortvals_format4all = NaN(1,cell4statistics_max_elems);
%                                                     sortvals_format4all(1:numel(sortvals)) = sortvals;
%                                                     S.type='()';
%                                                     extended_cell{pidx}=subsasgn(extended_cell{pidx},S,num2cell(sortvals_format4all));
%                                                     S.type='{}';
%                                                     extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(sortvals_format4all(1:no_datapoints)));
%                                          
                                                case {'Positive_Ratio','positive_ratio','POSITIVE_RATIO'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_positive_ratio(x,[],0),orig_cell4statistics{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    
                                                    
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_positive_ratio(cell4statistics,[],0));

                                                case {'Angle_Ratio','angle_ratio','ANGLE_RATIO'}
                                                    if length(plot_mode.statistics(pidx,:)) == 1
                                                        ratio_angle = pi/2;
                                                    else
                                                        ratio_angle = plot_mode.statistics{pidx,2};
                                                    end

                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_angle_ratio(x,[],ratio_angle),orig_cell4statistics{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_angle_ratio(cell4statistics,[],ratio_angle));

%                                                     extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_angle_ratio(cell4statistics,[],ratio_angle));
                                                    
                                                case {'Mode', 'mode', 'MODE'}
                                                    if isfield(plot_mode,'mode_edges') && ~isempty(plot_mode.mode_edges) && ...
                                                            isfield(plot_mode,'mode_bandwith') && ~isempty(plot_mode.mode_bandwith)
                                                        extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_distribution_mode(x,plot_mode.mode_edges,plot_mode.mode_bandwith),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                        S.type='{}';
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_distribution_mode(cell4statistics,plot_mode.mode_edges,plot_mode.mode_bandwith));
                                                    else
                                                        error([mfilename ': plot_mode.mode_edges and plot_mode.mode_bandwith for calculation of distribution and mode are undefined'])
                                                    end
                                                case {'Mode_Height', 'mode_height', 'MODE_HEIGHT'}
                                                    if isfield(plot_mode,'mode_edges') && ~isempty(plot_mode.mode_edges) && ...
                                                            isfield(plot_mode,'mode_bandwith') && ~isempty(plot_mode.mode_bandwith)
                                                        extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_distribution_modeHeight(x,plot_mode.mode_edges,plot_mode.mode_bandwith),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                        S.type='{}';
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_distribution_modeHeight(cell4statistics,plot_mode.mode_edges,plot_mode.mode_bandwith));
                                                    else
                                                        error([mfilename ': plot_mode.mode_edges and plot_mode.mode_bandwith for calculation of distribution and mode are undefined'])
                                                    end
                                                case {'Sum','sum','SUM'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                    
                                                case {'Min','min','MIN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) min(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,min(cell4statistics));
                                                case {'Max','max','MAX'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) max(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,max(cell4statistics));
                                                case {'Pool','pool','POOL'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,orig_cell{gr,dy,tr,tm,foc,nb});
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,cell4statistics);
                                                case {'Hist','hist','HIST'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                case {'Hist+Angle_Ratio','hist+angle_ratio','HIST+ANGLE_RATIO'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_angle_ratio(cell4statistics,[],pi/2));
                                                    
                                                case {'Hist+Mean','hist+mean','HIST+MEAN'}
                                                     S.type='()'; % RH 2017-03-29. This is suddenly necessary, but has not been all the centuries before. It kind of make sense, but should be set already. 
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(cell4statistics));
                                                case {'Hist+Median','hist+median','HIST+MEDIAN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                    
                                                case {'Hist+Mode','hist+mode','HIST+MODE'}
                                                    %                                         keyboard
                                                    %                                                 if dy==13 && tr==2; keyboard; end
                                                    hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    %                                                 hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    S.type='()';
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    
                                                    S.type='{}';
                                                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
%                                                         [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
                                                        [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,plot_mode.edges(locs));
                                                    end
                                                case {'Hist+1stMode','hist+1stMode','HIST+1STMODE','hist+1stMode'}
                                                    S.type='()';
                                                    %                                                 if dy==13 && tr==2; keyboard; end
                                                    if ~isfield(plot_mode,'dist_method') ||  strcmpi(plot_mode.dist_method,'hist')
                                                        hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    else % for ksdensity
                                                        hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    end
                                                    
                                                    %                                                 hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    
                                                    S.type='{}';
                                                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
%                                                         [~,locs] = findpeaks(hitemp);
                                                        [~,locs] = idSocial_auxiliaries_findpeaks(hitemp);
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,plot_mode.edges(min(locs)));
                                                    end
                                                    
                                                case {'Hist+LogNormalMode','hist+lognormalmode','HIST+LOGNORMALMODE'}
                                                     hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    %                                                 hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    S.type='()';
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    
                                                    S.type='{}';
                                                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                                                        parmhat = lognfit(cell4statistics(~isnan(cell4statistics)));
                                                        logf = lognpdf(plot_mode.edges,parmhat(1),parmhat(2));
                                                        [~,mxIdx] = max(logf);
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,plot_mode.edges(mxIdx));
                                                    end


                                                case {'Hist+Circvar','hist+circvar','HIST+CIRCVAR'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_var(cell4statistics));
                                                case {'Hist+Var','hist+var','HIST+VAR'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanvar(cell4statistics));
                                                case {'RankDist+Mean','rankdist+mean','RANKDIST+MEAN'}
                                                    if length(plot_mode.statistics(pidx,:)) == 1
                                                        no_datapoints = 5;
                                                    else
                                                        no_datapoints = plot_mode.statistics{pidx,2};
                                                    end
                                                    S.type='()';
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));

                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean([orig_cell{gr,dy,tr,tm,foc,nb}{1:no_datapoints}]));
                                                 case {'RankMean','rankmean','RANKMEAN'}
                                                    if length(plot_mode.statistics(pidx,:)) == 1
                                                        no_datapoints = 5;
                                                    else
                                                        no_datapoints = plot_mode.statistics{pidx,2};
                                                    end
                                                    S.type='()';
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S, cellfun(@(x) nanmean(x(1:min(no_datapoints,numel(x)))),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(cell4statistics(1:min(no_datapoints,numel(cell4statistics)))));
                                                    
                                            end
                                            
                                            if ~isempty(orig_cell4statistics{gr,dy,tr,tm,foc,nb})
                                                if iscell(orig_cell4statistics{gr,dy,tr,tm,foc,nb})
                                                    S.type='()';
                                                    extended_cell_numberDataPoints{pidx} = subsasgn(extended_cell_numberDataPoints{pidx},S,cellfun(@(x) sum(~isnan(x)),orig_cell4statistics{gr,dy,tr,tm,foc,nb}));
                                                    
                                                elseif isnumeric(orig_cell4statistics{gr,dy,tr,tm,foc,nb}) % For hists?!!!
                                                    S.type='()';
                                                    extended_cell_numberDataPoints{pidx} = subsasgn(extended_cell_numberDataPoints{pidx},S,cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb}));
                                            
                                                end
                                            end
                                            % Normalization for
                                            % histograms/distributions
                                            clear no_frames no_frames_effective edgArea edgAreaRadial
                                            if isfield(plot_mode,'normalization') && ~isempty(plot_mode.statistics{pidx}) && ...
                                                    ( ...
                                                    (~isempty(strfind(lower(plot_mode.statistics{pidx}),'hist')) && strcmpi(plot_mode.display_mode,'hist')) || ...
                                                    strcmpi(plot_mode.display_mode,'MapPolar') ||...
                                                    strcmpi(plot_mode.display_mode,'Map') ...
                                                    ) && ...
                                                     ~strcmpi(plot_mode.statistics{pidx},'pool')
                                                 
                                                 if strcmpi(plot_mode.display_mode,'hist') %&& strcmpi(plot_mode.normalization{pidx},'density')
                                                     if isfield(plot_mode,'edges')
                                                        edgArea = plot_mode.edges(2)-plot_mode.edges(1);
                                                        edgAreaRadial = pi*plot_mode.edges(2:end).^2 - pi*plot_mode.edges(1:end-1).^2;
                                                     else
                                                         edgArea = 1;
                                                     end
                                                     if all(iscell(orig_cell4statistics{gr,dy,tr,tm,foc,nb}))
                                                        no_frames = nansum([orig_cell4statistics{gr,dy,tr,tm,foc,nb}{:}]);
                                                     else
                                                        no_frames = sum(~isnan(orig_cell4statistics{gr,dy,tr,tm,foc,nb}));%numel(orig_cell4statistics{gr,dy,tr,tm,foc,nb});
                                                     end
                                                     
%                                                      no_frames_effective = nansum([extended_cell{pidx}{gr,dy,tr,tm,foc,nb,:}]);%sum(~isnan(orig_cell4statistics{gr,dy,tr,tm,foc,nb}));%
                                                     Sref = S; Sref.type = '{}';
                                                     extcell_temp = subsref(extended_cell{pidx},Sref);
                                                     no_frames_effective = nansum([extcell_temp]);%sum(~isnan(orig_cell4statistics{gr,dy,tr,tm,foc,nb}));%

                                                     clear extcell_temp
                                                 elseif strcmpi(plot_mode.display_mode,'Map') %&& strcmpi(plot_mode.normalization{pidx},'density')
                                                     no_frames = nansum([orig_cell4statistics{gr,dy,tr,tm,foc,nb}{:}]);
                                                     no_frames_effective = nansum([orig_cell4statistics{gr,dy,tr,tm,foc,nb}{:}]);
                                                     edgArea = (plot_mode.edges{1}(2)-plot_mode.edges{1}(1))*(plot_mode.edges{2}(2)-plot_mode.edges{2}(1));
                                                 elseif strcmpi(plot_mode.display_mode,'MapPolar')% && strcmpi(plot_mode.normalization{pidx},'density')
                                                     no_frames = nansum([orig_cell4statistics{gr,dy,tr,tm,foc,nb}{:}]);
                                                     no_frames_effective = nansum([orig_cell4statistics{gr,dy,tr,tm,foc,nb}{:}]);
                                                     ar = pi*plot_mode.edges{1}(2:end).^2 - pi*plot_mode.edges{1}(1:end-1).^2;
                                                     ft =     (plot_mode.edges{2}(2)-plot_mode.edges{2}(1))/(2*pi);
                                                     edgArea  = repmat(ar .* ft,[size(plot_mode.edges{2},2)-1,1]);
                                                     
                                                 end
                                                 
                                                Sref = S; Sref.type = '()';
                                                
                                                switch lower(plot_mode.normalization{pidx})
                                                    
                                                    
                                                    case 'none'
                                                    case 'no_frames'
                                                        extcell_temp = subsref(extended_cell{pidx},Sref);
                                                        extended_cell{pidx}=subsasgn(extended_cell{pidx},Sref,cellfun(@(x) x./no_frames,extcell_temp,'UniformOutput',false));

%                                                         extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = cellfun(@(x) x./no_frames,extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:),'UniformOutput',false);
                                                    case 'no_frames_effective'
                                                        extcell_temp = subsref(extended_cell{pidx},Sref);
                                                        extended_cell{pidx}=subsasgn(extended_cell{pidx},Sref,cellfun(@(x) x./no_frames_effective,extcell_temp,'UniformOutput',false));

%                                                         extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = cellfun(@(x) x./no_frames_effective,extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:),'UniformOutput',false);
                                                    case 'density'
                                                        extcell_temp = subsref(extended_cell{pidx},Sref);
                                                        if ~strcmpi(plot_mode.display_mode,'MapPolar')
%                                                             extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = cellfun(@(x) x./no_frames./edgArea,extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:),'UniformOutput',false);
                                                            extended_cell{pidx}=subsasgn(extended_cell{pidx},Sref,cellfun(@(x) x./no_frames./edgArea,extcell_temp,'UniformOutput',false));

                                                        else
                                                            extended_cell{pidx} = subsasgn(extended_cell{pidx},Sref,num2cell(cell2mat(squeeze(extcell_temp))./edgArea./no_frames));

                                                            % extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = num2cell(cell2mat(squeeze(extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:)))./edgArea./no_frames);
                                                        end
                                                    case 'density_effective'
                                                        extcell_temp = subsref(extended_cell{pidx},Sref);
                                                        if ~strcmpi(plot_mode.display_mode,'MapPolar')
                                                            extended_cell{pidx}=subsasgn(extended_cell{pidx},Sref,cellfun(@(x) x./no_frames_effective./edgArea,extcell_temp,'UniformOutput',false));
                                                            
                                                            % extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = cellfun(@(x) x./no_frames_effective./edgArea,extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:),'UniformOutput',false);
                                                        else
                                                             extended_cell{pidx} = subsasgn(extended_cell{pidx},Sref,num2cell(cell2mat(squeeze(extcell_temp))./edgArea./no_frames_effective));

%                                                             extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = num2cell(cell2mat(squeeze(extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:)))./edgArea./no_frames_effective);
                                                        end
                                                    case 'radial_density'
                                                        extcell_temp = subsref(extended_cell{pidx},Sref);
                                                         extended_cell{pidx} = subsasgn(extended_cell{pidx},Sref,num2cell(cell2mat(squeeze(extcell_temp))./edgAreaRadial'./no_frames));

%                                                         extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = num2cell(cell2mat(squeeze(extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:)))./edgAreaRadial'./no_frames);
                                                    case 'radial_density_effective'
                                                        extcell_temp = subsref(extended_cell{pidx},Sref);
                                                        extended_cell{pidx} = subsasgn(extended_cell{pidx},Sref,num2cell(cell2mat(squeeze(extcell_temp))./edgAreaRadial'./no_frames_effective));

%                                                         extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:) = num2cell(cell2mat(squeeze(extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:)))./edgAreaRadial'./no_frames_effective);

                                                end
                                                clear extcell_temp
                                            end
                                            
                                            
                                        end
                                    catch
                                        keyboard
                                    end
                                elseif strcmp(funcOutType,'double') && any(size(orig_cell{1})>1) % Core function output is numeric vector
                                    S.type='{}';
                                    try
                                        if ~isempty(orig_cell{gr,dy,tr,tm,foc,nb})
                                            switch plot_mode.statistics{pidx}
                                                case {'Mean','mean','MEAN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanmean(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(cell4statistics));
                                                case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanmean(abs(orig_cell{gr,dy,tr,tm,foc,nb})));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(abs(cell4statistics)));
                                                    
                                                case {'Circ_Mean','circ_mean','CIRC_MEAN'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_circ_mean(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_mean(cell4statistics));
                                                case {'Circ_Var','circ_var','CIRC_VAR'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_circ_var(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_var(cell4statistics));
                                                    
                                                case {'Std','std','STD'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanstd(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanstd(cell4statistics));
                                                case {'Var','var','VAR'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanstd(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanstd(cell4statistics));
                                                    
                                                case {'Median','median','MEDIAN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanmedian(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                case {'Positive_Ratio','positive_ratio','POSITIVE_RATIO'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_positive_ratio(orig_cell{gr,dy,tr,tm,foc,nb},[],0));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_positive_ratio(cell4statistics,[],0));
                                                case {'Angle_Ratio','angle_ratio','ANGLE_RATIO'}
                                                    if length(plot_mode.statistics(pidx,:)) == 1
                                                        ratio_angle = pi/2;
                                                    else
                                                        ratio_angle = plot_mode.statistics{pidx,2};
                                                    end
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_angle_ratio(orig_cell{gr,dy,tr,tm,foc,nb},[],ratio_angle));
                                                    S.type='{}';
                                                    
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_angle_ratio(cell4statistics,[],ratio_angle));
                                                case {'Sum','sum','SUM'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_auxiliaries_nansum(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_auxiliaries_nansum(cell4statistics));
                                                    
                                                case {'Min','min','MIN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,min(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,min(cell4statistics));
                                                case {'Max','max','MAX'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,max(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,max(cell4statistics));
                                                    
                                                case {'Pool','pool','POOL'}
                                                    S.type='()';
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,num2cell(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,cell4statistics);
                                                case {'Hist','hist','HIST'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,sum(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                case {'Hist+Mode','hist+mode','HIST+MODE'}
                                                    if ~isfield(plot_mode,'dist_method') ||  strcmpi(plot_mode.dist_method,'hist')
                                                        %                                                         hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                        hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    else % for ksdensity
                                                        %                                                         hitemp=cellfun(@(x) sum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                        hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    end
                                                    %                                                     hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));

                                                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                                                        %                                                         [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
                                                        [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
                                                        
                                                        
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,plot_mode.edges(locs));
                                                    end
                                                case {'Hist+1stMode','hist+1stMode','HIST+1STMODE','hist+1stmode'}
                                                    if ~isfield(plot_mode,'dist_method') ||  strcmpi(plot_mode.dist_method,'hist')
                                                        %                                                         hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                        hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    else % for ksdensity
                                                        %                                                         hitemp=cellfun(@(x) sum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                        hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    end
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    
                                                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
%                                                         [~,locs] = findpeaks(hitemp);
                                                        [~,locs] = idSocial_auxiliaries_findpeaks(hitemp);

                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,plot_mode.edges(min(locs)));
                                                    end
                                                case {'Hist+Circvar','hist+circvar','HIST+CIRCVAR'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_auxiliaries_nansum(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_var(cell4statistics));
                                                case {'Hist+Median','hist+median','HIST+MEDIAN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_auxiliaries_nansum(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                    
                                                case {'Hist+Var','hist+var','HIST+VAR'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_auxiliaries_nansum(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanvar(cell4statistics));
                                                    
                                                    
                                            end
                                            
                                            
                                        end
                                    catch
                                        keyboard
                                    end
                                elseif strcmp(funcOutType,'cell')  % Core function output is cell vector
                                    S.type='{}';
                                    try
                                        if ~isempty(orig_cell{gr,dy,tr,tm,foc,nb})
                                            switch plot_mode.statistics{pidx}
                                                case {'Mean','mean','MEAN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanmean(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(vertcat(orig_cell{gr,dy,tr,tm,foc,nb})));
                                                case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanmean(abs(orig_cell{gr,dy,tr,tm,foc,nb})));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmean(abs(vertcat(orig_cell{gr,dy,tr,tm,foc,nb}))));
                                                    
                                                case {'Positive_Ratio','positive_ratio','POSITIVE_RATIO'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_positive_ratio(orig_cell{gr,dy,tr,tm,foc,nb},[],0));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_positive_ratio(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb}),[],0));
                                                case {'Angle_Ratio','angle_ratio','ANGLE_RATIO'}
                                                    if length(plot_mode.statistics(pidx,:)) == 1
                                                        ratio_angle = pi/2;
                                                    else
                                                        ratio_angle = plot_mode.statistics{pidx,2};
                                                    end
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_angle_ratio(orig_cell{gr,dy,tr,tm,foc,nb},[],ratio_angle));
                                                    S.type='{}';
                                                    
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_angle_ratio(cell4statistics,[],ratio_angle));
                                                    
                                                    
                                                case {'Circ_Mean','circ_mean','CIRC_MEAN'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_circ_mean(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_mean(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                case {'Circ_Var','circ_var','CIRC_VAR'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_circ_var(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_var(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                    
                                                case {'Std','std','STD'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanstd(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanstd(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                case {'Var','var','VAR'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanvar(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanvar(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                    
                                                case {'Median','median','MEDIAN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,nanmedian(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                    
                                                case {'Min','min','MIN'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,min(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,min(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                case {'Max','max','MAX'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,max(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,max(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                case {'Sum','sum','SUM'}
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,idSocial_auxiliaries_nansum(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_auxiliaries_nansum(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb})));
                                                    %
                                                case {'Mode','mode','MODE'}
                                                    if ~all(isnan(orig_cell{gr,dy,tr,tm,foc,nb}))
                                                        
                                                        S.type='{}';
                                                        
                                                        if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                                                            if ~isfield(plot_mode,'dist_bandwidth') || isempty(plot_mode.dist_bandwidth)
                                                                warning([mfilename ': Kernel density width for mode calculation set = 1'])
                                                                plot_mode.dist_bandwidth = 1;
                                                            end
                                                            if ~isfield(plot_mode,'dist_method') || isempty(plot_mode.dist_method)
                                                                warning([mfilename ': Kernel density method for mode calculation set to "ksdensity_epanechnikov"'])
                                                                plot_mode.dist_method = 'ksdensity_epanechnikov';
                                                            end
                                                            
                                                            md = idSocial_distribution_mode(orig_cell{gr,dy,tr,tm,foc,nb},plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method);
                                                            
                                                            extended_cell{pidx}=subsasgn(extended_cell{pidx},S,md);
                                                            
                                                            md = idSocial_distribution_mode(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb}),plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method);
                                                            
                                                            extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,md);
                                                        end
                                                        
                                                        
                                                    end
                                                    
                                                case {'Mode_Height','mode_height','MODE_HEIGHT'}
                                                    if ~all(isnan(orig_cell{gr,dy,tr,tm,foc,nb}))
                                                        
                                                        S.type='{}';
                                                        
                                                        if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                                                            if ~isfield(plot_mode,'dist_bandwidth') || isempty(plot_mode.dist_bandwidth)
                                                                warning([mfilename ': Kernel density width for mode calculation set = 1'])
                                                                plot_mode.dist_bandwidth = 1;
                                                            end
                                                            if ~isfield(plot_mode,'dist_method') || isempty(plot_mode.dist_method)
                                                                warning([mfilename ': Kernel density method for mode calculation set to "ksdensity_epanechnikov"'])
                                                                plot_mode.dist_method = 'ksdensity_epanechnikov';
                                                            end
                                                            
                                                            md = idSocial_distribution_modeHeight(orig_cell{gr,dy,tr,tm,foc,nb},plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method);
                                                            
                                                            extended_cell{pidx}=subsasgn(extended_cell{pidx},S,md);
                                                            
                                                            md = idSocial_distribution_modeHeight(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb}),plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method);
                                                            
                                                            extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,md);
                                                        end
                                                        
                                                        
                                                    end
                                                case {'1stMode','1stmode','1STMODE'}
                                                   
                                                    if ~all(isnan(orig_cell{gr,dy,tr,tm,foc,nb}))
                                                        
                                                        S.type='{}';
%                                                          keyboard
                                                        if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                                                            if ~isfield(plot_mode,'dist_bandwidth') || isempty(plot_mode.dist_bandwidth)
                                                                warning([mfilename ': Kernel density width for mode calculation set = 1'])
                                                                plot_mode.dist_bandwidth = 1;
                                                            end
                                                            if ~isfield(plot_mode,'dist_method') || isempty(plot_mode.dist_method)
                                                                warning([mfilename ': Kernel density method for mode calculation set to "ksdensity_epanechnikov"'])
                                                                plot_mode.dist_method = 'ksdensity_epanechnikov';
                                                            end
                                                            
                                                            md = idSocial_distribution_mode(orig_cell{gr,dy,tr,tm,foc,nb},plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method,'first');
                                                            
%                                                             For tests:
%                                                             [hi,prob] = idSocial_auxiliaries_kerneldensity(orig_cell{gr,dy,tr,tm,foc,nb},plot_mode.edges,plot_mode.dist_bandwidth,plot_mode.dist_method);
                                                            
                                                            extended_cell{pidx}=subsasgn(extended_cell{pidx},S,md);
                                                            
                                                            md = idSocial_distribution_mode(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb}),plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method,'first');
                                                            
                                                            extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,md);
                                                        end
                                                        
                                                        
                                                    end
                                                    %%
                                                case {'1stMode_Height','1stmode_height','1STMODE_HEIGHT'}
                                                    if ~all(isnan(orig_cell{gr,dy,tr,tm,foc,nb}))
                                                        
                                                        S.type='{}';
                                                        
                                                        if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                                                            if ~isfield(plot_mode,'dist_bandwidth') || isempty(plot_mode.dist_bandwidth)
                                                                warning([mfilename ': Kernel density width for mode calculation set = 1'])
                                                                plot_mode.dist_bandwidth = 1;
                                                            end
                                                            if ~isfield(plot_mode,'dist_method') || isempty(plot_mode.dist_method)
                                                                warning([mfilename ': Kernel density method for mode calculation set to "ksdensity_epanechnikov"'])
                                                                plot_mode.dist_method = 'ksdensity_epanechnikov';
                                                            end
                                                            
                                                            md = idSocial_distribution_modeHeight(orig_cell{gr,dy,tr,tm,foc,nb},plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method,'first');
                                                            
                                                            extended_cell{pidx}=subsasgn(extended_cell{pidx},S,md);
                                                            
                                                            md = idSocial_distribution_modeHeight(vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb}),plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method,'first');
                                                            
                                                            extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,md);
                                                        end
                                                        
                                                        
                                                    end
                                                case {'Pool','pool','POOL'}
                                                    %                                     extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:)=ttt;%orig_cell{gr,dy,tr,tm,foc,nb};
                                                    S.type='()';
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,orig_cell(gr,dy,tr,tm,foc,nb));
                                                    S.type='{}';
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,vertcat(orig_cell4statistics{gr,dy,tr,tm,foc,nb}));
                                                    
                                                case {'Hist','hist','HIST'}
                                                    %                                                 extended_cell{pidx}=subsasgn(extended_cell{pidx},S,sum(~isnan(orig_cell{gr,dy,tr,tm,foc,nb})));
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,sum(orig_cell{gr,dy,tr,tm,foc,nb}));
                                                    
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                case {'Hist+Mode','hist+mode','HIST+MODE'}
                                                    if ~isfield(plot_mode,'dist_method') ||  strcmpi(plot_mode.dist_method,'hist')
                                                        %                                                         hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                        hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    else % for ksdensity
                                                        hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    end
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    
                                                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
%                                                         [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
                                                        [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
                                                        
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,plot_mode.edges(locs));
                                                    end
                                                case {'Hist+1stMode','hist+1stMode','HIST+1STMODE','hist+1stmode'}
                                                    if ~isfield(plot_mode,'dist_method') ||  strcmpi(plot_mode.dist_method,'hist')
%                                                         hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                        hitemp=cellfun(@(x) sum(~isnan(x)),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    else % for ksdensity
%                                                         hitemp=cellfun(@(x) sum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                        hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    end
%                                                     hitemp=cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb});
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    
                                                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
%                                                         [~,locs] = findpeaks(hitemp);
                                                        [~,locs] = idSocial_auxiliaries_findpeaks(hitemp);
                                                        
                                                        extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,plot_mode.edges(min(locs)));
                                                    end
                                                case {'Hist+Circvar','hist+circvar','HIST+CIRCVAR'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,idSocial_circ_var(cell4statistics));
                                                case {'Hist+Median','hist+median','HIST+MEDIAN'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanmedian(cell4statistics));
                                                case {'Hist+Var','hist+var','HIST+VAR'}
                                                    
                                                    extended_cell{pidx}=subsasgn(extended_cell{pidx},S,cellfun(@(x) idSocial_auxiliaries_nansum(x),orig_cell{gr,dy,tr,tm,foc,nb},'UniformOutput',false));
                                                    extended_cell_histstatistics{pidx}=subsasgn(extended_cell_histstatistics{pidx},S,nanvar(cell4statistics));
                                                    
                                                    
                                            end
                                            
                                            
                                        end
                                    catch
                                        keyboard
                                    end
                                else % if core output is a single numerical value: convert to cell so to be able to treat extended cell the same for all cases later on.
                                    try
                                        if ~isempty(orig_cell{gr,dy,tr,tm,foc,nb})
                                            extended_cell{pidx}(gr,dy,tr,tm,foc,nb,:,:,:)=num2cell(squeeze(orig_cell{gr,dy,tr,tm,foc,nb}));
                                        end
                                    catch
                                        keyboard
                                    end
                                end
                            end
                            catch
                                keyboard
                            end
                        end
                    end
                end
            end
        end
    end
    fprintf('\n')
end

% CORRECT:
% clear orig_cell

% Treat slices (parts/time from resulting from slicing big trajectories)
if parts_and_slices
    disp('Merging slices ...')
    slice_tic = tic;
    extended_cell_sliced = cell(extended_cell);
    extended_cell_sliced_histstats = cell(extended_cell);
    for pidx=1:no_filterentries
        fprintf('%s %d','Filter entry: ',pidx)
        extended_cell_sliced_size = size(extended_cell{pidx}); 
        extended_cell_sliced_size_histstats = size(extended_cell_histstatistics{pidx}); 
        extended_cell_sliced_size(timeidx) = 1;
        extended_cell_sliced_size_histstats(timeidx) = 1;
        extended_cell_sliced{pidx} = cell(extended_cell_sliced_size);
        extended_cell_sliced_histstats{pidx} = cell(extended_cell_sliced_size_histstats);
        for gr=gr_idces
%             fprintf('%d,',gr)
            for dy=day_idces
%                 fprintf('%d,',dy)
                for tr=tr_idces
%                    fprintf('%d,',tr)
                        for foc=foc_idces
%                             fprintf('%d,',foc)
                            for nb=nb_idces
                                fprintf('(%d,%d,%d,%d,%d): ',gr,dy,tr,foc,nb)
                                S.type='()';
                                S.subs=vertcat({find(gr==gr_idces);...
                                    find(dy==day_idces);...
                                    find(tr==tr_idces);...
                                    ':';...
                                    find(foc==foc_idces);...
                                    find(nb==nb_idces);},repmat({':'},ndms,1)); 
                                
%                                 cellfun(@(x) all(isempty(x),2),extended_cell{pidx}(gr,dy,tr,:,foc,nb,:,:))
                                parts_allEmpty = squeeze(all(cellfun(@(x) isempty(x), extended_cell{pidx}(gr,dy,tr,:,foc,nb,:)),7));
                                
                                ndp = extended_cell_numberDataPoints{pidx}(gr,dy,tr,~parts_allEmpty,foc,nb,:,:);
                                
                                mm = cell2mat(extended_cell{pidx}(gr,dy,tr,~parts_allEmpty,foc,nb,:,:));
                                mm2 = cell2mat(extended_cell_histstatistics{pidx}(gr,dy,tr,~parts_allEmpty,foc,nb,:,:));
                                
                                
                                
                                if ~any(size(mm)==0)
                                    siz = ones(1,ndims(mm));
                                    siz(4) = sum(~parts_allEmpty);
                                    weights = ndp./repmat(nansum(ndp,4),siz);
                                    
                                    try
                                        weights2 = nansum(ndp(1,1,1,:,1,1,:),7)./ ...
                                            repmat(nansum(nansum(ndp(1,1,1,:,1,1,:),7),4),siz);
                                    catch
                                        keyboard
                                    end
                                    try
                                        if any([strfind(lower(plot_mode.statistics{pidx}),'hist'),  ...
                                                strfind(lower(plot_mode.statistics{pidx}),'sum')])
                                            fprintf('%s %d %s \n','Sum over ', size(extended_cell{pidx}(gr,dy,tr,:,foc,nb,:),4),' slices.')
                                            extended_cell_sliced{pidx}(gr,dy,tr,1,foc,nb,:,:,:,:,:) = num2cell(idSocial_auxiliaries_nansum(mm,4));
                                            
                                        else
                                            extended_cell_sliced{pidx}(gr,dy,tr,1,foc,nb,:,:,:,:,:) = num2cell(nansum(mm.*weights,4)./nansum(weights,4));
                                            fprintf('%s %d %s \n','Weighted sum over ', size(extended_cell{pidx}(gr,dy,tr,:,foc,nb,:),4),' slices. ')
                                        end
                                    catch
                                        keyboard
                                    end
                                    clear mm
                                    extended_cell_sliced_histstats{pidx}(gr,dy,tr,1,foc,nb,:,:,:,:,:) = num2cell(nansum(mm2.*weights2 ,4)./nansum(weights2 ,4));
                                end
                            end
                        end
                 
                end
            end
        end
        fprintf('\n')
    end
    extended_cell = extended_cell_sliced;
    clear extended_cell_sliced
    disp(['Merging done. This took ' num2str(toc(slice_tic)) 's.'])
end
    
%% Pool data (merge cell dimensions)
pool_cell=cell(no_filterentries,1);
pool_cell_histstatistics=cell(no_filterentries,1);
dim_names_act=cell(no_filterentries,size(dim_names,1));
data_cell_size=cell(no_filterentries,1);

% Special treatment for dimension 'time': Since it counts datapoints into overlapping bins (spacing 2),
% we eliminate the first (t=0 to binwidth/2), the last (t=end-binwidth/2 to
% end) and every other bin in case time IS NOT xaxis. We do this in order to avoid data duplication.
% EXCEPTION: When time is filtered manually in 'filter_idx' or when parts =
% time come from slicing big trajecotries. 
if ~isempty(pool_idx)
    for pidx=1:no_filterentries
        if isempty(filter_idx) || ~any(vertcat(filter_idx{pidx,1:2:end})==timeidx) 
            if  ~isempty(pool_idx{pidx})
                if any(pool_idx{pidx}==timeidx) && parts_and_slices ~= 1
                    extsize=size(extended_cell{pidx},timeidx);
                    S.type='()';
                    S.subs=repmat({':'},ndims(extended_cell{pidx}),1);
                    S.subs{timeidx}=2:2:extsize-1;
                    extended_cell{pidx}=subsref(extended_cell{pidx},S);
                    extended_cell_histstatistics{pidx}=subsref(extended_cell_histstatistics{pidx},S);
                end
            end
        end
    end
end


for lg=1:no_filterentries
    data_cell_size{lg}=size(extended_cell{lg});
    %     for k=3:3:size(data_idx,2)
    %         if ~isempty(data_idx{lg,k})
    %             data_cell_size{lg}(data_idx{lg,k-1})=size(data_idx{lg,k},2);
    %         end
    %     end
end
for pidx=1:no_filterentries
    if ~isempty(pool_idx)
        
        if  ~isempty(pool_idx{pidx})
            
            
            cell_permute=permute(extended_cell{pidx}, [setdiff(1:ndims(extended_cell{pidx}),pool_idx{pidx}),pool_idx{pidx}]);
            cell_permute_histstatistics=permute(extended_cell_histstatistics{pidx},...
                [setdiff(1:ndims(extended_cell_histstatistics{pidx}),pool_idx{pidx}),pool_idx{pidx}]);
            
            %             dim_names_permute=dim_names([setdiff(1:size(dim_names,1),pool_idx{pidx}),pool_idx{pidx}]);
            S.type='()';
            S.subs=repmat({':'},ndims(extended_cell{pidx})-size(pool_idx{pidx},2)+1,1);
            cell_pool=subsref(cell_permute,S);
            cell_pool_histstatistics=subsref(cell_permute_histstatistics,S);
            
            % Insert some dimensions of length 1 where before have been dimensions which not longer exist:
            reshsiz=size(extended_cell{pidx});
            reshsiz=[reshsiz prod(reshsiz(pool_idx{pidx}))];
            reshsiz(pool_idx{pidx})=ones(1,size(pool_idx{pidx},2));
            cell_pool=reshape(cell_pool,reshsiz);
            reshsiz=size(extended_cell_histstatistics{pidx});
            reshsiz=[reshsiz prod(reshsiz(pool_idx{pidx}))];
            reshsiz(pool_idx{pidx})=ones(1,size(pool_idx{pidx},2));
            cell_pool_histstatistics=reshape(cell_pool_histstatistics,reshsiz);
            % Shift the newly merged dimension from the end of the array to the
            % place/the dimension on which the first statistical operation will
            % be applied a bit further down:
            permidx=[1:min(pool_idx{pidx})-1 ndims(cell_pool) min(pool_idx{pidx})+1:ndims(cell_pool)-1 min(pool_idx{pidx})];
            pool_cell{pidx}=permute(cell_pool,permidx);
            permidx=[1:first_statop_dim(pidx)-1 ndims(cell_pool_histstatistics) first_statop_dim(pidx)+1:ndims(cell_pool_histstatistics)-1 first_statop_dim(pidx)];
            pool_cell_histstatistics{pidx}=permute(cell_pool_histstatistics,permidx);
            
            
            
        end
        
    else
        pool_cell{pidx}=extended_cell{pidx};
        pool_cell_histstatistics{pidx}=extended_cell_histstatistics{pidx};
        %     for pidx=1:no_filterentries
        %         dim_names_act(pidx,1:size(dim_names,1))=dim_names';
        %     end
    end
end
clear cell_permute extended_cell

%% Pool rest data=dimensions which have not been assigned manually (merge cell dimensions)

% Special treatment for dimension 'time': Since it counts datapoints into overlapping bins (spacing 2),
% we eliminate the first (t=0 to binwidth/2), the last (t=end-binwidth/2 to
% end) and every other bin in case time IS NOT xaxis. We do this in order to avoid data duplication.
% EXCEPTION: When time is filtered manually in 'filter_idx'
for pidx=1:no_filterentries
    pool_rest_idx{pidx}=pool_rest_idx{pidx}(pool_rest_idx{pidx}<=ndims(pool_cell{pidx}));
    if any(pool_rest_idx{pidx}==timeidx) && ...
            ~(~isempty(pool_idx) && any(pool_idx{pidx}==timeidx)) && parts_and_slices ~= 1
        if isempty(filter_idx) || ~any(vertcat(filter_idx{pidx,1:2:end})==timeidx)
            if ~isempty(pool_rest_idx{pidx}) && size(pool_rest_idx{pidx},2)>1
                
                poolsize=size(pool_cell{pidx},timeidx);
                S.type='()';
                S.subs=repmat({':'},ndims(pool_cell{pidx}),1);
                S.subs{timeidx}=2:2:poolsize-1;
                pool_cell{pidx}=subsref(pool_cell{pidx},S);
                pool_cell_histstatistics{pidx}=subsref(pool_cell_histstatistics{pidx},S);
                
                
                data_cell_size{pidx}(timeidx)=size(pool_cell{pidx},timeidx);
            end
        end
    end
end

pool_rest_cell=cell(no_filterentries,1);
pool_rest_cell_histstatistics=cell(no_filterentries,1);
for pidx=1:no_filterentries
    if ~isempty(pool_rest_idx{pidx}) && size(pool_rest_idx{pidx},2)>1
        
        % Shift the dimensions we want to merge to the end of the cell:
        cell_permute=permute(pool_cell{pidx}, [setdiff(1:ndims(pool_cell{pidx}),pool_rest_idx{pidx}),pool_rest_idx{pidx}]);
        cell_permute_histstatistics=permute(pool_cell_histstatistics{pidx},[setdiff(1:ndims(pool_cell{pidx}),pool_rest_idx{pidx}),pool_rest_idx{pidx}]);
        % Merge the <size(pool_rest_idx{pidx},2)> last dimensions using
        % (:,:,..,:)
        S.type='()';
        S.subs=repmat({':'},ndims(pool_cell{pidx})-size(pool_rest_idx{pidx},2)+1,1);
        
        cell_pool_rest=subsref(cell_permute,S);
        cell_pool_rest_histstatistics=subsref(cell_permute_histstatistics,S);
        
        % Insert some dimensions of length 1 where before have been dimensions which no longer exist:
        reshsiz=size(pool_cell{pidx});
        reshsiz=[reshsiz prod(reshsiz(pool_rest_idx{pidx}(pool_rest_idx{pidx}<=length(reshsiz))))];
        reshsiz(pool_rest_idx{pidx}(pool_rest_idx{pidx}<length(reshsiz)))=ones(1,size(pool_rest_idx{pidx}(pool_rest_idx{pidx}<length(reshsiz)),2));
        cell_pool_rest=reshape(cell_pool_rest,reshsiz);
        
        reshsiz_median=size(pool_cell_histstatistics{pidx});
        reshsiz_median=[reshsiz_median prod(reshsiz_median(pool_rest_idx{pidx}(pool_rest_idx{pidx}<=length(reshsiz_median))))];
        reshsiz_median(pool_rest_idx{pidx}(pool_rest_idx{pidx}<length(reshsiz_median)))=ones(1,size(pool_rest_idx{pidx}(pool_rest_idx{pidx}<length(reshsiz_median)),2));
        cell_pool_rest_histstatistics=reshape(cell_pool_rest_histstatistics,reshsiz_median);
        % Shift the newly merged dimension from the end of the array to the
        % place/the dimension on which the first statistical operation will
        % be applied a bit further down:
        permidx=[1:first_statop_dim(pidx)-1 length(reshsiz) first_statop_dim(pidx)+1:length(reshsiz)-1 first_statop_dim(pidx)];
        pool_rest_cell{pidx}=permute(cell_pool_rest,permidx);
        permidx=[1:first_statop_dim(pidx)-1 length(reshsiz_median) first_statop_dim(pidx)+1:length(reshsiz_median)-1 first_statop_dim(pidx)];
        pool_rest_cell_histstatistics{pidx}=permute(cell_pool_rest_histstatistics,permidx);
        
        
    else
        
        pool_rest_cell{pidx}=pool_cell{pidx};
        pool_rest_cell_histstatistics{pidx}=pool_cell_histstatistics{pidx};
    end
end
clear cell_permute pool_cell

% Make dimension names for pooled data
for pidx=1:no_filterentries
    dim_names_act(pidx,:)=dim_names_short;
    % If pool and rest pool have dimensions in common and will end up in the same index:
    if ~isempty(pool_idx) && ~isempty(pool_rest_idx) && isempty(intersect(pool_rest_idx{pidx},pool_idx{pidx}))
        % For pool_rest:
        if ~isempty(pool_rest_idx{pidx}) && size(pool_rest_idx{pidx},2)>1
            sort_pool_rest_idx=sort(pool_rest_idx{pidx});
            
            
            dim_name_merged=dim_names_short(sort_pool_rest_idx);
            dim_name_merged2=[];
            for k=1:size(dim_name_merged,2)
                dim_name_merged2=[dim_name_merged2 dim_name_merged{k} ','];
                %             keyboard
            end
            siz_pool_cell=size(pool_cell{pidx});
            dim_name_merged2=dim_name_merged2(1:end-1);
            dim_name_merged2=[dim_name_merged2 '{' num2str(prod(siz_pool_cell(sort_pool_rest_idx))) '*'...
                num2str(max(cellfun(@(x) length(x),vertcat(pool_rest_cell{pidx}(:))))) '}'];
            dim_names_act{pidx,first_statop_dim(pidx)}=dim_name_merged2;
            
            if  ~isempty(setdiff(pool_rest_idx{pidx},first_statop_dim(pidx)))
                dim_names_act(pidx,setdiff(pool_rest_idx{pidx},first_statop_dim(pidx)))=repmat({''},[1,size(setdiff(pool_rest_idx{pidx},first_statop_dim(pidx)),2)]);
            end
            
        end
        % For pool:
        
        if  ~isempty(pool_idx) && ~isempty(pool_idx{pidx})
            sort_pool_idx=sort(pool_idx{pidx});
            sort_pool_idx=sort_pool_idx(sort_pool_idx<=length(filter_complete(pidx,:))/2);
            
            dim_name_merged=dim_names_short(sort_pool_idx);
            dim_name_merged2=[];
            for k=1:size(dim_name_merged,1)
                dim_name_merged2=[dim_name_merged2 dim_name_merged{k} '[' filter_complete{pidx,sort_pool_idx(k)*2-1+1} '],'];
                %             keyboard
            end
            dim_name_merged2=dim_name_merged2(1:end-1);
            dim_name_merged2=[dim_name_merged2 '{' num2str(prod(data_cell_size{pidx}(sort_pool_idx))) '}'];
            dim_names_act{pidx,min(pool_idx{pidx})}=dim_name_merged2;
            
            if  ~isempty(setdiff(pool_idx{pidx},min(pool_idx{pidx})))
                dim_names_act(pidx,setdiff(pool_idx{pidx},min(pool_idx{pidx})))=repmat({''},[1,size(setdiff(pool_idx{pidx},min(pool_idx{pidx})),2)]);
            end
            
        end
        
    else
        
        if  ~isempty(pool_idx) && ~isempty(pool_idx{pidx}) || ~isempty(pool_rest_idx) && ~isempty(pool_rest_idx{pidx})
            if isempty(pool_idx)
                merge_idces=pool_rest_idx{pidx};
            elseif isempty(pool_rest_idx)
                merge_idces=pool_idx{pidx};
            else
                merge_idces=unique([pool_idx{pidx} pool_rest_idx{pidx}]);
            end
            sort_pool_idx=sort(merge_idces);
            sort_pool_idx=sort_pool_idx(sort_pool_idx<=length(filter_complete(pidx,:))/2);
            dim_name_merged=dim_names_short(sort_pool_idx);
            dim_name_merged2=[];
            for k=1:size(dim_name_merged,1)
                %                 dim_name_merged2=[dim_name_merged2 dim_name_merged{k} '[' num2str(data_cell_size{pidx}(sort_pool_idx(k))) '],'];
                dim_name_merged2=[dim_name_merged2 dim_name_merged{k} '[' filter_complete{pidx,sort_pool_idx(k)*2-1+1} '],'];
                
                
                %             keyboard
            end
            dim_name_merged2=dim_name_merged2(1:end-1);
            dim_name_merged2=[dim_name_merged2 '{' num2str(prod(data_cell_size{pidx}(sort_pool_idx))) '*'...
                num2str(max(cellfun(@(x) length(x),vertcat(pool_rest_cell{pidx}(:))))) '}'];
            dim_names_act{pidx,first_statop_dim(pidx)}=dim_name_merged2;
            if  ~isempty(setdiff(merge_idces,first_statop_dim(pidx)))
                dim_names_act(pidx,setdiff(merge_idces,first_statop_dim(pidx)))=repmat({''},[1,size(setdiff(merge_idces,first_statop_dim(pidx)),2)]);
            end
            
            
            if  ~isempty(setdiff(merge_idces,first_statop_dim(pidx)))
                dim_names_act(pidx,setdiff(merge_idces,first_statop_dim(pidx)))=repmat({''},[1,size(setdiff(merge_idces,first_statop_dim(pidx)),2)]);
            end
        else
            
            dim_names_act(pidx,1:size(dim_names,1))=strcat(dim_names_short,'[', num2str(data_cell_size{pidx}'),']')';
        end
        
    end
    
end

%% Apply Mean/Median/etc. corresponding to data_idx
% data_cell=cell(no_filterentries,1);
data_array=cell(no_filterentries,1);
data_Median_array=cell(no_filterentries,1);
data_dev_array=cell(no_filterentries,1);
data_sign_array=cell(no_filterentries,1);
data_signMedian_array=cell(no_filterentries,1);
statistics_on_idx=NaN(no_filterentries,1);
statistics_type = cell(no_filterentries,1);
% First step: Take the first index in data_idx, apply statistics to that
% dimension, and transform the output cell into a double array.Then perform the rest of the
% operations on this array.

% Special treatment for dimension 'time': Since it counts datapoints into overlapping bins (spacing 2),
% we eliminate the first (t=0 to binwidth/2), the last (t=end-binwidth/2 to
% end) and every other bin in case time IS NOT xaxis. We do this in order to avoid data duplication.
% EXCEPTION: When time is filtered manually in 'filter_idx'
for pidx=1:no_filterentries
    if isempty(filter_idx) || ~any(vertcat(filter_idx{pidx,1:2:end})==timeidx)
        if ~isempty(data_idx{pidx}) && ~all(cellfun(@(x) isempty(x),data_idx(pidx,:)))
            
            if any(data_idx{pidx,2}==timeidx) && ...
                    ~(~isempty(pool_rest_idx) && any(pool_rest_idx{pidx}==timeidx)) && ...
                    ~(~isempty(pool_idx) && any(pool_idx{pidx}==timeidx))
                S.type='()';
                S.subs=repmat({':'},ndims(pool_rest_cell{pidx}),1);
                S.subs{timeidx}=2:2:size(pool_rest_cell{pidx},timeidx)-1;
                pool_rest_cell{pidx}=subsref(pool_rest_cell{pidx},S);
                pool_rest_cell_histstatistics{pidx}=subsref(pool_rest_cell_histstatistics{pidx},S);
            end
        end
    end
end

for pidx=1:no_filterentries
    if ~isempty(data_idx{pidx}) && ~all(cellfun(@(x) isempty(x),data_idx(pidx,:)))
        
        
        
        data_sign_cell_size=[size(pool_rest_cell{pidx}) ones(1,5)];
        
        data_Median_cell_size=[size(pool_rest_cell_histstatistics{pidx}) ones(1,5)];
        data_signMedian_cell_size=data_Median_cell_size;
        
        data_cell_size=data_sign_cell_size;
        data_cell_size(data_idx{pidx,2})=1;
        data_Median_cell_size(data_idx{pidx,2})=1;
        
        
        data_array{pidx}=NaN(data_cell_size);
        data_Median_array{pidx}=NaN(data_Median_cell_size);
        data_dev_array{pidx}=[];
        data_sign_array{pidx}=[];
        data_signMedian_array{pidx}=[];
        
        
        
        idxcombis=vectorsize2indexcombinations(data_cell_size);
        
        % Merge/pool the cells; these may have size 1 (when the core function returns a
        % double, or >1 when it returns values sorted into bins (=cell))
        merged_cells=cell(size(idxcombis,1),1);
        
        for idcnt=1:size(idxcombis,1)
            S.subs=num2cell(idxcombis(idcnt,:)');
            S.type='()';
            S.subs(data_idx{pidx,2})={':'};
            
            A=squeeze(subsref(pool_rest_cell{pidx},S));
            A(cellfun(@(x) isempty(x),A))={NaN};
            merged_cells{idcnt}=vertcat(A{:});
        end
        idxcombis_Median=vectorsize2indexcombinations(data_Median_cell_size);
        merged_Median_cells=cell(size(idxcombis_Median,1),1);
        for idcnt=1:size(idxcombis_Median,1)
            S.subs=num2cell(idxcombis_Median(idcnt,:)');
            S.type='()';
            S.subs(data_idx{pidx,2})={':'};
            B=squeeze(subsref(pool_rest_cell_histstatistics{pidx},S));
            B(cellfun(@(x) isempty(x),B))={NaN};
            merged_Median_cells{idcnt}=vertcat(B{:});
        end
        
        pool_data_cell_size=cellfun(@(x) length(x),merged_cells);
        max_pool_data_cell_size=max(pool_data_cell_size);
        pool_data_cell_size_Median=cellfun(@(x) length(x),merged_Median_cells);
        max_pool_data_cell_size_Median=max(pool_data_cell_size_Median);
        
        if last_statop_idx(pidx)==1
            data_sign_cell_size(data_idx{pidx,2})=max_pool_data_cell_size;
            data_signMedian_cell_size(data_idx{pidx,2})=max_pool_data_cell_size_Median;
            data_sign_array{pidx}=NaN(data_sign_cell_size);
            data_signMedian_array{pidx}=NaN(data_signMedian_cell_size);
        else
            data_sign_array{pidx}=NaN(data_cell_size);
            data_signMedian_array{pidx}=NaN(data_Median_cell_size);
        end
        
        % Convert into numeric array
        for idcnt=1:size(idxcombis,1)
            % The data_sign_array misses
            % the last statistical
            % operation:
            % If the first operation is
            % also the last AND the core
            % data has not simply been
            % pooled together (which could mean that pool_rest_cell still cotains cells
            % of various length) or if what I said in () is not true, store the
            % values from which the
            % mean/median/etc is calculated
            % above.
            % In the histogram case, we have allready calculated the medians
            if last_statop_idx(pidx)==1
                
                %                 tic
                %                 S.subs=num2cell(idxcombis(idcnt,:)');
                %                 S.type='()';
                %                 S.subs(data_idx{pidx,2})={1:pool_data_cell_size(idcnt)};
                %                                 data_sign_array{pidx} = ...
                %                     subsasgn(data_sign_array{pidx},S,merged_cells{idcnt});
                %                 toc %Elapsed time is 0.246172 seconds!!!!!!!!!!!!!!!!!
                %
                % There are some speed issues here. I tried various methods, so far
                % the less elegant is the fastest in my case: simple loop
                %                 tic
                %                 for nd=1:pool_data_cell_size(idcnt)
                %                     idx_vector=ones(1,12);
                %                     idx_vector(1,1:numel(idxcombis(idcnt,:)))=idxcombis(idcnt,:);
                %                     idx_vector(data_idx{pidx,2})=nd;
                %                     data_sign_array{pidx}(idx_vector(1),...
                %                         idx_vector(2),...
                %                         idx_vector(3),...
                %                         idx_vector(4),...
                %                         idx_vector(5),...
                %                         idx_vector(6),...
                %                         idx_vector(7),...
                %                         idx_vector(8),...
                %                         idx_vector(9),...
                %                         idx_vector(10),...
                %                         idx_vector(11),...
                %                         idx_vector(12))=...
                %                         merged_cells{idcnt}(nd);
                %
                %                 end
                %               for nd=1:pool_data_cell_size(idcnt)
                idx_vector=num2cell(ones(1,12));
                idx_vector(1,1:numel(idxcombis(idcnt,:)))=num2cell(idxcombis(idcnt,:));
                idx_vector(data_idx{pidx,2})={1:pool_data_cell_size(idcnt)};
                %                     try
                data_sign_array{pidx}(idx_vector{1},...
                    idx_vector{2},...
                    idx_vector{3},...
                    idx_vector{4},...
                    idx_vector{5},...
                    idx_vector{6},...
                    idx_vector{7},...
                    idx_vector{8},...
                    idx_vector{9},...
                    idx_vector{10},...
                    idx_vector{11},...
                    idx_vector{12})=...
                    merged_cells{idcnt};
                %                     catch
                %                         keyboard
                %                     end
                %                     %                 end
                %                 toc %Elapsed time is 0.009299 seconds.
                
                
                %
                %                 tic
                %                 substr=num2cell(idxcombis(idcnt,:)');
                %
                %                 substr=cellfun(@(x) num2str(x),substr,'UniformOutput',false);
                %                 substr{data_idx{pidx,2}}=['1:' num2str(pool_data_cell_size(idcnt))];
                %                 substr=strcat(substr,',')';
                %                 substr=[substr{:}];
                %                 substr=substr(1:end-1);
                %                 eval(['data_sign_array{pidx}(' substr ')=merged_cells{idcnt};']);
                %              toc
                
                
            elseif  last_statop_idx(pidx)>1 % if there are more operations to come, treat
                % data_sign_array the same as data_array.
                S.subs=num2cell(idxcombis(idcnt,:)');
                %                 S.subs(data_idx{pidx,2})={1:pool_data_cell_size(idcnt)};
                S.type='()';
                switch data_idx{pidx,1}
                    
                    case {'Mean','mean','MEAN'}
                        
                        data_sign_array{pidx} = subsasgn(data_sign_array{pidx},S,nanmean(merged_cells{idcnt}));
                    case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                        
                        data_sign_array{pidx} = subsasgn(data_sign_array{pidx},S,nanmean(abs(merged_cells{idcnt})));
                        
                    case {'Std','std','STD'}
                        
                        data_sign_array{pidx} = subsasgn(data_sign_array{pidx},S,nanstd(merged_cells{idcnt}));
                        
                    case {'Median','median','MEDIAN'}
                        data_sign_array{pidx} = subsasgn(data_sign_array{pidx},S,nanmedian(merged_cells{idcnt}));
                    case {'Hist','hist','HIST'}
                        data_sign_array{pidx} = subsasgn(data_sign_array{pidx},S,sum(~isnan(merged_cells{idcnt})));
                        
                end
            end
            
        end
        %         ...and for medians:
        for idcnt=1:size(idxcombis_Median,1)
            % The data_sign_array misses
            % the last statistical
            % operation:
            % If the first operation is
            % also the last AND the core
            % data has not simply been
            % pooled together (which could mean that pool_rest_cell still cotains cells
            % of various length) or if what I said in () is not true, store the
            % values from which the
            % mean/median/etc is calculated
            % above.
            % In the histogram case, we have allready calculated the medians
            if last_statop_idx(pidx)==1
                
                
                %                 for nd=1:pool_data_cell_size_Median(idcnt)
                %                     idx_vector=ones(1,12);
                %                     idx_vector(1,1:numel(idxcombis_Median(idcnt,:)))=idxcombis_Median(idcnt,:);
                %                     idx_vector(data_idx{pidx,2})=nd;
                %                     data_signMedian_array{pidx}(idx_vector(1),...
                %                         idx_vector(2),...
                %                         idx_vector(3),...
                %                         idx_vector(4),...
                %                         idx_vector(5),...
                %                         idx_vector(6),...
                %                         idx_vector(7))=...
                %                         merged_Median_cells{idcnt}(nd);
                %                 end
                %                 for nd=1:pool_data_cell_size_Median(idcnt)
                idx_vector=num2cell(ones(1,12));
                idx_vector(1,1:numel(idxcombis_Median(idcnt,:)))=num2cell(idxcombis_Median(idcnt,:));
                idx_vector(data_idx{pidx,2})={1:pool_data_cell_size_Median(idcnt)};
                data_signMedian_array{pidx}(idx_vector{1},...
                    idx_vector{2},...
                    idx_vector{3},...
                    idx_vector{4},...
                    idx_vector{5},...
                    idx_vector{6},...
                    idx_vector{7})=...
                    merged_Median_cells{idcnt};
                %                 end
                %                 tic
                %                 S.subs=num2cell(idxcombis_Median(idcnt,:)');
                %                 S.type='()';
                %                 S.subs(data_idx{pidx,2})={1:pool_data_cell_size_Median(idcnt)};
                %
                %                 data_signMedian_array{pidx} = ...
                %                      subsasgn(data_signMedian_array{pidx},S,merged_Median_cells{idcnt});
                %                 toc
            elseif  last_statop_idx(pidx)>1 % if there are more operations to come, treat
                % data_sign_array the same as data_array.
                S.subs=num2cell(idxcombis_Median(idcnt,:)');
                S.type='()';
                switch data_idx{pidx,1}
                    
                    case {'Mean','mean','MEAN'}
                        data_signMedian_array{pidx} = subsasgn(data_signMedian_array{pidx},S,nanmean(merged_Median_cells{idcnt}));
                    case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                        data_signMedian_array{pidx} = subsasgn(data_signMedian_array{pidx},S,nanmean(abs(merged_Median_cells{idcnt})));
                        
                    case {'Std','std','STD'}
                        data_signMedian_array{pidx} = subsasgn(data_signMedian_array{pidx},S,nanstd(merged_Median_cells{idcnt}));
                        
                    case {'Median','median','MEDIAN'}
                        data_signMedian_array{pidx} = subsasgn(data_signMedian_array{pidx},S,nanmedian(merged_Median_cells{idcnt}));
                    case {'Hist','hist','HIST'}
                        data_signMedian_array{pidx} = subsasgn(data_signMedian_array{pidx},S,nanmedian(merged_Median_cells{idcnt}));
                        
                end
            end
            
        end
        clear idxcombis_Median
        
        if last_statop_idx(pidx)==1
            switch lower(data_idx{pidx,1})
                case 'mean'
                    data_array{pidx} = nanmean(data_sign_array{pidx},data_idx{pidx,2});
                    data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,2});
                    data_Median_array{pidx} = nanmean(data_signMedian_array{pidx},data_idx{pidx,2});
                case 'meanabs'
                    data_array{pidx} = nanmean(abs(data_sign_array{pidx},data_idx{pidx,2}));
                    data_dev_array{pidx} = nanstd(abs(data_sign_array{pidx},[],data_idx{pidx,2}));
                    data_Median_array{pidx} = nanmean(abs(data_signMedian_array{pidx},data_idx{pidx,2}));
                case 'std'
                    data_array{pidx} = nanstd(data_sign_array{pidx},data_idx{pidx,2});
                    data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,2});
                    data_Median_array{pidx} = nanstd(data_signMedian_array{pidx},data_idx{pidx,2});
                case 'var'
                    data_array{pidx} = nanvar(data_sign_array{pidx},data_idx{pidx,2});
                    data_dev_array{pidx} = nanvar(data_sign_array{pidx},[],data_idx{pidx,2});
                    data_Median_array{pidx} = nanvar(data_signMedian_array{pidx},data_idx{pidx,2});
                    
                case 'median'
                    data_dev_array{pidx} = prctile(data_sign_array{pidx}, [prctl 100-prctl],data_idx{pidx,2});
                    data_array{pidx} =  nanmedian(data_sign_array{pidx},data_idx{pidx,2});
                    data_Median_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,2});
                case 'sum'
                    data_dev_array{pidx} = idSocial_auxiliaries_nansum(data_sign_array{pidx},data_idx{pidx,2});
                    data_array{pidx} =  idSocial_auxiliaries_nansum(data_sign_array{pidx},data_idx{pidx,2});
                    data_Median_array{pidx} = idSocial_auxiliaries_nansum(data_signMedian_array{pidx},data_idx{pidx,2});
                case 'angle_ratio'
                    data_dev_array{pidx} = prctile(data_sign_array{pidx}, [prctl 100-prctl],data_idx{pidx,2});
                    data_array{pidx} =  idSocial_angle_ratio(data_sign_array{pidx},data_idx{pidx,2},pi/2);
                    data_Median_array{pidx} = idSocial_angle_ratio(data_signMedian_array{pidx},data_idx{pidx,2},pi/2);
%                 case 'rank'
%                    
%                     sortvals = data_sign_array{pidx};
%                     sortvals(isnan(sortvals)) = -inf;
%                     sortvals = sort(sortvals,1,'descend');
%                     sortvals(isinf(sortvals)) = NaN;
%                     
%                     sortvals_format4all = NaN(1,cell4statistics_max_elems);
%                     sortvals_format4all(1:numel(sortvals)) = sortvals;
%                     data_array{pidx} =  sortvals_format4all;
%                     data_Median_array{pidx} = nanmean(sortvals_format4all(1:5));
                    
                case 'positive_ratio'
                    data_dev_array{pidx} = prctile(data_sign_array{pidx}, [prctl 100-prctl],data_idx{pidx,2});
                    data_array{pidx} =  idSocial_positive_ratio(data_sign_array{pidx},data_idx{pidx,2},0);
                    data_Median_array{pidx} = idSocial_positive_ratio(data_signMedian_array{pidx},data_idx{pidx,2},0);
                case 'mode'
                    if isfield(plot_mode,'mode_edges') && ~isempty(plot_mode.mode_edges) && ...
                            isfield(plot_mode,'mode_bandwith') && ~isempty(plot_mode.mode_bandwith)
                        
                        szref=size(data_array{pidx});
                        szref(data_idx{pidx,2})=1;
                        peak_cmbs=vectorsize2indexcombinations(szref);
                        
                        S.type='()';
                        for act=1:size(peak_cmbs,1)
                            S.subs=num2cell(peak_cmbs(act,:)');
                            Sref.subs=S.subs;
                            Sref.type='()';
                            Sref.subs(data_idx{pidx,2})={':'};
                            
                            temp=idSocial_distribution_mode(subsref(data_sign_array{pidx},Sref),plot_mode.mode_edges,plot_mode.mode_bandwith);
                            Sref.subs(data_idx{pidx,2})={1};
                            data_array{pidx}=subsasgn(data_array{pidx},Sref,temp);
                            
                            
                        end
                        szrefMedian= size(data_signMedian_array{pidx});
                        szrefMedian(data_idx{pidx,2})=1;
                        peak_cmbs_Median=vectorsize2indexcombinations(szrefMedian);
                        S.type='()';
                        for act=1:size(peak_cmbs_Median,1)
                            S.subs=num2cell(peak_cmbs_Median(act,:)');
                            Sref.subs=S.subs;
                            Sref.type='()';
                            Sref.subs(data_idx{pidx,2})={':'};
                            
                            temp=idSocial_distribution_mode(subsref(data_signMedian_array{pidx},Sref),plot_mode.mode_edges,plot_mode.mode_bandwith);
                            data_Median_array{pidx}=subsasgn(data_Median_array{pidx},Sref,temp);
                        end
                        
                        
                    else
                        error([mfilename ': plot_mode.mode_edges and plot_mode.mode_bandwith for calculation of distribution and mode are undefined'])
                    end
                    data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,2});
                case '1stmode'
                    if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges) && ...
                            isfield(plot_mode,'dist_bandwidth') && ~isempty(plot_mode.dist_bandwidth)
                        
                        szref=size(data_array{pidx});
                        szref(data_idx{pidx,2})=1;
                        peak_cmbs=vectorsize2indexcombinations(szref);
                        
                        S.type='()';
                        for act=1:size(peak_cmbs,1)
                            S.subs=num2cell(peak_cmbs(act,:)');
                            Sref.subs=S.subs;
                            Sref.type='()';
                            Sref.subs(data_idx{pidx,2})={':'};
                            
                            temp=idSocial_distribution_mode(subsref(data_sign_array{pidx},Sref),plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method,'first');
                            Sref.subs(data_idx{pidx,2})={1};
                            data_array{pidx}=subsasgn(data_array{pidx},Sref,temp);
                            
                            
                        end
                        szrefMedian= size(data_signMedian_array{pidx});
                        szrefMedian(data_idx{pidx,2})=1;
                        peak_cmbs_Median=vectorsize2indexcombinations(szrefMedian);
                        S.type='()';
                        for act=1:size(peak_cmbs_Median,1)
                            S.subs=num2cell(peak_cmbs_Median(act,:)');
                            Sref.subs=S.subs;
                            Sref.type='()';
                            Sref.subs(data_idx{pidx,2})={':'};
                            
                            temp=idSocial_distribution_mode(subsref(data_signMedian_array{pidx},Sref),plot_mode.edges, ...
                                                                plot_mode.dist_bandwidth,plot_mode.dist_method,'first');
                            data_Median_array{pidx}=subsasgn(data_Median_array{pidx},Sref,temp);
                        end
                        
                        
                    else
                        error([mfilename ': plot_mode.edges and plot_mode.dist_bandwith for calculation of distribution and mode are undefined'])
                    end
                    data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,2});
                        
                case 'mode_height'
                    if isfield(plot_mode,'mode_edges') && ~isempty(plot_mode.mode_edges) && ...
                            isfield(plot_mode,'mode_bandwith') && ~isempty(plot_mode.mode_bandwith)
                        
                        szref=size(data_array{pidx});
                        szref(data_idx{pidx,2})=1;
                        peak_cmbs=vectorsize2indexcombinations(szref);
                        
                        S.type='()';
                        for act=1:size(peak_cmbs,1)
                            S.subs=num2cell(peak_cmbs(act,:)');
                            Sref.subs=S.subs;
                            Sref.type='()';
                            Sref.subs(data_idx{pidx,2})={':'};
                            
                            temp=idSocial_distribution_modeHeight(subsref(data_sign_array{pidx},Sref),plot_mode.mode_edges,plot_mode.mode_bandwith);
                            Sref.subs(data_idx{pidx,2})={1};
                            data_array{pidx}=subsasgn(data_array{pidx},Sref,temp);
                            
                            
                        end
                        szrefMedian= size(data_signMedian_array{pidx});
                        szrefMedian(data_idx{pidx,2})=1;
                        peak_cmbs_Median=vectorsize2indexcombinations(szrefMedian);
                        S.type='()';
                        for act=1:size(peak_cmbs_Median,1)
                            S.subs=num2cell(peak_cmbs_Median(act,:)');
                            Sref.subs=S.subs;
                            Sref.type='()';
                            Sref.subs(data_idx{pidx,2})={':'};
                            
                            temp=idSocial_distribution_modeHeight(subsref(data_signMedian_array{pidx},Sref),plot_mode.mode_edges,plot_mode.mode_bandwith);
                            data_Median_array{pidx}=subsasgn(data_Median_array{pidx},Sref,temp);
                        end
                        
                        
                    else
                        error([mfilename ': plot_mode.mode_edges and plot_mode.mode_bandwith for calculation of distribution and mode are undefined'])
                    end
                    data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,2});
                    
                case 'hist'
                    data_array{pidx}= idSocial_auxiliaries_nansum(data_sign_array{pidx},data_idx{pidx,2});
                    data_Median_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,2});
                    data_dev_array{pidx} = prctile(data_sign_array{pidx}, [prctl 100-prctl],data_idx{pidx,2});
                case 'hist+mean'
                    data_array{pidx}= idSocial_auxiliaries_nansum(data_sign_array{pidx},data_idx{pidx,2});
                    data_Median_array{pidx} = nanmean(data_signMedian_array{pidx},data_idx{pidx,2});
                    data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,2});
                case 'hist+mode'
                    t2=tic;
                    data_array{pidx}= idSocial_auxiliaries_nansum(data_sign_array{pidx},data_idx{pidx,2});
                    data_Median_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,2});
                    data_dev_array{pidx} = prctile(data_sign_array{pidx}, [prctl 100-prctl],data_idx{pidx,2});
                    
                    
                    % The following is quite 'cutre':
                    
                    try
                        if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                            szref=size(data_array{pidx});
                            szref(7)=1;
                            peak_cmbs=vectorsize2indexcombinations(szref);
                            locs=NaN(szref);
                            S.type='()';
                            for act=1:size(peak_cmbs,1)
                                S.subs=num2cell(peak_cmbs(act,:)');
                                Sref.subs=S.subs;
                                Sref.type='()';
                                Sref.subs(7)={':'};
                                hitemp=squeeze(subsref(data_array{pidx},Sref));
                                
                                
%                                 [~,locs_temp] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
                                [~,locs_temp] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);

                                if ~isempty(locs_temp)
                                    locs=subsasgn(locs,S, plot_mode.edges(locs_temp));
                                else
                                    locs = NaN;
                                end
                                
                            end
                            data_Median_array{pidx} = locs;
                        end
                    catch
                        keyboard
                    end
                    toc(t2)
                case 'hist+1stmode'
                    
                    
                    data_array{pidx}= idSocial_auxiliaries_nansum(data_sign_array{pidx},data_idx{pidx,2});
                    data_Median_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,2});
                    data_dev_array{pidx} = prctile(data_sign_array{pidx}, [prctl 100-prctl],data_idx{pidx,2});
                    
                    
                    % The following is quite 'cutre':
                    
                    try
                        if isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
                            szref=size(data_array{pidx});
                            szref(7)=1;
                            peak_cmbs=vectorsize2indexcombinations(szref);
                            locs=NaN(szref);
                            S.type='()';
                            for act=1:size(peak_cmbs,1)
                                S.subs=num2cell(peak_cmbs(act,:)');
                                Sref.subs=S.subs;
                                Sref.type='()';
                                Sref.subs(7)={':'};
                                hitemp=squeeze(subsref(data_array{pidx},Sref));
                                
                                
%                                 [~,locs_temp] = findpeaks(hitemp);
                                [~,locs_temp] = idSocial_auxiliaries_findpeaks(hitemp);
                                if ~isempty(locs_temp)
                                    locs=subsasgn(locs,S, plot_mode.edges(min(locs_temp)));
                                else
                                    locs = NaN;
                                end
                                
                            end
                            data_Median_array{pidx} = locs;
                        end
                    catch
                        keyboard
                    end
                case 'rankdist+mean'
                    warning([mfilename ': ''rankdist+mean'' with fixed number of points = 5. There might be more issues!'])
                    pooldata = data_sign_array{pidx};
                    ndim = ndims(pooldata);
                    permvec = [setxor(1:ndim,data_idx{pidx,2}) data_idx{pidx,2}];
                    pooldata = permute(pooldata,permvec);
                    pooldata = pooldata(:,:,:,:,:,:);
                    pooldata(isnan(pooldata)) = -inf;
                    pooldata = sort(pooldata,6,'descend');
                    pooldata(isinf(pooldata)) = NaN;
                    pooldata = pooldata(:,:,:,:,:,1:size(data_sign_array{pidx},7));
                    [~,idxs]=sort(permvec);
                    pooldata = permute(pooldata,idxs);
                    data_array{pidx}= pooldata;
                    data_Median_array{pidx} = nanmean( pooldata(:,:,:,:,:,1:5),6);
                    data_dev_array{pidx} = nanstd( pooldata(:,:,:,:,:,1:5),[],6);
                    
            end
            statistics_on_idx(pidx) = data_idx{pidx,2};
            statistics_type{pidx} = data_idx{pidx,1};
          
            % Normalization for
            % histograms/distributions
            if isfield(plot_mode,'normalization') && ~isempty(plot_mode.statistics{pidx}) &&  ...
                    strcmpi(plot_mode.statistics{pidx},'pool') && ...
                    ((~isempty(strfind(lower(data_idx{pidx,1}),'hist')) && strcmpi(plot_mode.display_mode,'hist')) || ...
                    strcmpi(plot_mode.display_mode,'MapPolar') ||...
                    strcmpi(plot_mode.display_mode,'Map') ...
                    ) 
                
                clear no_frames no_frames_effective edgArea
                if strcmpi(plot_mode.display_mode,'hist') %&& strcmpi(plot_mode.normalization{pidx},'density')
                    edgArea = plot_mode.edges(2)-plot_mode.edges(1);
                    no_frames = sum(squeeze(data_array{pidx}(:)));
                    no_frames_effective = sum(squeeze(data_array{pidx}));
                    
                    
                elseif strcmpi(plot_mode.display_mode,'Map') %&& strcmpi(plot_mode.normalization{pidx},'density')
                    no_frames = sum(squeeze(data_array{pidx}(:)));
                    no_frames_effective = sum(squeeze(data_array{pidx}(:)));
                    edgArea = (plot_mode.edges{1}(2)-plot_mode.edges{1}(1))*(plot_mode.edges{2}(2)-plot_mode.edges{2}(1));
                elseif strcmpi(plot_mode.display_mode,'MapPolar')% && strcmpi(plot_mode.normalization{pidx},'density')
                    no_frames = sum(squeeze(data_array{pidx}(:)));
                    no_frames_effective = sum(squeeze(data_array{pidx}(:)));
                    ar = pi*plot_mode.edges{1}(2:end).^2 - pi*plot_mode.edges{1}(1:end-1).^2;
                    ft =     (plot_mode.edges{2}(2)-plot_mode.edges{2}(1))/(2*pi);
                    edgArea  = repmat(ar .* ft,[size(plot_mode.edges{2},2)-1,1]);
                    
                end
                
                
                switch lower(plot_mode.normalization{pidx})
                    case 'none'
                    case 'no_frames'
                        data_array{pidx} = data_array{pidx}/no_frames;
                    case 'no_frames_effective'
                        data_array{pidx} = data_array{pidx}/no_frames_effective;
                    case 'density'
                       
                        if ~strcmpi(plot_mode.display_mode,'MapPolar')
                            data_array{pidx} = data_array{pidx}/no_frames/no_frames./edgArea;
                        else
                            data_array{pidx}(1,1,1,1,1,1,:,:) = squeeze(data_array{pidx})./edgArea./no_frames;
                        end
                        
                end
            end
            
        elseif last_statop_idx(pidx)>1                      % if there are more operations to come,
            data_array{pidx}=data_sign_array{pidx};      % data_sign_array == data_array at this moment.
            data_Median_array{pidx} =data_signMedian_array{pidx};
        end
        
        
        switch data_idx{pidx,1}
            case {'Mean','mean','MEAN'}
                dim_names_act{pidx,data_idx{pidx,2}(1)}=['MEAN(' dim_names_act{pidx,data_idx{pidx,2}(1)} ')'];
            case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                dim_names_act{pidx,data_idx{pidx,2}(1)}=['MEAN(ABS(' dim_names_act{pidx,data_idx{pidx,2}(1)} '))'];
                
            case {'Std','std','STD'}
                dim_names_act{pidx,data_idx{pidx,2}(1)}=['STD(' dim_names_act{pidx,data_idx{pidx,2}(1)} ')'];
            case {'Median','median','MEDIAN'}
                dim_names_act{pidx,data_idx{pidx,2}(1)}=['MEDIAN(' dim_names_act{pidx,data_idx{pidx,2}(1)} ')'];
        end
        
    end
    
end
% ...and the rest of the operations on this array:
data_string=cell(no_filterentries,1);

for pidx=1:no_filterentries
    if ~isempty(data_idx{pidx}) && ~all(cellfun(@(x) isempty(x),data_idx(pidx,3:end))) && size(data_idx(pidx,:),2)>3
        for op_ind=3:2:size(data_idx(pidx,:),2)
            if ~isempty(data_idx{pidx,op_ind})
                
                switch data_idx{pidx,op_ind}
                    case {'Mean','mean','MEAN'}
                        
                        
                        dim_names_act{pidx,data_idx{pidx,2}(1)}=[ 'MEAN(' dim_names_act{pidx,data_idx{pidx,op_ind+1}}  '[' filter_complete{pidx,2*data_idx{pidx,op_ind+1}} '],' dim_names_act{pidx,data_idx{pidx,2}(1)}  ')'];
                        dim_names_act{pidx,data_idx{pidx,op_ind+1}}='';
                        
                        data_array{pidx} = nanmean( data_array{pidx},data_idx{pidx,op_ind+1});
                        data_Median_array{pidx} = nanmean( data_Median_array{pidx},data_idx{pidx,op_ind+1});
                     case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                        
                        
                        dim_names_act{pidx,data_idx{pidx,2}(1)}=[ 'MEAN(ABS(' dim_names_act{pidx,data_idx{pidx,op_ind+1}}  '[' filter_complete{pidx,2*data_idx{pidx,op_ind+1}} '],' dim_names_act{pidx,data_idx{pidx,2}(1)}  '))'];
                        dim_names_act{pidx,data_idx{pidx,op_ind+1}}='';
                        
                        data_array{pidx} = nanmean( abs(data_array{pidx},data_idx{pidx,op_ind+1}));
                        data_Median_array{pidx} = nanmean(abs( data_Median_array{pidx},data_idx{pidx,op_ind+1}));
                           
                    case {'Std','std','STD'}
                        
                        
                        dim_names_act{pidx,data_idx{pidx,2}(1)}=[ 'STD(' dim_names_act{pidx,data_idx{pidx,op_ind+1}}  '[' filter_complete{pidx,2*data_idx{pidx,op_ind+1}} '],' dim_names_act{pidx,data_idx{pidx,2}(1)}  ')'];
                        dim_names_act{pidx,data_idx{pidx,op_ind+1}}='';
                        
                        data_array{pidx} = nanstd( data_array{pidx},data_idx{pidx,op_ind+1});
                        data_Median_array{pidx} = nanstd( data_Median_array{pidx},data_idx{pidx,op_ind+1});
                        
                        
                    case {'Median','median','MEDIAN'}
                        dim_names_act{pidx,data_idx{pidx,2}(1)}=[ 'MEDIAN(' dim_names_act{pidx,data_idx{pidx,op_ind+1}} '[' filter_complete{pidx,2*data_idx{pidx,op_ind+1}} '],' dim_names_act{pidx,data_idx{pidx,2}(1)}  ')'];
                        dim_names_act{pidx,data_idx{pidx,op_ind+1}}='';
                        
                        data_array{pidx} = nanmedian( data_array{pidx},data_idx{pidx,op_ind+1});
                        data_Median_array{pidx} = nanmedian(data_Median_array{pidx},data_idx{pidx,op_ind+1});
                        
                    case {'Hist','hist','HIST'}
                        dim_names_act{pidx,data_idx{pidx,2}(1)}=[ 'HIST(' dim_names_act{pidx,data_idx{pidx,op_ind+1}} '[' filter_complete{pidx,2*data_idx{pidx,op_ind+1}} '],' dim_names_act{pidx,data_idx{pidx,2}(1)}  ')'];
                        dim_names_act{pidx,data_idx{pidx,op_ind+1}}='';
                        
                        data_array{pidx}= sum(~isnan(data_sign_array{pidx}),data_idx{pidx,op_ind+1});
                        data_Median_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,op_ind+1});
                    case {'Hist+1stMode','hist+1stmode','HIST+1STMODE'}
                        dim_names_act{pidx,data_idx{pidx,2}(1)}=[ 'HIST+1STMODE(' dim_names_act{pidx,data_idx{pidx,op_ind+1}} '[' filter_complete{pidx,2*data_idx{pidx,op_ind+1}} '],' dim_names_act{pidx,data_idx{pidx,2}(1)}  ')'];
                        dim_names_act{pidx,data_idx{pidx,op_ind+1}}='';
                        
                        data_array{pidx}= sum(~isnan(data_sign_array{pidx}),data_idx{pidx,op_ind+1});
                        data_Median_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,op_ind+1});
                        
                end
                % Leave out last statistical operation on the
                % data_sign_array/execute statistical op. only if it is not the last one to come:
                
                if op_ind+1~=(last_statop_idx(pidx)-1)*2+2
                    switch data_idx{pidx,op_ind}
                        case {'Mean','mean','MEAN'}
                            data_sign_array{pidx} = nanmean(data_sign_array{pidx},data_idx{pidx,op_ind+1});
                            data_signMedian_array{pidx} = nanmean( data_signMedian_array{pidx},data_idx{pidx,op_ind+1});
                            
                        case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                            data_sign_array{pidx} = nanmean(abs(data_sign_array{pidx},data_idx{pidx,op_ind+1}));
                            data_signMedian_array{pidx} = nanmean( abs(data_signMedian_array{pidx},data_idx{pidx,op_ind+1}));
                            
                        case {'Std','std','STD'}
                            data_sign_array{pidx} = nanstd(data_sign_array{pidx},data_idx{pidx,op_ind+1});
                            data_signMedian_array{pidx} = nanstd( data_signMedian_array{pidx},data_idx{pidx,op_ind+1});
                            
                        case {'Median','median','MEDIAN'}
                            data_sign_array{pidx} = nanmedian( data_sign_array{pidx},data_idx{pidx,op_ind+1});
                            data_signMedian_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,op_ind+1});
                        case {'Hist','hist','HIST'}
                            data_sign_array{pidx}= sum(~isnan(data_sign_array{pidx}),data_idx{pidx,op_ind+1});
                            data_signMedian_array{pidx} = nanmedian(data_signMedian_array{pidx},data_idx{pidx,op_ind+1});
                            
                            
                    end
                else % if the last operation is reached, calulate confidence intervals/standard error
                    switch data_idx{pidx,op_ind}
                        case {'Mean','mean','MEAN'}
                            
                            data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,op_ind+1});
                        case {'MeanAbs','meanabs','MEANABS','Meanabs','meanAbs'}
                            
                            data_dev_array{pidx} = nanstd(abs(data_sign_array{pidx},[],data_idx{pidx,op_ind+1}));
                            
                        case {'Std','std','STD'}
                            % This is actually not the
                            % deviation or anything...what
                            % is the error of the standard
                            % deviation? Just put this here
                            % to have something
                            data_dev_array{pidx} = nanstd(data_sign_array{pidx},[],data_idx{pidx,op_ind+1});
                            
                        case {'Median','median','MEDIAN'}
                            
                            data_dev_array{pidx} = prctile(data_sign_array{pidx}, [prctl 100-prctl],data_idx{pidx,op_ind+1});
                        case  {'Hist','hist','HIST'}
                            
                            data_dev_array{pidx} = prctile(data_signMedian_array{pidx}, [prctl 100-prctl],data_idx{pidx,op_ind+1});
                            
                    end
                    statistics_type{pidx} = data_idx{pidx,op_ind};
                    statistics_on_idx(pidx) = data_idx{pidx,op_ind+1};
                end
            end
        end
    end
    
    data_string{pidx}=dim_names_act{pidx,data_idx{pidx,2}(1)};
    
end
clear pool_rest_cell
% Sizes of data_sign_array and data_array should only diverge in one
% dimension (the dimension of the last statistical operation applied to the data);
% only exception: the dimension of the last statistical operation is one.
for lg=1:no_filterentries
    if ~isempty(data_sign_array{lg})
        ndm=ndims(data_array{lg});
        ndm_sign=ndims(data_sign_array{lg});
        sz=[size(data_array{lg}) ones(1,ndm_sign-ndm)];
        sz_sign=size(data_sign_array{lg});
        if  sum(sz~=sz_sign)~=1 && ...
                size(data_array{lg},data_idx{lg,2*(last_statop_idx(lg)-1)+2})>1
            error([mfilename ': SizeDifference'],[mfilename ': Sizes of data_sign_array and data_array should diverge in one dimension.']);
        end
    end
end
%% 'Chapter II': Selecting data for figures, subplots and filter
% filter data depends on subplot data depends on figure data, or: First
% select, for which dimensions of data array to open a new figure, then for
% which of the dimensions of the thus restricted data create a new subplot,
% and again for filters.
% Note that each subplot and figure will have the same filter items, i.e.,
% while for example a new subplot is created for each trial of each group,
% all the subplots show data for focals (it is not possible to show data
% for focals in one and data for neighbors in the other subplot.

%% Select filter indices
S.type='()';
autocomb_idx_final=[];
filter_idx_final=[];


% data_array_final with dimensions [total_no_autocombentries,1]
data_array_final={};%cell(1,1);
data_Median_array_final={};%cell(1,1);
data_sign_array_final={};
data_signMedian_array_final={};
data_dev_array_final={};
dim_names_final={};
data_string_final={};
autocomb_string_final={};
statistics_on_idx_final=[];
statistics_type_final={};
lg_count=1;
for lg=1:no_filterentries
    if ~isempty(autocomb_idx{lg,1})
        lg_combi2=[];
        no_lgs=1;
        if ~isempty(autocomb_idx{lg,1})
            no_lgs_sep=NaN(size(autocomb_idx{lg,1},2));
        else
            no_lgs_sep=1;
        end
        for k=1:size(autocomb_idx{lg,1},2)
            no_lgs=no_lgs*size(data_array{lg},autocomb_idx{lg,1}(k));
            no_lgs_sep(k)=size(data_array{lg},autocomb_idx{lg,1}(k));
        end
        if ~isempty(no_lgs_sep)
            sets=cell(size(no_lgs_sep,1),1);
            for k=1:size(no_lgs_sep,2)
                sets{k}=1:no_lgs_sep(k);
            end
            c = cell(1, numel(sets));
            [c{:}] = ndgrid( sets{end:-1:1} );
            lg_combi = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
            lg_combi = lg_combi(:,end:-1:1);
            for lt=1:size(lg_combi,1)
                lg_combi2{lt}=lg_combi(lt,:);
            end
            lg_combi=lg_combi2';
            no_combis=size(lg_combi,1);
            
            autocomb_idx_final=vertcat(autocomb_idx_final,...
                [repmat(autocomb_idx(lg),no_combis,1) lg_combi]);
            if  ~isempty(filter_idx) && ~isempty(filter_idx(lg,:))
                filter_idx_final=vertcat(filter_idx_final,...
                    [repmat(filter_idx(lg,:),no_combis,1)]);
            end
            if ~isempty(plot_mode.autocombs) && ~isempty(plot_mode.autocombs(lg,1))
                autocomb_string_final=vertcat(autocomb_string_final,...
                    [repmat(plot_mode.autocombs(lg,1),no_combis,1) lg_combi]);
            end
        else
            lg_combi=1;
        end
        %         %%%%%%%%%%%%%%
        %
        for cmb=1:no_combis
            S.subs=repmat({':'},max(ndims(data_sign_array{lg}),ndims(data_array{lg})),1);
            for k=1:size(autocomb_idx_final{lg_count,1},2)
                S.subs(autocomb_idx_final{lg_count,1}(k))={autocomb_idx_final{lg_count,2}(k)};
            end
            try
                data_array_final=vertcat(data_array_final,...
                    subsref(data_array{lg},S));
                data_Median_array_final=vertcat(data_Median_array_final,...
                    subsref(data_Median_array{lg},S));
                if ~isempty(data_sign_array{lg})
                    data_sign_array_final=vertcat(data_sign_array_final,...
                        subsref(data_sign_array{lg},S));
                end
                if ~isempty(data_signMedian_array{lg})
                    data_signMedian_array_final=vertcat(data_signMedian_array_final,...
                        subsref(data_signMedian_array{lg},S));
                end
                if ~isempty(data_dev_array{lg})
                    data_dev_array_final=vertcat(data_dev_array_final,...
                        subsref(data_dev_array{lg},S));
                end
                statistics_on_idx_final=vertcat(statistics_on_idx_final,...
                    statistics_on_idx(lg));
                statistics_type_final=vertcat(statistics_type_final,...
                    statistics_type{lg});
            catch
                keyboard
            end
            lg_count=lg_count+1;
        end
        dim_names_final=vertcat(dim_names_final,...
            repmat(dim_names_act(lg,:),no_combis,1));
        data_string_final=vertcat(data_string_final,...
            repmat(data_string(lg),no_combis,1));
    else
        data_array_final=data_array;%cell(1,1);
        data_Median_array_final=data_Median_array;%cell(1,1);
        data_sign_array_final=data_sign_array;
        data_signMedian_array_final=data_signMedian_array;
        data_dev_array_final=data_dev_array;
        dim_names_final=dim_names_act;
        data_string_final=data_string;
        autocomb_string_final=[];
        statistics_on_idx_final=statistics_on_idx;
        statistics_type_final=statistics_type;
    end
end

% If autocombs contains dimensions which also appear in filter, update idices in
% filter_idx_final:
no_rows=size(autocomb_idx_final,1);

for row=1:no_rows
    if ~isempty(filter_idx_final) && ~isempty(autocomb_idx_final) && ...
            ~isempty(filter_idx_final(row,:)) && ~isempty(autocomb_idx_final{row,1})
        fdim=filter_idx_final(row,2:2:end);
        fid=[filter_idx_final{row,1:2:end}];
        
        combid=autocomb_idx_final{row,1};
        combdim=autocomb_idx_final{row,2};
        
        [~,ia,ib]=intersect(combid,fid);
        
        filter_idx_final(row,ib*2-1 + 1)=num2cell(combdim(ia));
    end
end



%% Select subplot indices
S.type='()';
subplot_idx_final=[];

% Generate possible combinations of indices for subplots
if size(subplot_idx,2)>1 && ~isempty(subplot_idx{1,2}) % Indices are given explicitly
    subplot_idx_final=subplot_idx;
elseif isempty(subplot_idx) || isempty(subplot_idx{1})
    subplot_idx_final=[];
else
    for sp=1:size(subplot_idx,1)
        sp_combi2=[];
        no_sps=1;
        if ~isempty(subplot_idx{sp,1})
            no_sps_sep=NaN(size(subplot_idx{sp,1},2));
        else
            no_sps_sep=1;
        end
        for k=1:size(subplot_idx{sp,1},2)
            no_sps_sep(k)=size(data_array{sp},subplot_idx{sp,1}(k));
        end
        
        if ~isempty(no_sps_sep)
            sets=cell(size(no_sps_sep,1),1);
            for k=1:size(no_sps_sep,2)
                sets{k}=1:no_sps_sep(k);
            end
            c = cell(1, numel(sets));
            [c{:}] = ndgrid( sets{end:-1:1} );
            sp_combi = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
            sp_combi = sp_combi(:,end:-1:1);
            for lt=1:size(sp_combi,1)
                sp_combi2{lt}=sp_combi(lt,:);
            end
            sp_combi=sp_combi2';
            no_combis=size(sp_combi,1);
            subplot_idx_final=vertcat(subplot_idx_final,...
                [repmat(subplot_idx(sp),no_combis,1) sp_combi]);
            
        else
            sp_combi=1;
        end
    end
end

%% Data description string

no_filterentries_final=size(data_array_final,1);
% no_filterentries_final=no_filterentries;
filter_string_final=cell(size(filter_idx_final,1),1);


for lg=1:no_filterentries_final
    autocomb_string2=[];
    if ~isempty(autocomb_idx_final) && ~isempty(autocomb_idx_final{lg})
        
        for k=1:size(autocomb_idx_final{lg,1},2)
            autocomb_string2=[autocomb_string2 dim_names_short{autocomb_idx_final{lg,1}(k)} ' ' num2str(autocomb_idx_final{lg,2}(k)) '; '];
        end
    end
    
    
    if ~isempty(autocomb_string2) && ~isempty(filter_string_final) &&...
            (size(filter_string_final,1)>=lg && ~isempty(filter_string_final{lg}))
        filter_string_final{lg,1}=[autocomb_string2 filter_string_final{lg}(1:end-2)];
    elseif ~isempty(autocomb_string2) && (isempty(filter_string_final) || (size(filter_string_final,1)<lg || isempty(filter_string_final{lg})))
        filter_string_final{lg,1}=autocomb_string2;
    elseif ~isempty(filter_idx)
        
        filter_string_final{lg,1}='';
        for fe = 1:2:size(filter_idx{lg},2)
            filter_string_final{lg,1} = [filter_string_final{lg,1} dim_names_short{filter_idx{lg,fe}} ' ' num2str(filter_idx{lg,fe+1}) ' '];
        end
        filter_string_final{lg,1}=filter_string_final{lg,1}(1:end-1);
        filter_string_final{lg,1} = [filter_string_final{lg,1} ':'];
    else
        filter_string_final{lg,1}='';
    end
end
% clear cell_permute stats_cell

data_string_final=strcat(filter_string_final,{' '}, data_string_final);

% if ~isempty(filter_string_final{lg})
% else
%     data_string_final=strcat(filter_string_final,{' '},'-',{' '}, data_string_final);
% end

%%
% Generate (no_subplots,no_filterentries_final)-cells with all the data and
% info correpsonding to each subplot and filter.
no_subplots=size(subplot_idx_final,1);
spdata=cell(no_subplots,no_filterentries_final);
spdata_Median=cell(no_subplots,no_filterentries_final);
spdata_sign=cell(no_subplots,no_filterentries_final);
spdata_signMedian=cell(no_subplots,no_filterentries_final);
spdata_dev=cell(no_subplots,no_filterentries_final);
sp_filter_string=cell(no_subplots,no_filterentries_final);
sp_legend_string=cell(no_subplots,no_filterentries_final);
sp_data_string=cell(no_subplots,no_filterentries_final);
sp_statistics_on_idx=NaN(no_subplots,no_filterentries_final);
sp_statistics_type=cell(no_subplots,no_filterentries_final);
sp_dim_names=cell(no_subplots,no_filterentries_final);
sp_filter_idx=cell(no_subplots,1);
if  ~isempty(subplot_idx_final)
    for sp=1:no_subplots
        
        for lg=1:no_filterentries_final
            S.subs=repmat({':'},max(ndims(data_array_final{lg}),ndims(data_sign_array_final{lg})),1);
            for k=1:size(subplot_idx_final{sp,1},2)
                S.subs(subplot_idx_final{sp,1}(k))={subplot_idx_final{sp,2}(k)};
            end
            spdata{sp,lg}=subsref(data_array_final{lg},S);
            spdata_Median{sp,lg}=subsref(data_Median_array_final{lg},S);
            if ~isempty(data_sign_array_final)
                spdata_sign{sp,lg}=subsref(data_sign_array_final{lg},S);
            end
            if ~isempty(data_signMedian_array_final)
                spdata_signMedian{sp,lg}=subsref(data_signMedian_array_final{lg},S);
            end
            if ~isempty(data_dev_array_final)
                spdata_dev{sp,lg}=subsref(data_dev_array_final{lg},S);
            end
            sp_filter_string{sp,lg}=filter_string_final{lg,:};
            sp_filter_idx{sp}=filter_idx_final;
            
            subplot_string2=[];
            for k=1:size(subplot_idx_final{sp,1},2)
                subplot_string2=[subplot_string2 dim_names_short{subplot_idx_final{sp,1}} ' ' num2str(subplot_idx_final{sp,2}) '; '];
            end
            
            
            sp_data_string{sp,lg}=[legend_identifiers{lg}, ': ',  subplot_string2,data_string_final{lg,:}];
            sp_dim_names{sp,lg}=dim_names_final(lg,:);
            sp_statistics_on_idx(sp,lg)=statistics_on_idx_final(lg);
            sp_statistics_type{sp,lg}=statistics_type_final{lg};
        end
        if ~isfield(plot_mode,'legendstring')
            plot_mode.legendstring=[];
        end
        if isfield(plot_mode,'legendstring')  && ~isempty(plot_mode.legendstring) && size(plot_mode.legendstring,1)==1 && no_subplots>1
            plot_mode.legendstring=repmat(plot_mode.legendstring,[no_subplots 1 1 1]);
        end
        
        if isempty(plot_mode.legendstring) || size(plot_mode.legendstring,1) < sp || size(plot_mode.legendstring(sp,:),2)<no_filterentries_final
            if isfield(plot_mode,'legendstring') && numel(plot_mode.legendstring)<no_filterentries_final
                sp_legend_string(sp,:)=legend_identifiers(1:no_filterentries_final)';
            elseif numel(plot_mode.legendstring)>=no_filterentries_final && any(size(plot_mode.legendstring))==1
                sp_legend_string(sp,:)=plot_mode.legendstring(1:no_filterentries_final);
            end
        else
            sp_legend_string(sp,:)=plot_mode.legendstring(sp,:);
        end
    end
else
    sp=1;
    
    for lg=1:no_filterentries_final
        
        spdata{sp,lg}=data_array_final{lg};
        spdata_Median{sp,lg}=data_Median_array_final{lg};
        if ~isempty(data_sign_array_final)
            spdata_sign{sp,lg}=data_sign_array_final{lg};
        end
        if ~isempty(data_signMedian_array_final)
            spdata_signMedian{sp,lg}=data_signMedian_array_final{lg};
        end
        if ~isempty(data_dev_array_final)
            spdata_dev{sp,lg}=data_dev_array_final{lg};
        else
            
        end
        
        sp_filter_string{sp,lg}=filter_string_final{lg,:};
        sp_filter_idx{sp}=filter_idx_final;
        
        sp_data_string{sp,lg}=[legend_identifiers{lg}, ': ',data_string_final{lg,:}];
        sp_dim_names{sp,lg}=dim_names_final(lg,:);
        sp_statistics_on_idx(sp,lg)=statistics_on_idx_final(lg);
        sp_statistics_type{sp,lg}=statistics_type_final{lg};
    end
    if ~isfield(plot_mode,'legendstring')
        plot_mode.legendstring=[];
    end
    if isfield(plot_mode,'legendstring')  && ~isempty(plot_mode.legendstring) && size(plot_mode.legendstring,1)==1 && no_subplots>1
        plot_mode.legendstring=repmat(plot_mode.legendstring,[no_subplots 1 1 1]);
    end
    
    if isempty(plot_mode.legendstring) || size(plot_mode.legendstring(sp,:),2)<no_filterentries_final
        if isfield(plot_mode,'legendstring') && numel(plot_mode.legendstring)<no_filterentries_final
            sp_legend_string(sp,:)=legend_identifiers(1:no_filterentries_final)';
        elseif numel(plot_mode.legendstring)>=no_filterentries_final && any(size(plot_mode.legendstring))==1
            sp_legend_string(sp,:)=plot_mode.legendstring(1:no_filterentries_final);
        end
    else
        sp_legend_string(sp,:)=plot_mode.legendstring(sp,:);
    end
end

% Determine number of data points used in last statistical operation:
data_no_datapoints=cell(no_subplots,no_filterentries_final);
if ~isempty(spdata_sign)
    for sp=1:max(1,no_subplots)
        for lg=1:no_filterentries_final
            if ~strcmpi(statistics_type{pidx},'Hist')
                data_no_datapoints{sp,lg}=sum(~isnan(spdata_sign{sp,lg}),sp_statistics_on_idx(sp,lg));
            else
                data_no_datapoints{sp,lg}=sum(~isnan(spdata_signMedian{sp,lg}),sp_statistics_on_idx(sp,lg));
            end
        end
    end
end
% Gernerate subplot titles
subplotstring=cell(no_subplots,1);
for sp=1:no_subplots
    
    prts=dim_names(subplot_idx_final{sp,1});
    % Order field names according to dim_names:
    prt_orden=NaN(size(prts,1),1);
    for k1=1:size(prts,1)
        prt_orden(k1)=find(strcmp(dim_names,prts(k1)));
    end
    [~,prt_orden_idx]=sort(prt_orden);
    prts=prts(prt_orden_idx);
    for k=1:size(subplot_idx_final{sp,2},2)
        
        subplotstring{sp}=[subplotstring{sp} prts{k} ' ' num2str(subplot_idx_final{sp,2}(k)) ', '];
    end
    subplotstring{sp}=subplotstring{sp}(1:end-2);
    
end
% Determine dimension of xaxis/yaxis:
xaxis_length=cellfun(@(x) size(x,xaxis_idx{1}),spdata);
xaxis_length=max(xaxis_length(:));
if ~isempty(yaxis_idx) && ~isempty(yaxis_idx{1})
    yaxis_length=cellfun(@(x) size(x,yaxis_idx{1}),spdata);
    yaxis_length=max(yaxis_length(:));
else
    yaxis_length=0;
    yaxis_idx={[]};
end
disp('Done.')
%% Generate output structure
output.cell_size=cell_size_total;
output.dim_names_original=dim_names_orig;
output.dim_names_inplot=extra_dim_names;
output.data=spdata;
output.data_sign=spdata_sign;
output.data_Median=spdata_Median;
output.data_signMedian=spdata_signMedian;
output.data_dev=spdata_dev;
output.statistics_on_idx=sp_statistics_on_idx;
output.statistics_type=sp_statistics_type;
output.data_no_datapoints=data_no_datapoints;
output.data_string=sp_data_string;
output.filterstring=sp_filter_string;
output.legendstring=sp_legend_string;
output.dim_names=sp_dim_names;
output.subplotstring=subplotstring;
output.filter_idx=sp_filter_idx;
output.subplot_idx=subplot_idx_final;
output.data_idx=data_idx;
output.figure_idx=figure_idx;
output.xaxis_idx=xaxis_idx{1};
output.xaxis_length=xaxis_length;
output.yaxis_idx=yaxis_idx{1};
output.yaxis_length=yaxis_length;
output.display_mode=plot_mode.display_mode;
if isfield(plot_mode,'edges')
    output.dist_edges = plot_mode.edges;
end
if isfield(plot_mode,'dist_method')
    output.dist_method = plot_mode.dist_method;
end
if isfield(plot_mode,'dist_bandwidth')
    output.dist_bandwidth = plot_mode.dist_bandwidth;
end

%%% Prepare for plot [RCH 20170120]
output.xaxis = plot_mode.xaxis;
output.yaxis = plot_mode.yaxis;
if isfield(plot_mode,'label') && ~isempty(plot_mode.label)
    label = plot_mode.label;
else
    label = '';
end
% Xaxis
if (strcmpi(plot_mode.xaxis,'EdgesX')|| strcmpi(plot_mode.xaxis,'EdgeX')) && isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
    if iscell(plot_mode.edges) && numel(plot_mode.edges) == 1 && isnumeric(plot_mode.edges{1}) % distributions
        output.XTick{1} = plot_mode.edges{1};
        output.xlabel = label;
    elseif iscell(plot_mode.edges) && numel(plot_mode.edges) == 2 % 2d plot & maps
        output.XTick{1} = plot_mode.edges{1};
    else
        output.XTick = plot_mode.edges;
        output.xlabel = label;

    end
elseif strcmpi(plot_mode.xaxis,'Time') && isfield(plot_mode,'timeintervals_in_min') && ~isempty(plot_mode.timeintervals_in_min)
    output.XTick = (0:xaxis_length-1)*plot_mode.timeintervals_in_min/2;
    output.xlabel = 'Time (min)';
    output.ylabel = label;
elseif strcmpi(plot_mode.xaxis,'Time') && isfield(plot_mode,'time_intervals_min') && ~isempty(plot_mode.time_intervals_min)
    output.XTick = (0:xaxis_length-1)*plot_mode.time_intervals_min/2;
    output.xlabel = 'Time (min)';
    output.ylabel = label;

else % no edges, no histogram
    switch plot_mode.xaxis{1}
        case 'Group'
            output.XTick = 1:cell_size_total(1);
        case 'Subset'
            output.XTick = 1:cell_size_total(2);
        case 'Trial'
            output.XTick = 1:cell_size_total(3);
        case 'Focal'
            output.XTick = 1:cell_size_total(5);
        case 'Neighbor'
            output.XTick = 1:cell_size_total(6);
    end
    output.xlabel = plot_mode.xaxis{1};
    output.ylabel = label;
end
if isfield(output,'XTick') && ~isempty(output.XTick)
    output.XTickLabel = output.XTick;
end

% Yaxis
if (strcmpi(plot_mode.yaxis,'EdgesY')|| strcmpi(plot_mode.yaxis,'EdgeY')) && isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)
    if iscell(plot_mode.edges) && numel(plot_mode.edges) == 1 && isnumeric(plot_mode.edges{1}) % distributions
        output.YTick{1} = plot_mode.edges{2};
    elseif iscell(plot_mode.edges) && numel(plot_mode.edges) == 2 % 2d plot & maps
        output.YTick{1} = plot_mode.edges{2};
    else
        output.YTick = plot_mode.edges;
    end
elseif strcmpi(plot_mode.yaxis,'Time') && isfield(plot_mode,'timeintervals_in_min') && ~isempty(plot_mode.timeintervals_in_min)
    output.YTick = (0:yaxis_length-1)*plot_mode.timeintervals_in_min/2;
    output.ylabel = 'Time (min)';
elseif strcmpi(plot_mode.yaxis,'Time') && isfield(plot_mode,'time_intervals_min') && ~isempty(plot_mode.time_intervals_min)
    output.YTick = (0:xaxis_length-1)*plot_mode.time_intervals_min/2;
    output.ylabel = 'Time (min)';
elseif ~isempty(plot_mode.yaxis) % no edges, no histogram
    switch plot_mode.yaxis{1}
        case 'Group'
            output.YTick = 1:cell_size_total(1);
        case 'Subset'
            output.YTick = 1:cell_size_total(2);
        case 'Trial'
            output.YTick = 1:cell_size_total(3);
        case 'Focal'
            output.YTick = 1:cell_size_total(5);
        case 'Neighbor'
            output.YTick = 1:cell_size_total(6);
    end
    output.ylabel = plot_mode.yaxis{1};
else
    output.YTick = [];
end
if isfield(output,'YTick') && ~isempty(output.YTick)
    output.YTickLabel = output.YTick;
end


if isfield(plot_mode,'ylabel') && ~isempty(plot_mode.ylabel)
    output.ylabel = plot_mode.ylabel;
end
if isfield(plot_mode,'xlabel')&& ~isempty(plot_mode.xlabel)
    output.xlabel = plot_mode.xlabel;
end

if ((strcmpi(plot_mode.xaxis,'EdgesX')|| strcmpi(plot_mode.xaxis,'EdgeX')) && isfield(plot_mode,'edges') && ~isempty(plot_mode.edges)) && ...
        isfield(plot_mode,'normalization') && ~any(strcmpi(plot_mode.normalization,'none')) && ...
        isfield(plot_mode,'ylabel') && strcmpi(plot_mode.ylabel,'counts')
    plot_mode.ylabel = 'P';
end
        




end


function [plot_mode_idx, dim_name_act]=idSocial_function_wrapper_cellstr2idx(plot_mode,strcols,dim_names,dim_name_act)
% Substitutes strings in the strcol-th column of plot_mode by their index in dim_names.
if nargin<4 || isempty(dim_name_act)
    dim_name_act=repmat(dim_names',max(1,size(plot_mode,1)),1);
end

% if size(plot_mode,1) ~= size(dim_name_act,1)
%
% end

if size(plot_mode,1)==1 && size(dim_name_act,1)>1 % If combinations for filter has been calculated automatically
    plot_mode=repmat(plot_mode,size(dim_name_act,1),1);
end

if size(plot_mode,1)~=size(dim_name_act,1) && ~isempty(plot_mode)
    error([mfilename ': SizeMismatch'],[mfilename ': Size Mismatch']);
end
% if isempty(plot_mode)
%    	warning([mfilename ': Empty mode options'],[mfilename ':  Empty mode options']);
% end

if ~isempty(plot_mode) && ~isempty(strcols)
    plot_mode_idx=plot_mode;
    for strcol=strcols
        
        if  ~(strcol>size(plot_mode,2))
            % Find 'isempties'
            isemp=cellfun(@(x) isempty(x),plot_mode(:,strcol));
            plot_mode(isemp,strcol)=repmat({''},sum(isemp),1);
            % Substitute
            try
                LGcell = cellfun(@(x) strfind(plot_mode(:,strcol),x), dim_names,'UniformOutput',false);
            catch
                keyboard
            end
            LGcell=[LGcell{:}];
            LGarray=NaN(size(LGcell));
            for row=1:size(LGcell,1)
                for col=1:size(LGcell,2)
                    LGarray(row,col)=~isempty(LGcell{row,col});
                end
                
                idces=find(LGarray(row,:));
                if isempty(idces)
                    plot_mode_idx{row,strcol}=[];
                else
                    allready_gone=strcmp(dim_name_act(row,:),'');
                    if ~isempty(intersect(find(allready_gone),idces))
                        affected_dims=intersect(find(allready_gone),idces);
                        aff_names=[];
                        for k=1:size(affected_dims,1)
                            aff_names=[aff_names ', ' dim_names{affected_dims(k)}];
                        end
                        error([mfilename ': Field(s) ''' aff_names(3:end) ''' appear(s) twice!'])
                    end
                    dim_name_act(row,idces)=repmat({''},length(idces),1);
                    plot_mode_idx{row,strcol}=idces;
                end
            end
            
        end
        
    end
else
    plot_mode_idx=[];
end
end

function [plot_mode_idx, dim_names_act]=idSocial_function_wrapper_idxcell2combinations(plot_mode,cell_size,dim_names)
%

if size(plot_mode,2)<2 % No indices given, size(plot_mode)=1; calculate possible combinations
    no_lgs=1;
    if ~isempty(plot_mode{1})
        no_lgs_sep=NaN(size(plot_mode{1}));
    else
        no_lgs_sep=1;
    end
    for lg=1:size(plot_mode{1},2)
        no_lgs=no_lgs*cell_size(plot_mode{1}(lg)); % {1} because if size(filter_idx,1)==1, size(yval) also ==1.
        no_lgs_sep(lg)=cell_size(plot_mode{1}(lg));
    end
    if ~isempty(no_lgs_sep)
        sets=cell(size(no_lgs_sep,1),1);
        for k=1:size(no_lgs_sep,2)
            sets{k}=1:no_lgs_sep(k);
        end
        c = cell(1, numel(sets));
        [c{:}] = ndgrid( sets{end:-1:1} );
        lg_combi = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
        lg_combi = lg_combi(:,end:-1:1);
        lg_combi2 = cell(size(lg_combi,1),1);
        for lt=1:size(lg_combi,1)
            lg_combi2{lt}=lg_combi(lt,:);
        end
        filter_idx=repmat({plot_mode{1}},size(lg_combi2,1),1);
        filter_idx(:,2)=lg_combi2;
        plot_mode_idx=filter_idx;
        dim_names_act=repmat(dim_names,size(plot_mode_idx,1),1);
    end
    
else
    plot_mode_idx=plot_mode;
    dim_names_act=dim_names;
end
end
