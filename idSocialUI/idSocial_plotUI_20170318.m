function idSocial_plotUI(input_data,funcstring,act_set,figure_windows)

if nargin < 3 || isempty(act_set)
    act_set = 1; % Change this to last set!
end
if nargin < 4 || isempty(figure_windows)
    figure_windows = []; % Change this to last set!
end

sign_handle = [];




act_data = input_data.(funcstring).results(:,act_set+1);
pm = input_data.(funcstring).(['Set' num2str(act_set)]).plot_mode;
dm = pm.display_mode;
[~,project_name]=idSocial_recursiveGetOptionsFromOptionsCell(input_data.options,'project_name');

dat = act_data{3};
plotResultsIndivData = act_data{4};

% For distributions:
distr_statisics = pm.statistics{1};
distribution_flag = ~isempty(strfind(lower(distr_statisics),'hist'));
distr_statistics = strrep(lower(distr_statisics), 'hist+','');

labels = act_data(9);

xticks = act_data{8}{1};
yticks = act_data{8}{2};

xticklabel = act_data{9}{1};
yticklabel = act_data{9}{2};

if iscell(xticklabel) && numel(xticklabel)==1; xticklabel = xticklabel{1}; end
if iscell(xticks); xticks = xticks{1}; end

if iscell(yticklabel)&& numel(yticklabel)==1; yticklabel = yticklabel{1}; end
if iscell(yticks); yticks = yticks{1}; end


xaxislabel = act_data{7}{1};
yaxislabel = act_data{7}{2};

errH = [];

% For maps:
comp_set_idx = -1;
Xctrs = [];
Yctrs = [];

figure_windows = figure_windows(ishandle(figure_windows));

%% For statistics


if isempty(figure_windows)
    no_plots = 0;
    lineH = NaN(1,1);
else
    lineHT = findall(figure_windows(end),'Type','Line');
    check_idx = false(1,size(lineHT,1));
    for k=1:size(lineHT,1)
        udt = get(lineHT,'UserData');
        if isempty(udt) || ~isstruct(udt) || ~isfield(udt,'project_name')
            if ~isempty(lineHT)
                check_idx(k)=true;
            end
        end
    end
    lineHT(check_idx)=[];
    no_plots = numel(lineHT);
    if ~isempty(lineHT)
        lineH = NaN(1,size(lineHT,2)+1);
        lineH(1,1:size(lineHT,2)) = lineHT;
    end
end
linewidth = NaN(no_plots+1,1);
linewidth_DistrStats = NaN(no_plots+1,1);
linecolor = NaN(no_plots+1,3);
xdataCell = cell(no_plots+1,1);
ydataCell = cell(no_plots+1,1);
errH = NaN(no_plots+1,1);
distrStatsH = NaN(no_plots+1,1);
errorbar_type = repmat({'None'},[no_plots+1,1]);


for k = 1:no_plots
    linewidth(k) = get(lineH(k),'LineWidth');
    linewidth_DistrStats(k) = linewidth(k);
    linecolor(k,:) = get(lineH(k),'Color');
    set(lineH(k),'DeleteFcn',{@LineDelete,k});
    xdataCell{k} = get(lineH(k),'XData');
    ydataCell{k} = get(lineH(k),'YData');
    addlistener(lineH(k), 'LineWidth', 'PostSet', @(hAxes, eventData) setLineWidth(hAxes, eventData,k,lineH(k)));
    addlistener(lineH(k), 'Color', 'PostSet', @(hAxes, eventData) setColor(hAxes, eventData,k,lineH(k)));
    
end



marked_lines = zeros(1,no_plots+1);
marked_DistrStats = zeros(1,no_plots+1);

marked_points = zeros(max(cellfun(@(x) numel(x),xdataCell)),no_plots+1);
markh = NaN(max(cellfun(@(x) numel(x),xdataCell)),no_plots+1);

no_possiblemarks = inf;

%%


if  strcmpi(dm,'dynamicMapPolar')
    plotfh = figure('Units','normalized','Position',[.45 .3 .5 .5],'Name','Plot Results','NumberTitle','off','Color','w');
    figure_windows = [figure_windows plotfh];
    figure(figure_windows(end))
    if strcmpi(pm.statistics,'positive_ratio') || ...
            (strcmpi(pm.statistics,'pool') && strcmpi(pm.data,'positive_ratio'))
        limval = max(1-min(dat(:)), max(dat(:)));
        clims=[1-limval limval];
    else
        limval = max(abs(dat(:)));
        clims=[-limval  limval];
    end
    
    [h,ax_lines,Xctrs,Yctrs]= ...
        idSocial_imagescPolar(xticks,yticks,dat',clims,[12 12]);
    info.project_name = project_name;
    info.funcstring = funcstring;
    info.act_set = act_set;
    set(h,'UserData',info)
    
    set(gca,'YDir','normal')
    set(h, 'EdgeColor', 'none');
    axis equal  tight;% off
    hold on
    xlabel(xaxislabel)
    ylabel(yaxislabel)
    set(gca,'XTick',[-xticks(end:-1:2) xticks],'XTickLabel',[-xticks(end:-1:2) xticks])
    set(gca,'YTick',[-xticks(end:-1:2) xticks],'YTickLabel',[-xticks(end:-1:2) xticks])
    colormap4InteractionMapsRedBlue = [];
    load('colormap4InteractionMapsRedBlue.mat')
    colormap(colormap4InteractionMapsRedBlue)
    idSocial_improveGraphics(figure_windows(end))
    
    mh = uimenu(gcf,'Label','Statistics');
    eh2 = uimenu(mh,'Label','Statistical Signifcance');
    
    eh21 = uimenu(eh2,'Label','Apply test','Enable','off');
    
    test_list = idSocial_auxiliaries_statisticalSignificance;
    for tl = 1:numel(test_list)
        spec_options = idSocial_auxiliaries_statisticalSignificance([],test_list{tl});
        eh231 = uimenu(eh21,'Label',test_list{tl},'Callback',{@setStatTest,spec_options});
    end
    
    eh22 = uimenu(eh2,'Label','Select data');
    eh24 = uimenu(eh2,'Label','Remove significance markers','Callback',{@remSignMarkers,sign_handle});
    
    ticks_actSet = input_data.(funcstring).results{8,act_set+1};
    no_sets = size(input_data.(funcstring).results,2)-1;
    sign_sets = false(1,no_sets);
    for as = 1:no_sets
        if as~=act_set
            xt = input_data.(funcstring).results{8,as+1};
            sign_sets(as) = isequal(ticks_actSet{1},xt{1}) & isequal(ticks_actSet{2},xt{2});
        end
    end
    
    eh25 = NaN(1,no_sets);
    
    % If random controls exist:
    if size(input_data.(funcstring).results(:,act_set+1),1)>9 && ~isempty(input_data.(funcstring).results(10,act_set+1))
        eh25(no_sets+1) = uimenu(eh22,'Label','Control','Callback',{@setSetSelect,no_sets+1});
    end
    
    for as = 1:no_sets
        if sign_sets(as)
            eh25(as) = uimenu(eh22,'Label',['Set ' num2str(as)],'Callback',{@setSetSelect,as});
        end
    end
    
    eh24b = uimenu(eh22,'Label','Clear selection','Callback',@clearSetSelection);
    
    
elseif  strcmpi(dm,'MapPolar')
    plotfh = figure('Units','normalized','Position',[.45 .3 .5 .5],'Name','Plot Results','NumberTitle','off','Color','w');
    figure_windows = [figure_windows plotfh];
    figure(figure_windows(end))
    %     if strcmpi(pm.statistics,'positive_ratio')
    %         limval = max(1-min(dat(:)), max(dat(:)));
    %         clims=[1-limval limval];
    %     else
    limval = max(abs(dat(:)));
    clims=[0  limval];
    %     end
    
    [h,ax_lines,Xctrs,Yctrs]= ...
        idSocial_imagescPolar(xticks,yticks,dat',clims,[12 12]);
    
    info.project_name = project_name;
    info.funcstring = funcstring;
    info.act_set = act_set;
    set(h,'UserData',info)
    set(gca,'YDir','normal')
    
    set(h, 'EdgeColor', 'none');
    axis equal  tight;% off
    hold on
    xlabel(labels{1})
    ylabel(labels{2})
    set(gca,'XTick',[-xticks(end:-1:2) xticks],'XTickLabel',[-xticks(end:-1:2) xticks])
    set(gca,'YTick',[-xticks(end:-1:2) xticks],'YTickLabel',[-xticks(end:-1:2) xticks])
    %         load('colormap4InteractionMapsRedBlue.mat')
    colormap(jet)
    idSocial_improveGraphics(figure_windows(end))
    
else
    if isempty(figure_windows)
        plotfh = figure('Units','normalized','Position',[.45 .3 .5 .5],'Name','Plot Results','NumberTitle','off','Color','w');
        figure_windows = [figure_windows plotfh];
        
        
    else
        mh = findobj(gcf,'Label','Statistics');
        delete(mh);
    end
    figure(figure_windows(end))
    
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@datacursorupdatefcn)
    
    mh = uimenu(gcf,'Label','Statistics');
    eh1 = uimenu(mh,'Label','Show statistics');
    eh11 = uimenu(eh1,'Label','Standard deviation','Callback',@addSTD);
    eh12 = uimenu(eh1,'Label','SEM','Callback',@addSEM);
    if distribution_flag
        eh12b = uimenu(eh1,'Label',['Distribution ' distr_statistics],'Callback',@addDistrStats);
    end
    eh13 = uimenu(eh1,'Label','Remove statistics markers','Callback',@remError);
    eh2 = uimenu(mh,'Label','Statistical Signifcance');
    
    eh21 = uimenu(eh2,'Label','Apply test');
    eh22 = uimenu(eh2,'Label','Select data');
    eh22a = uimenu(eh22,'Label','Select lines','Callback',@setLineSelect);
    eh22b = uimenu(eh22,'Label','Select data points','Callback',@setPointSelect);
    if distribution_flag
        eh23b = uimenu(eh22,'Label',['Select distribution ' distr_statistics],'Callback',@setDistrStatsSelect,'Enable','off');
    end
    eh24 = uimenu(eh22,'Label','Clear selection','Callback',@clearSelection);
    eh24 = uimenu(eh2,'Label','Remove significance markers','Callback',{@remSignMarkers,sign_handle});
    
    test_list = idSocial_auxiliaries_statisticalSignificance;
    
    for tl = 1:numel(test_list)
        spec_options = idSocial_auxiliaries_statisticalSignificance([],test_list{tl});
        eh231 = uimenu(eh21,'Label',test_list{tl},'Callback',{@setStatTest,spec_options});
    end
    
    hold on
    try
        ph=plot(xticks,dat,'-o');
        info.project_name = project_name;
        info.funcstring = funcstring;
        info.act_set = act_set;
        set(ph,'UserData',info)
        
        
        linewidth(end) =  get(ph,'LineWidth');
        linewidth_DistrStats(end) =  get(ph,'LineWidth');
        linecolor(end,:) = get(ph,'Color');
        set(ph,'DeleteFcn',{@LineDelete,no_plots+1});
        
        
        xdataCell{end} = get(ph,'XData');
        ydataCell{end} = get(ph,'YData');
        lineH(end) = ph;
        addlistener(ph, 'LineWidth', 'PostSet', @(hAxes, eventData) setLineWidth(hAxes, eventData,no_plots+1,lineH(no_plots+1)));
        addlistener(ph, 'Color', 'PostSet', @(hAxes, eventData) setColor(hAxes, eventData,no_plots+1,lineH(no_plots+1)));
        %         addlistener(gca, 'YLim', 'PostSet', @(hAxes, eventData) afterYLimChange(hAxes, eventData));
        
        errorbar_type{end} = 'None';
        
        
        
        udat = get(gcf,'UserData');
        udat = [udat act_set];
        set(gcf,'UserData',udat);
        xlabel(xaxislabel)
        ylabel(yaxislabel)
        set(gca,'XTick',xticks,'XTickLabel',xticklabel)
        idSocial_improveGraphics(figure_windows(end))
        ylim = get(gca,'YLim');
    catch
        close(gcf)
        warndlg(['Sorry, I cannot plot the data.'])
        
    end
end

    function setStatTest(src,event,act_options)
        
        out_opts = idSocialUI_propertyTable([],act_options);
        
        
        if ~isempty(out_opts)
            
            if strcmpi(dm,'dynamicMapPolar')
            else
                act_info = get(lineH,'UserData');

                if any(marked_points(:))
                    [~,point_idces] = sort(marked_points(:),'descend');
                    point_idces = point_idces(1:sum(marked_points(:)>0));
                    lidx = ceil(point_idces/size(marked_points,1));
                    pidx = mod(point_idces-1,size(marked_points,1))+1;
                elseif any(marked_DistrStats(:))
                    [~,point_idces] = sort(marked_DistrStats(:),'descend');
                    point_idces = point_idces(1:sum(marked_DistrStats(:)>0));
                    lidx = ceil(point_idces/size(marked_DistrStats,1));
                    
                elseif any(marked_lines)
                    [~, line_idces] = sort(marked_lines,'descend');
                else
                    line_idces = 1:size(act_info,1);
                    
                end
            end
            
            if strcmpi(dm,'dynamicMapPolar')
                display_mode = dm;
                val1 = input_data.(funcstring).results{4,act_set+1};
                if comp_set_idx>0 && comp_set_idx<no_sets+1
                    val2 = input_data.(funcstring).results{4,comp_set_idx+1};
                elseif comp_set_idx==no_sets+1
                    val2 = input_data.(funcstring).results{11,act_set+1};
                end
                val = NaN(2,max(size(val1,1),size(val1,2)),size(val1,2),size(val1,3));
                val(1,1:size(val1,1),:,:) =val1;
                val(2,1:size(val2,1),:,:) =val2;
                lineHReorder = [];
                
            else
                display_mode = 'plot2d';
                if any(marked_points(:)) || any(marked_DistrStats(:))
                    lineHReorder = lineH;
                    datasetIndiv = cell(size(lidx,1),1);
                    xticks_act = cell(size(lidx,1),1);
                    count = 1;
                    for lh = lidx'
                        
                        if any(marked_points(:))
                            datasetIndiv{count} = input_data.(act_info{lh}.funcstring).results{4,act_info{lh}.act_set+1}(:,pidx(count));
                            xtdum = get(lineH(lh),'XData');
                            xticks_act{count} = xtdum(pidx(count));
                        elseif any(marked_DistrStats(:))
                            datasetIndiv{count} = input_data.(act_info{lh}.funcstring).results{6,act_info{lh}.act_set+1};
                            xtdum = get(distrStatsH(lh),'XData');
                            xticks_act{count} = xtdum(1);
                        end
                        
                        
                        count = count + 1;
                    end
                    no_max_dpts = max(cellfun(@(x) size(x,1),datasetIndiv));
                    val = NaN(size(lidx,1),no_max_dpts,1);
                    for lh = 1:size(lidx,1)
                        
                        val(lh,1:size(datasetIndiv{lh},1),1) = datasetIndiv{lh};
                    end
                    display_mode = [xticks_act{:}];
                    
                else% any(marked_lines)
                    datasetIndiv = cell(size(act_info,1),1);
                    xticks_act = cell(size(act_info,1),1);
                    lineHReorder = lineH(line_idces); %Re-order for idSocial_auxiliaries_plotSignificance
                    for lh = line_idces
                        
                        datasetIndiv{lh} = input_data.(act_info{lh}.funcstring).results{4,act_info{lh}.act_set+1};
                        xticks_act{lh} = get(lineH(lh),'XData');
                        
                    end
                    min_x = min(cellfun(@(x) min(x),xticks_act));
                    max_x = max(cellfun(@(x) max(x),xticks_act));
                    no_max_dpts = max(cellfun(@(x) size(x,1),datasetIndiv));
                    no_parts = max_x - min_x + 1;
                    val = NaN(numel(line_idces),no_max_dpts,no_parts);
                    count = 1;
                    for lh = line_idces
                        start_idx = find(xticks_act{lh}(1)==(min_x:max_x));
                        val(count,1:size(datasetIndiv{lh},1),start_idx:start_idx+no_parts-1) = datasetIndiv{lh};
                        count =  count + 1;
                    end
                end
            end
            
            if ~isempty(out_opts)
                if any(any(marked_lines)) || any(any(marked_points)) || any(any(marked_DistrStats))
                    out_opts.one_way_only = true;
                else
                    out_opts.one_way_only = false;
                end
                
                [H, p, CI, stat, sampstat]= ...
                    idSocial_auxiliaries_statisticalSignificance(val,get(src,'Label'),out_opts);
                if iscell(out_opts.tail(end)) && iscell(out_opts.tail{end})
                    tail =  out_opts.tail{end}{1};
                else
                    tail =  out_opts.tail{1};
                end
                sign_handle_temp = idSocial_auxiliaries_plotSignificance(lineHReorder,p,out_opts.alpha,tail,display_mode,get(src,'Label'),Xctrs,Yctrs);
                sign_handle = [sign_handle sign_handle_temp];
                
                set(eh24,'Callback',{@remSignMarkers,sign_handle});
                if ~strcmpi(dm,'dynamicMapPolar')
                    for ml=1:numel(lineH)
                        set(lineH(ml), 'LineWidth', linewidth(ml));
                    end
                end
                clearSelection
                
            end
        end
    end

    function setLineSelect(src,event)
        plotedit off
        datacursormode off
        %         % Remove points
        %         for mp = 1:size(marked_points,1)
        %             for mp2 = 1:size(marked_points,2)
        %                 if ishandle(markh(mp,mp2))
        %                     delete(markh(mp,mp2));
        %                     marked_points(mp,mp2) = 0;
        %                 end
        %             end
        %         end
        clearSelection
        
        for k2 = 1:no_plots+1
            
            set(lineH(k2),'ButtonDownFcn',{@LineSelected,k2});
            
        end
        
    end


    function setDistrStatsSelect(src,event)
        plotedit off
        datacursormode off
        
        clearSelection
        
        for k2 = 1:no_plots+1
            
            set(distrStatsH(k2),'ButtonDownFcn',{@DistrStatsSelected,k2});
            
        end
        
    end

    function setPointSelect(~,event)
        plotedit off
        datacursormode off
        % Remove line selection
        %         marked_lines = zeros(size(marked_lines));
        %
        %         for ml=find(~marked_lines)
        %             set(lineH(ml), 'LineWidth', linewidth(ml));
        %         end
        
        clearSelection
        
        for k2 = 1:no_plots+1
            
            set(lineH(k2),'ButtonDownFcn',{@PointSelected,k2});
            
        end
        
    end

    function clearSelection(src,event)
        
        for mp = 1:size(marked_points,1)
            for mp2 = 1:size(marked_points,2)
                if ishandle(markh(mp,mp2))
                    delete(markh(mp,mp2));
                    marked_points(mp,mp2) = 0;
                end
            end
        end
        
        marked_lines = zeros(size(marked_lines));
        
        for ml=find(~marked_lines)
            if ishandle(lineH(ml)) && ishghandle(lineH(ml)) % && isgraphics(lineH(ml))
                set(lineH(ml), 'LineWidth', linewidth(ml));
            end
        end
        marked_DistrStats = zeros(size(marked_DistrStats));
        for ml=find(~marked_DistrStats)
            if ishandle(distrStatsH(ml)) && ishghandle(distrStatsH(ml)) %%&& isgraphics(distrStatsH(ml))
                set(distrStatsH(ml), 'LineWidth', linewidth_DistrStats(ml));
            end
        end
        
    end

    function clearSetSelection(src,event)
        for as1 = 1:size(eh25,2)
            if ~isnan(eh25(as1)) && as1<=no_sets
                set(eh25(as1),'Label',['Set ' num2str(as1)]);
            elseif ~isnan(eh25(as1)) && as1>no_sets
                set(eh25(as1),'Label','Control');
            end
        end
        set(eh21,'Enable','off')
    end

    function setSetSelect(src,event,comp_set)
        comp_set_idx = comp_set;
        set(eh21,'Enable','on')
        if comp_set<=no_sets
            set(eh25(comp_set),'Label',['Set ' num2str(comp_set) ' ' char(hex2dec('2713'))]);
            
        else
            set(eh25(comp_set),'Label',['Control' char(hex2dec('2713'))]);
        end
        
        for as1 = 1:size(eh25,2)
            if as1~=comp_set && ~isnan(eh25(as1)) && as1<=no_sets
                set(eh25(as1),'Label',['Set ' num2str(as1)]);
            elseif as1~=comp_set && as1>no_sets && ~isnan(eh25(as1))
                set(eh25(as1),'Label','Control');
            end
        end
        
    end

    function addSTD(src,event)
        
        act_info = get(lineH,'UserData');
        for lh = 1:size(act_info,1)
            if ~isempty(errH(lh)) && ishandle(errH(lh))
                delete(errH(lh));
                errH(lh) = NaN;
            end
            if ~strcmpi(errorbar_type{lh},'std');
                dataset = input_data.(act_info{lh}.funcstring).results{3,act_info{lh}.act_set+1};
                datasetIndiv = input_data.(act_info{lh}.funcstring).results{4,act_info{lh}.act_set+1};
                xticks_act = input_data.(act_info{lh}.funcstring).results{8,act_info{lh}.act_set+1};
                if iscell(xticks_act); xticks_act = xticks_act{1}; end
                errH(lh) = errorbar(xticks_act,dataset ,nanstd(datasetIndiv),'Color',get(lineH(lh),'Color'));
                uistack(errH(lh),'bottom')
                errorbar_type{lh} = 'std';
            end
            if isnan(errH(lh))
                errorbar_type{lh} = 'None';
            end
            
        end
    end

    function addSEM(src,event)
        act_info = get(lineH,'UserData');
        for lh = 1:size(act_info,1)
            if ~isempty(errH(lh)) && ishandle(errH(lh))
                delete(errH(lh));
                errH(lh) = NaN;
                
            end
            if ~strcmpi(errorbar_type{lh},'sem');
                dataset = input_data.(act_info{lh}.funcstring).results{3,act_info{lh}.act_set+1};
                datasetIndiv = input_data.(act_info{lh}.funcstring).results{4,act_info{lh}.act_set+1};
                xticks_act = input_data.(act_info{lh}.funcstring).results{8,act_info{lh}.act_set+1};
                if iscell(xticks_act); xticks_act = xticks_act{1}; end
                errH(lh) = errorbar(xticks_act,dataset ,nanstd(datasetIndiv)./sqrt(sum(~isnan(datasetIndiv))),'Color',get(lineH(lh),'Color'));
                uistack(errH(lh),'bottom')
                errorbar_type{lh} = 'sem';
            end
            if isnan(errH(lh))
                errorbar_type{lh} = 'None';
            end
            
        end
    end

    function addDistrStats(src,event)
        
        
        act_info = get(lineH,'UserData');
        for lh = 1:size(act_info,1)
            
            dataset = input_data.(act_info{lh}.funcstring).results{5,act_info{lh}.act_set+1};
            %                 datasetIndiv = input_data.(act_info{lh}.funcstring).results{6,act_info{lh}.act_set+1};
            distrStatsH(lh) = plot([dataset dataset], ylim,'Color',get(lineH(lh),'Color'));
            set(distrStatsH(lh),'DeleteFcn',{@distrStatDelete,no_plots+1})
            uistack(distrStatsH(lh),'bottom')
            
        end
        set(eh23b,'Enable','on');
    end


    function remError(src,event)
        act_info = get(lineH,'UserData');
        for lh = 1:size(act_info,1)
            if ~isempty(errH(lh)) && ishandle(errH(lh))
                delete(errH(lh));
                errH(lh) = NaN;
                
            end
            if ~isempty(distrStatsH) && ishandle(distrStatsH(lh)) && ishghandle(distrStatsH(lh)) %%&& isgraphics(distrStatsH(lh))
                delete(distrStatsH(lh))
                distrStatsH(lh)=NaN;
                set(eh23b,'Enable','off');
            end
            if isnan(errH(lh))
                errorbar_type{lh} = 'None';
            end
            
        end
    end


    function addStatSign(src,event)
        errorbar(xticks,dat,nanstd(plotResultsIndivData)./sqrt(sum(~isnan(plotResultsIndivData))),'Color',get(ph,'Color'))
    end

    function LineSelected(ObjectH, EventData,act_plot)
        
        
        
        if marked_lines(act_plot)>0
            %             for hc = 1:numel(H)
            %                 set(lineH(hc), 'LineWidth', linewidtlineH(hc));
            %             end
            set(lineH(act_plot), 'LineWidth', linewidth(act_plot));
            marked_lines(act_plot) = 0;
        else
            marked_lines(marked_lines>0) = marked_lines(marked_lines>0) + 1;
            marked_lines(act_plot) = 1;
            marked_lines(marked_lines>no_possiblemarks) = 0;
            for ml=find(~marked_lines)
                set(lineH(ml), 'LineWidth', linewidth(ml));
            end
            for ml=find(marked_lines)
                set(lineH(ml), 'LineWidth', linewidth(ml)+2);
            end
        end
    end

    function DistrStatsSelected(ObjectH, EventData,act_plot)
        
        
        
        if marked_DistrStats(act_plot)>0
            %             for hc = 1:numel(H)
            %                 set(lineH(hc), 'LineWidth', linewidtlineH(hc));
            %             end
            set(distrStatsH(act_plot), 'LineWidth', linewidth_DistrStats(act_plot));
            marked_DistrStats(act_plot) = 0;
            
        else
            marked_DistrStats(marked_DistrStats>0) = marked_DistrStats(marked_DistrStats>0) + 1;
            marked_DistrStats(act_plot) = 1;
            marked_DistrStats(marked_DistrStats>no_possiblemarks) = 0;
            for ml=find(~marked_DistrStats)
                set(distrStatsH(ml), 'LineWidth', linewidth_DistrStats(ml));
            end
            for ml=find(marked_DistrStats)
                set(distrStatsH(ml), 'LineWidth', linewidth_DistrStats(ml)+2);
            end
        end
    end



    function PointSelected(ObjectH, EventData,act_plot)
        
        C = get (get(ObjectH,'Parent'), 'CurrentPoint');
        xlim = get(get(ObjectH,'Parent'),'XLim');
        ylim = get(get(ObjectH,'Parent'),'YLim');
        Cnorm = [(C(1,1) - xlim(1))/(xlim(2)-xlim(1)) (C(1,2) - ylim(1))/(ylim(2)-ylim(1))];
        xnorm = (xdataCell{act_plot} - xlim(1))/(xlim(2)-xlim(1));
        ynorm = (ydataCell{act_plot} - ylim(1))/(ylim(2)-ylim(1));
        dst = sqrt((Cnorm(1)-xnorm).^2+(Cnorm(2)-ynorm).^2);
        
        %         dst = sqrt((C(1,1)-xdataCell{act_plot}).^2+(C(1,2)-ydataCell{act_plot}).^2);
        [~, idx] = min(dst);
        
        if marked_points(idx,act_plot)>0
            delete(markh(idx,act_plot));
            marked_points(idx,act_plot) = 0;
        else
            
            marked_points(marked_points>0) = marked_points(marked_points>0) + 1;
            
            marked_points(idx,act_plot) = 1;
            for mp = 1:size(marked_points,1)
                for mp2 = 1:size(marked_points,2)
                    if ishandle(markh(mp,mp2)) && marked_points(mp,mp2)>no_possiblemarks
                        delete(markh(mp,mp2));
                        marked_points(mp,mp2) = 0;
                    end
                end
            end
            
            markh(idx,act_plot) = plot(xdataCell{act_plot}(idx),ydataCell{act_plot}(idx),'s','Color',linecolor(act_plot,:),'MarkerSize',8); %get(ObjectH,'MarkerSize')+4
            uistack(markh(idx,act_plot),'bottom')
        end
    end

    function LineDelete(ObjectH, EventData,act_plot)
        
        xd = get(ObjectH,'XData');
        yd = get(ObjectH,'YData');
        
        xd2 = NaN; yd2 = NaN;
        if ~isempty(gco) && isprop(gco,'XData') && isprop(gco,'YData')
            xd2 = get(gco,'XData');
            yd2 = get(gco,'YData');
        end
        
        if isequal(xd,xd2) && isequal(yd,yd2)
            set(gco,'UserData',get(ObjectH,'UserData'));
            set(gco,'DeleteFcn',{@LineDelete,act_plot});
            lineH(act_plot) = gco;
            addlistener(lineH(act_plot), 'LineWidth', 'PostSet', @(hAxes, eventData) setLineWidth(hAxes, eventData,act_plot,lineH(act_plot)));
            addlistener(lineH(act_plot), 'Color', 'PostSet', @(hAxes, eventData) setColor(hAxes, eventData,act_plot,lineH(act_plot)));
            if isprop(gco,'Color') && isprop(ObjectH,'Color')
                set(gco,'Color',get(ObjectH,'Color'));
            elseif isprop(gco,'FaceColor') && isprop(ObjectH,'Color')
                set(gco,'FaceColor',get(ObjectH,'Color'));
            end
        else
            
            line_ud = get(ObjectH,'UserData');
            
            
            if ~isempty(line_ud) && isstruct(line_ud) && isfield(line_ud,'significanceMarkerH')
                delete(line_ud.significanceMarkerH)
            end
            if ~isempty(distrStatsH) && ishandle(distrStatsH(act_plot))
                delete(distrStatsH(act_plot))
                distrStatsH(act_plot)=NaN;
            end
            if any(marked_points(:,act_plot)>0)
                delete(markh(ishandle(markh(:,act_plot)),act_plot));
                %             marked_points(:,act_plot) = [];
            end
        end
        
    end

    function  distrStatDelete(obj,event_obj,act_plot)
        if ~any(ishandle(distrStatsH) & ishghandle(distrStatsH))%% & isgraphics(distrStatsH))
            set(eh23b,'Enable','off');
        end
    end

    function remSignMarkers(obj,event_obj,sign_handle)
        delete(sign_handle);
        ud = get(gca,'Userdata');
        if ~isempty(ud) && isstruct(ud)
            ud.ypos_signmarker = [];
            set(gca,'Userdata',ud);
        else
            set(gca,'Userdata',[]);
        end
        ylim1 = get(gca,'YLim');
        set(gca, 'ylimmode', 'auto')
        %         ylim2 = get(gca,'YLim');
        %         set(gca,'YLim',[ylim1(1) ylim2(2)]);
        %         ylim1 = get(gca,'YLim');
        %         set(gca,'YLim',[ylim1(1) ylim1(2)+(diff(ylim1)*0.05)]);
    end

    function txt = datacursorupdatefcn(empt,event_obj)
        % Customizes text of data tips
        
        pos = get(event_obj,'Position');
        pval = get(get(event_obj,'Target'),'UserData');
        if ~isempty(pval) && isstruct(pval) && isfield(pval,'p')
            if isfield(pval,'test_type') && ~isempty(pval.test_type)
                txt = {pval.test_type;['p = ',num2str(pval.p(pos(1)==xticks))]};
            else
                txt = ['p = ',num2str(pval.p(pos(1)==xticks))];
            end
        else
            txt = {['X: ',num2str(pos(1))],...
                ['Y: ',num2str(pos(2))]};
            
        end
    end

    function setLineWidth(hAxes, eventData,act_ln,lineHandle)
        if plotedit(gcf,'isactive')
            linewidth(act_ln) = get(lineHandle,'LineWidth');
        end
        
    end

    function setColor(hAxes, eventData,act_ln,lineHandle)
        ud = get(lineHandle,'UserData');
        if ~isempty(ud) && isfield(ud,'significanceMarkerH')
            c1 = get(ud.significanceMarkerH,'MarkerEdgeColor');
            c2 = get(ud.significanceMarkerH,'MarkerFaceColor');
            if isequal(linecolor(act_ln,:),c1)
                set(ud.significanceMarkerH,'MarkerEdgeColor',get(lineHandle,'Color'))
            elseif isequal(linecolor(act_ln,:),c2)
                set(ud.significanceMarkerH,'MarkerFaceColor',get(lineHandle,'Color'))
            end
        end
        if ~isempty(distrStatsH) && ishandle(distrStatsH(act_ln))
            set(distrStatsH(act_ln),'Color',get(lineHandle,'Color'))
        end
        linecolor(act_ln,:) = get(lineHandle,'Color');
    end

%     function afterYLimChange(hAxes, eventData)
%         act_lims = get(gca,'YLim');
%         act_pos = NaN(numel(sign_handle),2);
%         for sh = 1:numel(sign_handle)
%             act_pos(sh,:) = get(sign_handle,'Position');
%         end
%         max_pos = max(act_pos(:,2));
%     end

end
